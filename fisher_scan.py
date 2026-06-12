"""Self-contained Fisher's-exact 3D neighborhood enrichment test.

No SIR_REPO dependency — only BioPython + scipy + numpy. This mirrors the
PyMOL-pipeline method used in the ATP2B2 / HMGCR burden work:

  * residue-residue distance = minimum all-atom distance, inflated to 1000 A
    when *both* directions of the PAE exceed the cutoff (i.e. the relative
    geometry of the two residues is poorly predicted),
  * residues with mean pLDDT <= cutoff are dropped,
  * one-sided Fisher's exact test per neighborhood center (case enrichment),
  * per-residue score = -log10(min p over every neighborhood containing it).

The center residue is excluded from its own neighborhood (diagonal = inf),
matching the reference implementation.

Returns a per-residue DataFrame plus a small metadata dict; the caller colors
the structure by `neglog10_min_p` (white -> red) and renders the table.
"""
import gzip
import json
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import fisher_exact
from Bio.PDB import PDBParser

PLDDT_CUTOFF_DEFAULT = 50.0


# ---------------------------------------------------------------------------
# PAE (AlphaFold DB JSON) -> n x n numpy array
# ---------------------------------------------------------------------------
def load_pae(pae_path):
    """Parse an AlphaFold-DB PAE JSON file into an (n, n) float array.

    Handles the v4 style ([{"predicted_aligned_error": [[...]]}]) and the older
    sparse style ({"residue1": [...], "residue2": [...], "distance": [...]}).
    """
    with open(pae_path) as fh:
        data = json.load(fh)
    if isinstance(data, list):
        data = data[0] if data else {}
    if isinstance(data, dict):
        for key in ("predicted_aligned_error", "pae", "pae_matrix"):
            if key in data:
                return np.asarray(data[key], dtype=float)
        if {"residue1", "residue2", "distance"} <= set(data):
            r1 = np.asarray(data["residue1"], dtype=int)
            r2 = np.asarray(data["residue2"], dtype=int)
            dist = np.asarray(data["distance"], dtype=float)
            n = int(max(r1.max(), r2.max()))
            mat = np.zeros((n, n), dtype=float)
            mat[r1 - 1, r2 - 1] = dist
            return mat
    raise ValueError("Unrecognized PAE JSON format.")


# ---------------------------------------------------------------------------
# Structure parsing
# ---------------------------------------------------------------------------
def _parse_structure(pdb_gz_path):
    """Return (resnums, plddt, atom_coords, atom_res) for the first model.

    resnums   : (n,) residue numbers in chain order
    plddt     : (n,) mean per-residue B-factor (AlphaFold pLDDT)
    atom_coords: (A, 3) all atom coordinates, grouped by residue in order
    atom_res  : (A,) residue index (0..n-1) for each atom (sorted/contiguous)
    """
    parser = PDBParser(QUIET=True)
    opener = gzip.open if str(pdb_gz_path).endswith(".gz") else open
    with opener(pdb_gz_path, "rt") as fh:
        structure = parser.get_structure("p", fh)

    resnums, plddt, atom_coords, atom_res = [], [], [], []
    idx = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":          # skip hetero atoms / waters
                    continue
                coords = [a.get_coord() for a in residue]
                if not coords:
                    continue
                resnums.append(int(residue.id[1]))
                plddt.append(float(np.mean([a.bfactor for a in residue])))
                for c in coords:
                    atom_coords.append(c)
                    atom_res.append(idx)
                idx += 1
        break                                      # first model only
    return (np.asarray(resnums, dtype=int),
            np.asarray(plddt, dtype=float),
            np.asarray(atom_coords, dtype=float),
            np.asarray(atom_res, dtype=int))


def _residue_distance_matrix(atom_coords, atom_res, n):
    """Min all-atom distance between every pair of residues. Diagonal = inf."""
    # atom_res is contiguous and sorted by construction -> segment starts:
    starts = np.searchsorted(atom_res, np.arange(n))
    d = np.full((n, n), np.inf, dtype=float)
    for i in range(n):
        ai = atom_coords[starts[i]: (starts[i + 1] if i + 1 < n else len(atom_res))]
        amin = cdist(ai, atom_coords).min(axis=0)        # nearest atom of i to each atom
        d[i, :] = np.minimum.reduceat(amin, starts)      # reduce to per-residue min
    np.fill_diagonal(d, np.inf)                          # exclude self from own nbhd
    return d


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def run_fisher_3dnt(df_counts, pdb_gz_path, pae_path,
                    radius=15.0, pae_cutoff=15.0,
                    plddt_cutoff=PLDDT_CUTOFF_DEFAULT):
    """Run the Fisher 3D neighborhood test for one protein.

    df_counts : columns aa_pos, ac_case, ac_control (one row per variant residue).
    Returns (scores_df, meta).
      scores_df columns: aa_pos, center_p, min_p_containing, n_containing,
                         nbhd_case, nbhd_control, neglog10_min_p
                         (one row per pLDDT-valid residue).
      meta: n_tested, bonferroni_p, neglog10_bonferroni, n_case_total,
            n_ctrl_total, vmax.
    """
    resnums, plddt, atom_coords, atom_res = _parse_structure(pdb_gz_path)
    n = len(resnums)
    if n == 0:
        raise ValueError("No protein residues parsed from the structure.")

    d = _residue_distance_matrix(atom_coords, atom_res, n)

    # PAE inflation (skip when cutoff == 0, i.e. "off")
    if pae_cutoff and pae_cutoff > 0:
        pae = load_pae(pae_path)
        if pae.shape == (n, n):
            bad = (pae > pae_cutoff) & (pae.T > pae_cutoff)
            d[bad] = 1000.0
        # if PAE shape disagrees with the model, skip inflation rather than crash

    resnum_to_idx = {int(r): i for i, r in enumerate(resnums)}
    valid_mask = plddt > plddt_cutoff
    valid_resnums = {int(resnums[i]) for i in range(n) if valid_mask[i]}

    # per-residue case/control counts, restricted to valid residues in structure
    counts = {}
    for row in df_counts.itertuples(index=False):
        pos = int(row.aa_pos)
        if pos in valid_resnums:
            counts[pos] = (int(row.ac_case), int(row.ac_control))
    if not counts:
        raise ValueError(
            "No variant residues pass the pLDDT filter or map onto the structure. "
            "Check numbering/isoform, or lower the pLDDT cutoff.")

    n_case_total = sum(c for c, _ in counts.values())
    n_ctrl_total = sum(k for _, k in counts.values())
    if n_case_total == 0 or n_ctrl_total == 0:
        raise ValueError("Need at least one case and one control allele after filtering.")

    var_idx = np.array([resnum_to_idx[p] for p in counts], dtype=int)
    var_case = np.array([counts[p][0] for p in counts], dtype=float)
    var_ctrl = np.array([counts[p][1] for p in counts], dtype=float)

    tested_centers = sorted(counts.keys())          # resnums carrying a variant
    center_p, center_ab = {}, {}
    for r in tested_centers:
        ci = resnum_to_idx[r]
        within = d[ci, var_idx] <= radius
        a = int(var_case[within].sum())
        b = int(var_ctrl[within].sum())
        c = n_case_total - a
        dd = n_ctrl_total - b
        _, p = fisher_exact([[a, b], [c, dd]], alternative="greater")
        center_p[r] = float(p)
        center_ab[r] = (a, b)

    tested_idx = np.array([resnum_to_idx[r] for r in tested_centers], dtype=int)
    tested_p = np.array([center_p[r] for r in tested_centers], dtype=float)

    rows = []
    for i in range(n):
        if not valid_mask[i]:
            continue
        rn = int(resnums[i])
        dists = d[i, tested_idx]
        m = dists <= radius
        if m.any():
            min_p = float(tested_p[m].min())
            n_cont = int(m.sum())
        else:
            min_p, n_cont = np.nan, 0
        a, b = center_ab.get(rn, (np.nan, np.nan))
        rows.append({
            "aa_pos": rn,
            "center_p": center_p.get(rn, np.nan),
            "min_p_containing": min_p,
            "n_containing": n_cont,
            "nbhd_case": a,
            "nbhd_control": b,
        })
    scores = pd.DataFrame(rows)

    floor = 1e-300
    scores["neglog10_min_p"] = -np.log10(scores["min_p_containing"].clip(lower=floor))

    n_tested = len(tested_centers)
    bonf = 0.05 / n_tested if n_tested else np.nan
    vmax = float(np.nanmax(scores["neglog10_min_p"].values)) if len(scores) else 0.0
    n_kept = int(valid_mask.sum())
    _cp = np.array([center_p[r] for r in tested_centers]) if n_tested else np.array([])
    n_sig_bonf = int((_cp < bonf).sum()) if n_tested else 0
    n_sig_p01 = int((_cp < 0.01).sum()) if n_tested else 0
    min_p = float(_cp.min()) if n_tested else float("nan")
    meta = {
        "n_tested": n_tested,                       # variant residues tested (centers)
        "min_p": min_p,                             # smallest neighborhood p-value
        "n_sig_bonferroni": n_sig_bonf,             # centers below Bonferroni
        "n_sig_p01": n_sig_p01,                     # centers below p < 0.01
        "bonferroni_p": bonf,
        "neglog10_bonferroni": (-np.log10(bonf) if n_tested else np.nan),
        "n_case_total": n_case_total,
        "n_ctrl_total": n_ctrl_total,
        "vmax": vmax,
        "n_residues_total": n,
        "n_residues_kept": n_kept,
        "n_residues_omitted": n - n_kept,
        "frac_omitted": ((n - n_kept) / n) if n else 0.0,
        "mean_plddt": float(plddt.mean()) if n else float("nan"),
        "mean_plddt_kept": float(plddt[valid_mask].mean()) if n_kept else float("nan"),
    }
    return scores, meta
