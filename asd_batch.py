"""Batch-run the web Fisher 3D neighborhood test over the ASD proteins to build
the genome-wide Manhattan summary the site renders.

For each protein it fetches the AlphaFold structure (so it needs outbound network
— run it on your machine or a server, not an offline sandbox), maps that
protein's case/control variants, and records its top neighborhood. It then adds
two genome-wide, gene-level corrections:

  * BH-FDR across the per-protein top p-values, and
  * a pooled-permutation FWER: labels are permuted within every protein each
    replicate, and the global minimum p across all proteins is the null
    max-statistic (the genome-scale analogue of the within-protein permutation
    correction, which handles the correlated neighborhoods).

Output: asd_data/asd_web_results.tsv  (one row per protein).

Usage:
    python asd_batch.py                       # all proteins, 1000 sims
    python asd_batch.py --sims 0              # skip FWER (FDR only, fast)
    python asd_batch.py --limit 20            # quick test on 20 proteins
"""
import os
import argparse
import numpy as np
import pandas as pd

import sir_fetch as F
import fisher_scan as FS

HERE = os.path.dirname(os.path.abspath(__file__))
DEF_VARIANTS = os.path.join(HERE, "asd_data", "ASD_variants.tsv")
DEF_OUT = os.path.join(HERE, "asd_data", "asd_web_results.tsv")


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def pooled_fwer(models, obs_minp, n_sims, seed=0):
    """Gene-level FWER via a pooled min-P permutation across all proteins.

    models: list of dicts from FS.build_neighborhood_model (with 'pos_of_allele'
            precomputed). obs_minp: observed min-p per protein (same order).
    Returns (fwer_adj per protein, neglog10 threshold at FWER 0.05).
    """
    from scipy.stats import hypergeom
    rng = np.random.default_rng(seed)
    global_min = np.empty(int(n_sims), dtype=float)
    for s in range(int(n_sims)):
        gm = 1.0
        for m in models:
            ca = rng.choice(m["N"], size=m["n_case"], replace=False)
            cpp = np.bincount(m["pos_of_allele"][ca], minlength=len(m["positions"])).astype(float)
            a = m["M"] @ cpp
            p = np.where(m["nbhd_total"] > 0,
                         hypergeom.sf(a - 1, m["N"], m["n_case"], m["nbhd_total"]),
                         1.0).min()
            if p < gm:
                gm = p
        global_min[s] = gm
    global_min.sort()
    ranks = np.searchsorted(global_min, np.asarray(obs_minp), side="right")
    fwer = (1.0 + ranks) / (int(n_sims) + 1.0)
    thresh_p = float(np.quantile(global_min, 0.05))
    return fwer, float(-np.log10(max(thresh_p, 1e-300)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", default=DEF_VARIANTS)
    ap.add_argument("--out", default=DEF_OUT)
    ap.add_argument("--radius", type=float, default=15.0)
    ap.add_argument("--pae-cutoff", type=float, default=15.0)
    ap.add_argument("--plddt-cutoff", type=float, default=0.0,
                    help="0 matches the paper (keeps disordered residues in the background)")
    ap.add_argument("--sims", type=int, default=1000, help="pooled FWER replicates (0 = skip)")
    ap.add_argument("--limit", type=int, default=0, help="only first N proteins (testing)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    df = pd.read_csv(args.variants, sep="\t")
    genes = df[["uniprot_id", "gene_name"]].drop_duplicates().values.tolist()
    if args.limit:
        genes = genes[: args.limit]
    print(f"{len(genes)} proteins to process")

    records, models = [], []
    for i, (acc, gene) in enumerate(genes, 1):
        sub = df[df.uniprot_id == acc]
        counts = (sub.groupby("aa_pos")[["ac_case", "ac_control"]].sum()
                  .reset_index().rename(columns={"aa_pos": "aa_pos"}))
        try:
            pdb_gz, pae = F.fetch_alphafold(acc, seq_len=None)
            m = FS.build_neighborhood_model(
                counts, pdb_gz, pae, radius=args.radius,
                pae_cutoff=args.pae_cutoff, plddt_cutoff=args.plddt_cutoff)
        except Exception as e:
            print(f"  [{i}/{len(genes)}] {gene} ({acc}) skipped: {e}")
            continue
        minp = float(m["p_obs"].min())
        top_pos = int(m["positions"][int(np.argmin(m["p_obs"]))])
        records.append({
            "uniprot_id": acc, "gene_name": gene, "top_aa_pos": top_pos,
            "min_p": minp, "neglog10_min_p": -np.log10(max(minp, 1e-300)),
            "n_tested": len(m["positions"]),
            "n_case": m["n_case"], "n_control": m["n_ctrl"],
        })
        if args.sims > 0:
            m["pos_of_allele"] = np.repeat(np.arange(len(m["positions"])),
                                           m["tot"].astype(int))
            m["_minp"] = minp
            models.append(m)
        if i % 25 == 0:
            print(f"  [{i}/{len(genes)}] {gene}: min p = {minp:.2e}")

    res = pd.DataFrame(records)
    if res.empty:
        print("No proteins produced results."); return

    res["fdr"] = bh_fdr(res["min_p"].values)
    if args.sims > 0 and models:
        print(f"pooled permutation FWER: {args.sims} replicates over {len(models)} proteins…")
        obs = [m["_minp"] for m in models]
        fwer, thr = pooled_fwer(models, obs, args.sims, seed=args.seed)
        # models were appended in lockstep with `records`, so they share order
        fwer_by_acc = {rec["uniprot_id"]: fwer[j] for j, rec in enumerate(records)}
        res["fwer"] = res["uniprot_id"].map(fwer_by_acc)
        res["neglog10_fwer_thresh"] = thr
    else:
        res["fwer"] = np.nan
        res["neglog10_fwer_thresh"] = np.nan

    res = res.sort_values("min_p").reset_index(drop=True)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    res.to_csv(args.out, sep="\t", index=False)
    n_fdr = int((res["fdr"] < 0.05).sum())
    n_fwer = int((res["fwer"] < 0.05).sum()) if res["fwer"].notna().any() else 0
    print(f"wrote {args.out}: {len(res)} proteins, {n_fdr} at FDR<0.05, {n_fwer} at FWER<0.05")


if __name__ == "__main__":
    main()
