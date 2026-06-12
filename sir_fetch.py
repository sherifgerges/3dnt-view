"""Resolve a gene/protein name to a UniProt entry and fetch its AlphaFold model."""
import os
import re
import gzip
import requests

CACHE = os.environ.get("AFDB_CACHE", os.path.join(os.path.dirname(__file__), "afdb_cache"))
os.makedirs(CACHE, exist_ok=True)

ACCESSION_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")


def resolve_uniprot(query):
    """Return a list of {accession, gene, sequence} matches (human, reviewed).

    If `query` already looks like a UniProt accession, fetch that entry directly.
    Otherwise search by gene symbol. The caller can let the user pick if >1.
    """
    query = query.strip()
    if ACCESSION_RE.match(query.upper()):
        url = (f"https://rest.uniprot.org/uniprotkb/{query.upper()}"
               f"?fields=accession,gene_primary,sequence&format=json")
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        d = r.json()
        return [{
            "accession": d["primaryAccession"],
            "gene": (d.get("genes", [{}])[0].get("geneName", {}) or {}).get("value", query),
            "sequence": d["sequence"]["value"],
        }]

    url = ("https://rest.uniprot.org/uniprotkb/search"
           f"?query=gene:{query}+AND+organism_id:9606+AND+reviewed:true"
           "&fields=accession,gene_primary,sequence&format=json&size=10")
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    out = []
    for d in r.json().get("results", []):
        out.append({
            "accession": d["primaryAccession"],
            "gene": (d.get("genes", [{}])[0].get("geneName", {}) or {}).get("value", ""),
            "sequence": d["sequence"]["value"],
        })
    return out


def fetch_alphafold(accession, seq_len=None):
    """Download AlphaFold model (.pdb) + PAE (.json) for a protein.

    Uses the AlphaFold DB API to discover the current file URLs, so we don't
    hardcode a model version. The API can return several entries per protein
    (different model versions / AF2 vs AF3), so we pick the best FULL-LENGTH one
    using each entry's uniprotStart/uniprotEnd. Only a genuine multi-fragment
    protein (no single entry covering the whole sequence) raises.
    Returns (pdb_gz_path, pae_path); the PDB is gzipped to match the SIR convention.
    """
    pdb_gz = os.path.join(CACHE, f"AF-{accession}-model.pdb.gz")
    pae_path = os.path.join(CACHE, f"AF-{accession}-pae.json")
    if os.path.exists(pdb_gz) and os.path.exists(pae_path):
        return pdb_gz, pae_path

    api = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    r = requests.get(api, timeout=60)
    if r.status_code != 200:
        raise FileNotFoundError(
            f"AlphaFold has no entry for {accession} (API status {r.status_code}).")
    entries = r.json()
    if not entries:
        raise FileNotFoundError(f"AlphaFold returned no models for {accession}.")

    def start(e): return int(e.get("uniprotStart", 1) or 1)
    def end(e):   return int(e.get("uniprotEnd", 0) or 0)

    # prefer entries that begin at residue 1, then the one covering the most;
    # break ties by latest version / newest model
    candidates = [e for e in entries if start(e) == 1] or entries
    best = max(candidates, key=lambda e: (end(e), str(e.get("latestVersion", "")),
                                          str(e.get("modelCreatedDate", ""))))

    # only complain if even the best entry doesn't reach (near) the protein's end
    if seq_len and end(best) and end(best) < seq_len - 1:
        raise FileNotFoundError(
            f"{accession}: AlphaFold's best model only covers residues "
            f"{start(best)}-{end(best)} of {seq_len} (multi-fragment); "
            "the MVP handles single-fragment proteins only.")

    pdb_url = best.get("pdbUrl")
    pae_url = best.get("paeDocUrl") or best.get("paeUrl")
    if not pdb_url or not pae_url:
        raise FileNotFoundError(f"AlphaFold entry for {accession} is missing pdb/PAE URLs.")

    rp = requests.get(pdb_url, timeout=120); rp.raise_for_status()
    with gzip.open(pdb_gz, "wt") as fh:
        fh.write(rp.text)
    ra = requests.get(pae_url, timeout=120); ra.raise_for_status()
    with open(pae_path, "w") as fh:
        fh.write(ra.text)

    return pdb_gz, pae_path
