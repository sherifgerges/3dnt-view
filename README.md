# 3DNT web tool

A small Streamlit app: enter a gene/protein, upload case + control variant lists
(`E457K`- or `p.Glu457Lys`-style, one per line), and run a **3D neighborhood
enrichment test** on the protein's AlphaFold structure. Residues are colored on
the structure by how strongly the neighborhoods around them are enriched for
case variants.

The app is **self-contained** — pure Python, no external repo or reference data.
It resolves the protein via UniProt and downloads the AlphaFold model + PAE at
runtime (cached in `afdb_cache/`).

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

Then open the URL Streamlit prints (usually http://localhost:8501).

## How it works

```
gene/UniProt ─► UniProt REST ─► accession + sequence
                     └─► AlphaFold DB ─► structure (.pdb.gz) + PAE (cached in afdb_cache/)
case.txt / control.txt ─► parse variants ─► per-residue ac_case / ac_control
              │
              ▼
   fisher_scan.run_fisher_3dnt
     • all-atom min residue–residue distances, inflated to 1000 Å when both
       directions of the PAE exceed the cutoff
     • residues with mean pLDDT ≤ cutoff dropped
     • one-sided Fisher's exact test per neighborhood center (case enrichment)
     • per-residue score = −log10(min p over neighborhoods containing it)
              ▼
   results table + TSV  ·  structure colored white→red by enrichment
```

Files:
- `app.py` — Streamlit UI and orchestration
- `sir_fetch.py` — gene→UniProt resolution + AlphaFold download (cached)
- `variants.py` — parse/validate variant files, build per-residue counts
- `fisher_scan.py` — the self-contained Fisher 3D neighborhood test
- `examples/` — bundled ATP2B2 (Q01814) case/control lists for the demo button

## Deploy to Streamlit Community Cloud

1. Put this `3dnt_webtool/` folder in a GitHub repo (it can be the repo root, or
   a subfolder).
2. Go to https://share.streamlit.io, sign in with GitHub, and click **New app**.
3. Select the repo/branch and set **Main file path** to `app.py` (include the
   subfolder if it isn't the repo root, e.g. `3dnt_webtool/app.py`).
4. Deploy. Streamlit installs `requirements.txt` and gives you a
   `https://<name>.streamlit.app` URL.
5. To keep it private: in the app's **Settings → Sharing**, set viewers to
   specific emails (e.g. Hilary's) instead of public.

Notes:
- Outbound network to UniProt and the AlphaFold DB is allowed on Streamlit
  Cloud, so structure fetching works out of the box.
- `afdb_cache/` is ephemeral on the cloud; the bundled ATP2B2 model is committed
  so the example button is instant, and other proteins are fetched on demand.

## Notes / limitations

- **Significance is within-protein** (Bonferroni / p-value thresholds across this
  protein's tested centers), not genome-wide — labeled in the UI.
- **Single-fragment proteins only** — very large (multi-fragment AlphaFold)
  proteins aren't handled yet.
- The center residue is excluded from its own neighborhood, matching the
  reference PyMOL pipeline.

## Citation

See the **Citations** page in the app (bioRxiv 2025.08.25.672202 and
medRxiv 2026.05.29.26354366).
