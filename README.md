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

Enter a gene/UniProt ID and upload `case.txt` / `control.txt` (variant lists).
The method is three steps:

1. **Map case and control variants to the protein structure** (AlphaFold model,
   fetched from the AlphaFold DB via UniProt).
2. **For the neighborhood of radius R centered around each residue**, compare the
   case:control ratio in the neighborhood to the case:control ratio in the rest
   of the protein with a one-sided **Fisher's exact test**.
3. **Score each residue** as −log10 of the smallest p-value across the
   neighborhoods that contain it.

Neighborhoods use all-atom residue–residue distances, with a pair pushed apart
(distance inflated) when both directions of the AlphaFold **PAE** exceed the
cutoff — i.e. their relative position is poorly predicted. Residues with mean
pLDDT below the cutoff are dropped.

Output: a results table + downloadable TSV, and the structure colored white→red
by enrichment. Multiple-testing significance uses a permutation FWER (see the UI).

Files:
- `app.py` — Streamlit UI and orchestration
- `sir_fetch.py` — gene→UniProt resolution + AlphaFold download (cached)
- `variants.py` — parse/validate variant files, build per-residue counts
- `fisher_scan.py` — the self-contained Fisher 3D neighborhood test
- `examples/` — bundled ATP2B2 (Q01814) case/control lists for the demo button


## Notes / limitations

- **Significance is within-protein** (Bonferroni / p-value thresholds across this
  protein's tested centers), not genome-wide — labeled in the UI.
- **Single-fragment proteins only** — very large (multi-fragment AlphaFold)
  proteins aren't handled yet.
- The center residue is excluded from its own neighborhood, matching the
  reference PyMOL pipeline.

## Citations

If you use this tool, please cite:

**Genetic and structural evidence links Ca²⁺ dysregulation and ATP2B2 to neuropsychiatric illness.**
Sherif Gerges, Nikolaj Catois Straarup, Mohamed A. El-Brolosy, F. Kyle Satterstrom,
Nolan Kamitaki, Jiayi Yuan, Emi Ling, Carmen Gelze, Raozhou Lin, Melissa Goldman,
Curtis Mello, Tarjinder Singh, The Autism Sequencing Consortium, Jonathan S. Weissman,
Sabina Berretta, Jen Q. Pan, Hilary Finucane, Charlott Stock, Poul Nissen,
Steven A. McCarroll, Mark Daly. *bioRxiv* 2025.08.25.672202.
https://www.biorxiv.org/content/10.1101/2025.08.25.672202v5

**Systematic identification of disease-associated 3D neighborhoods in protein structures.**
Emily Nason, Sherif Gerges, F. Kyle Satterstrom, Bram L. Gorissen, Ruqi Liao,
Georgia Panagiotaropoulou, Jeremy Guez, The Autism Sequencing Consortium,
Konrad J. Karczewski, Mark Daly, Hilary Finucane. *medRxiv* 2026.05.29.26354366.
https://www.medrxiv.org/content/10.64898/2026.05.29.26354366v1
