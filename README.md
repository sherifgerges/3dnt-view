# 3DNT web tool

A Streamlit app to run 3DNT in a browser.

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


## ASD genome-wide mode

The **"ASD (genome-wide)"** page reproduces the paper's ASD analysis with this
site's own Fisher engine. It has two parts:

- a **Manhattan plot** over the 669 proteins the 3DNT paper tested — one point per
  protein (its most significant neighborhood), with FDR 0.05 (dashed) and FWER
  0.05 (dotted) threshold lines and labeled significant genes;
- an **individual protein** selector — pick any of the 669 genes and its 3DNT
  result (table + enrichment-colored structure) is computed live.

Bundled data lives in `asd_data/`:
- `ASD_variants.tsv` — per-residue case/control variants for the 669 proteins
  (shipped with the app so the drill-down works anywhere).
- `asd_web_results.tsv` — the Manhattan summary. **Generate it once** with:

  ```bash
  python asd_batch.py            # all 669 proteins, 1000 permutations
  python asd_batch.py --sims 0   # FDR only, faster (skip permutation FWER)
  python asd_batch.py --limit 20 # quick test on 20 proteins
  ```

  `asd_batch.py` fetches each protein's AlphaFold structure, so run it on a
  machine with internet (structures are cached in `afdb_cache/`, so it's
  resumable). It computes genome-wide, gene-level BH-FDR and a **pooled-permutation
  FWER** (labels permuted within every protein each replicate; the global minimum
  p across proteins is the null max-statistic). If `asd_web_results.tsv` is
  missing, the page still works — it just shows the individual-protein tool and a
  note to run the batch.

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
