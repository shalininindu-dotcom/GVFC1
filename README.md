# Genomic Variant Frequency Calculator

A small web app that takes a VCF file (or raw genotype counts) and reports, for each variant:

- allele frequencies (reference `p`, alternate `q`)
- genotype frequencies (HomRef, Het, HomAlt)
- Minor Allele Frequency (MAF) and a common / low-frequency / rare label
- a Hardy-Weinberg Equilibrium (HWE) chi-square test

Built with Python and Streamlit.

## Run it on your computer

```bash
git clone <your-repo-url>
cd <your-repo-folder>
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

Then open the address Streamlit prints (usually http://localhost:8501).

## Run the tests

```bash
pip install pytest
python -m pytest -q
```

## Project layout

| File | Purpose |
|---|---|
| `variants.py` | All the calculation logic (frequencies, MAF, HWE, VCF parsing) |
| `app.py` | The Streamlit web interface |
| `test_variants.py` | Unit tests for the calculations |
| `sample_data/sample.vcf` | Small test VCF (10 samples, 8 variants) |
| `requirements.txt` | Python packages the app needs |

## How the maths works

For N individuals with counts n_AA, n_Aa, n_aa:

- p = (2·n_AA + n_Aa) / 2N, q = (2·n_aa + n_Aa) / 2N, MAF = min(p, q)
- Expected genotype counts under HWE: p²N, 2pqN, q²N
- Chi-square = Σ (observed − expected)² / expected, with 1 degree of freedom

## Limitations

- Only bi-allelic variants are analysed (multi-allelic sites are skipped).
- Missing (`./.`) and haploid calls are excluded from the counts.
- The HWE test is a simple chi-square; for very small samples or rare variants an exact test is more accurate.

## Team

- Person A: _name_ , Person B: _name_
