"""
app.py - Streamlit web app for the Genomic Variant Frequency Calculator.

Run locally:   streamlit run app.py
"""

import io

import pandas as pd
import streamlit as st

from variants import analyze_counts, analyze_vcf

st.set_page_config(page_title="Genomic Variant Frequency Calculator", page_icon="🧬", layout="wide")

st.title("🧬 Genomic Variant Frequency Calculator")
st.write(
    "Compute allele frequencies, genotype frequencies, Minor Allele Frequency (MAF) "
    "and a Hardy-Weinberg Equilibrium (HWE) check for genetic variants."
)

mode = st.sidebar.radio("Input mode", ["Upload a VCF file", "Enter genotype counts"])

with st.sidebar.expander("What do these terms mean?"):
    st.markdown(
        """
- **Allele frequency**: how common an allele is among all copies in the population.
  p = reference allele, q = alternate allele (p + q = 1).
- **Genotype frequency**: share of people who are HomRef (AA), Het (Aa) or HomAlt (aa).
- **MAF**: frequency of the *less* common allele. Variants are often called
  common (>= 5%), low-frequency (1-5%) or rare (< 1%).
- **HWE**: compares observed genotype counts with the p², 2pq, q² expected in a
  randomly mating population. A small p-value (< 0.05) is a warning sign for
  genotyping errors, selection or population structure.
        """
    )


# ---------------------------------------------------------------------------
# Mode 1: VCF upload
# ---------------------------------------------------------------------------
if mode == "Upload a VCF file":
    st.header("Upload a VCF file")

    maf_threshold = st.sidebar.slider("MAF filter threshold", 0.0, 0.5, 0.05, 0.01)

    uploaded = st.file_uploader("Choose a .vcf file", type=["vcf"])
    use_sample = st.checkbox("Use the included sample file (sample_data/sample.vcf)")

    source = None
    if uploaded is not None:
        source = io.StringIO(uploaded.getvalue().decode("utf-8", errors="replace"))
    elif use_sample:
        with open("sample_data/sample.vcf", encoding="utf-8") as fh:
            source = io.StringIO(fh.read())

    if source is not None:
        try:
            df = analyze_vcf(source, maf_threshold=maf_threshold)
        except Exception as exc:  # show a friendly message instead of a traceback
            st.error(f"Could not read this VCF file: {exc}")
            st.stop()

        if df.empty:
            st.warning(
                "No usable bi-allelic variants with genotype (GT) calls were found. "
                "Check that the file is a valid VCF with sample columns."
            )
            st.stop()

        # Summary metrics
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Variants analysed", len(df))
        c2.metric("Pass MAF filter", int(df["Passes_MAF_filter"].sum()))
        c3.metric("Rare variants (MAF < 1%)", int((df["Category"] == "rare").sum()))
        c4.metric("Fail HWE (p < 0.05)", int((~df["HWE_pass"]).sum()))

        # Filters
        st.subheader("Results table")
        f1, f2 = st.columns(2)
        only_pass_maf = f1.checkbox("Show only variants passing the MAF filter")
        only_pass_hwe = f2.checkbox("Show only variants passing HWE")
        view = df
        if only_pass_maf:
            view = view[view["Passes_MAF_filter"]]
        if only_pass_hwe:
            view = view[view["HWE_pass"]]
        st.dataframe(view, hide_index=True)

        st.download_button(
            "Download results as CSV",
            view.to_csv(index=False).encode("utf-8"),
            file_name="variant_frequencies.csv",
            mime="text/csv",
        )

        # Charts
        st.subheader("Charts")
        left, right = st.columns(2)
        with left:
            st.caption("MAF distribution (number of variants per MAF bin)")
            bins = pd.cut(df["MAF"], bins=[0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5], include_lowest=True)
            hist = bins.value_counts().sort_index()
            hist.index = hist.index.astype(str)
            st.bar_chart(hist)
        with right:
            st.caption("Allele frequency (q) per variant")
            labels = df["ID"].where(df["ID"] != ".", df["CHROM"].astype(str) + ":" + df["POS"].astype(str))
            st.bar_chart(pd.Series(df["Alt_freq (q)"].values, index=labels.values))
    else:
        st.info("Upload a VCF file, or tick the sample file box, to get started.")


# ---------------------------------------------------------------------------
# Mode 2: manual genotype counts
# ---------------------------------------------------------------------------
else:
    st.header("Enter genotype counts")
    c1, c2, c3 = st.columns(3)
    n_hom_ref = c1.number_input("Homozygous reference (AA)", min_value=0, value=60, step=1)
    n_het = c2.number_input("Heterozygous (Aa)", min_value=0, value=30, step=1)
    n_hom_alt = c3.number_input("Homozygous alternate (aa)", min_value=0, value=10, step=1)

    if n_hom_ref + n_het + n_hom_alt == 0:
        st.warning("Enter at least one individual.")
        st.stop()

    res = analyze_counts(int(n_hom_ref), int(n_het), int(n_hom_alt))

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Reference allele freq (p)", f"{res['p_ref']:.4f}")
    m2.metric("Alternate allele freq (q)", f"{res['q_alt']:.4f}")
    m3.metric("MAF", f"{res['maf']:.4f}")
    m4.metric("Individuals", res["n_samples"])

    st.subheader("Genotype frequencies")
    geno = pd.DataFrame(
        {
            "Genotype": ["HomRef (AA)", "Het (Aa)", "HomAlt (aa)"],
            "Observed count": res["hwe_observed"],
            "Expected count (HWE)": [round(x, 2) for x in res["hwe_expected"]],
            "Observed frequency": [round(res["hom_ref"], 4), round(res["het"], 4), round(res["hom_alt"], 4)],
        }
    )
    st.dataframe(geno, hide_index=True)
    st.bar_chart(geno.set_index("Genotype")[["Observed count", "Expected count (HWE)"]])

    st.subheader("Hardy-Weinberg Equilibrium test")
    st.write(f"Chi-square = **{res['hwe_chi_square']:.4f}**, p-value = **{res['hwe_p_value']:.4f}**")
    if res["in_hwe"]:
        st.success("No significant deviation from Hardy-Weinberg Equilibrium (p >= 0.05).")
    else:
        st.error("Significant deviation from Hardy-Weinberg Equilibrium (p < 0.05).")

st.divider()
st.caption("Built for the Biopython project course. Multi-allelic sites and missing/haploid calls are skipped.")
