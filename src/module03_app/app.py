"""
Module 03: Genome Firewall decision-report app (Streamlit).

Upload a K. pneumoniae genome FASTA -> Module 01 (BLAST vs. the ResFinder
acquired-resistance-gene DB) detects genes -> Module 02 (deterministic
gate + gene-evidence rules) produces a per-antibiotic call with a
calibrated confidence, an evidence type, and a forced no-call whenever
evidence is weak or confidence is low.

DEFENSIVE-USE NOTE: this tool only reads and classifies existing genome
sequence. It does not design, modify, or synthesize anything. All output
is decision support only, per the fixed disclaimer shown on every result.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import streamlit as st

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src" / "module01_genome_reader"))
sys.path.insert(0, str(REPO_ROOT / "src" / "module02_predictor"))

from detect_amr_genes import (  # noqa: E402
    run_blast, collapse_to_best_hit_per_locus, load_gene_phenotypes,
    MIN_COVERAGE, PHENOTYPES_TSV,
)
from predict import predict_all, DISCLAIMER  # noqa: E402

BLAST_DB_DIR = REPO_ROOT / "data" / "blast_db"
BLAST_DB_PREFIX = BLAST_DB_DIR / "resfinder"


@st.cache_resource
def ensure_blast_db_built():
    if (BLAST_DB_DIR / "resfinder.nsq").exists() or (BLAST_DB_DIR / "resfinder.nal").exists():
        return
    BLAST_DB_DIR.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["makeblastdb", "-in", str(REPO_ROOT / "reference_db" / "resfinder_all.fsa"),
         "-dbtype", "nucl", "-out", str(BLAST_DB_PREFIX),
         "-title", "ResFinder acquired AMR gene DB"],
        check=True, capture_output=True,
    )


ensure_blast_db_built()

st.set_page_config(page_title="Genome Firewall - K. pneumoniae AMR Pilot", layout="wide")

st.title("Genome Firewall: K. pneumoniae AMR Decision Support (Pilot)")

st.error(
    "**This is a preliminary pilot built on n=25 lab-confirmed genomes due to time "
    "constraints of a 24-hour hackathon. Results should be treated as a proof-of-concept "
    "for the pipeline architecture, not as clinically meaningful predictions. Production "
    "deployment would require the full BV-BRC AMR Panel cohort (~6,700 K. pneumoniae genomes).**"
)
st.warning(
    "Of the 25 labeled genomes in this pilot's cohort, only **1** had usable matching "
    "nucleotide sequence available for gene detection (see README for the full data audit). "
    "This app's gene-detection path is fully functional and genome-count-agnostic, but the "
    "statistical model behind it was validated on a real cross-antibiotic phenotype baseline, "
    "**not** on genome-derived features -- see 'Evidence type' below for what that means for "
    "any single prediction."
)

uploaded_file = st.file_uploader("Upload a genome FASTA (.fna/.fasta/.fa)", type=["fna", "fasta", "fa"])
organism = st.text_input("Organism (only Klebsiella pneumoniae is in-scope for this pilot)",
                          value="Klebsiella pneumoniae")

if uploaded_file is not None:
    with tempfile.NamedTemporaryFile(suffix=".fna", delete=False) as tmp:
        tmp.write(uploaded_file.read())
        tmp_path = Path(tmp.name)

    if organism.strip().lower() != "klebsiella pneumoniae":
        st.error(
            "This pilot is scoped to Klebsiella pneumoniae only (see project scope rules). "
            "Predictions for other organisms are not supported and will not be produced."
        )
    else:
        with st.spinner("Running BLAST against the ResFinder acquired-resistance-gene database..."):
            raw_hits = [h for h in run_blast(tmp_path) if h["coverage"] >= MIN_COVERAGE]
            best_hits = collapse_to_best_hit_per_locus(raw_hits)
            gene_phenotypes = load_gene_phenotypes(PHENOTYPES_TSV)

            gene_hits = []
            for hit in best_hits:
                gene_acc = hit["sseqid"]
                pheno = gene_phenotypes.get(gene_acc, {"class": "unknown", "antibiotics": []})
                gene_hits.append({
                    "gene_accession": gene_acc,
                    "gene_name": gene_acc.split("_")[0],
                    "identity_pct": round(hit["pident"], 2),
                    "coverage_pct": round(hit["coverage"], 2),
                    "drug_class": pheno["class"],
                    "antibiotics": ";".join(pheno["antibiotics"]),
                })

        st.subheader("Detected acquired resistance genes (Module 01)")
        if gene_hits:
            st.dataframe(gene_hits, use_container_width=True)
        else:
            st.info("No acquired resistance genes detected above thresholds "
                     "(>=90% identity, >=60% coverage).")

        st.subheader("Per-antibiotic decision (Module 02 + Module 03)")
        predictions = predict_all(gene_hits, organism)

        for pred in predictions:
            call = pred["call"]
            if call == "no-call":
                badge = st.warning
            elif "fail" in call:
                badge = st.error
            else:
                badge = st.success

            with st.container(border=True):
                st.markdown(f"### {pred['antibiotic']}")
                col1, col2, col3 = st.columns(3)
                col1.metric("Call", pred["call"])
                col2.metric("Confidence", f"{pred['confidence']:.2f}")
                col3.metric("Evidence type", pred["evidence_type"])
                st.markdown(f"**Explanation:** {pred['explanation']}")

    Path(tmp_path).unlink(missing_ok=True)

st.divider()
st.info(DISCLAIMER)
