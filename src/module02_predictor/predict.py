"""
Module 02 step C: inference-time prediction for a freshly uploaded genome.

This is the function Module 03 (the Streamlit app) calls. It is
deliberately NOT the same code path as train_loocv.py's cross-antibiotic
statistical baseline, because that baseline needs OTHER real lab
susceptibility results for the same genome as its input features -- which
a freshly uploaded, never-before-tested genome does not have by
definition. Using it here would be silently fabricating input data.

So for a fresh genome upload, only two evidence sources are actually
available, and this module is honest about that:
  1. Deterministic gate (gate.py) -- always checked first, wins if it fires.
  2. Direct resistance-gene detection (Module 01 BLAST/ResFinder hits) --
     if a gene mapped to this antibiotic was found, that's evidence type
     "known_resistance_gene".
  3. Otherwise: "no_signal" -- this pilot has no genome-derived statistical
     model (see data_prep.py/train_loocv.py docstrings for why), so the
     honest output is a forced no-call, not a guess.

The retrospective LOOCV cross-antibiotic model from train_loocv.py is
reported separately in the evaluation section as a real, LOOCV-validated
statistical baseline over the labeled cohort -- it is not wired into this
upload-time inference path.
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gate import check_gate  # noqa: E402

TARGET_ANTIBIOTICS = [
    "trimethoprim/sulfamethoxazole",
    "gentamicin",
    "cefepime",
]

# Minimum confidence to avoid a forced no-call.
NO_CALL_CONFIDENCE_THRESHOLD = 0.70

DISCLAIMER = (
    "This is a decision-support estimate only. Results MUST be confirmed "
    "by standard laboratory antimicrobial susceptibility testing before "
    "any clinical or treatment decision is made."
)


def predict_antibiotic(antibiotic: str, gene_hits: list[dict], organism: str = "Klebsiella pneumoniae") -> dict:
    gate_result = check_gate(antibiotic, organism)
    if gate_result is not None:
        return {
            "antibiotic": antibiotic,
            "call": gate_result["call"],
            "confidence": gate_result["confidence"],
            "evidence_type": "deterministic_gate",
            "explanation": gate_result["explanation"],
            "no_call": False,
        }

    matching_genes = [
        h for h in gene_hits
        if antibiotic.lower() in [a.strip().lower() for a in h.get("antibiotics", "").split(";")]
    ]
    if matching_genes:
        gene_list = ", ".join(f"{h['gene_name']} ({h['identity_pct']}% id)" for h in matching_genes)
        confidence = min(0.95, 0.75 + 0.05 * len(matching_genes))
        return {
            "antibiotic": antibiotic,
            "call": "likely to fail",
            "confidence": confidence,
            "evidence_type": "known_resistance_gene",
            "explanation": f"Acquired resistance gene(s) detected with a mapped phenotype "
                            f"to this antibiotic: {gene_list}.",
            "no_call": False,
        }

    return {
        "antibiotic": antibiotic,
        "call": "no-call",
        "confidence": 0.0,
        "evidence_type": "no_signal",
        "explanation": (
            "No deterministic rule applies and no acquired resistance gene mapped to this "
            "antibiotic was detected. This pilot has no validated genome-derived statistical "
            "model (see README/limitations: only 1 of 25 labeled genomes had usable sequence "
            "for training). Absence of a detected gene is not proof of susceptibility -- "
            "defaulting to no-call rather than guessing."
        ),
        "no_call": True,
    }


def predict_all(gene_hits: list[dict], organism: str = "Klebsiella pneumoniae") -> list[dict]:
    results = [predict_antibiotic(ab, gene_hits, organism) for ab in TARGET_ANTIBIOTICS]
    for r in results:
        if not r["no_call"] and r["confidence"] < NO_CALL_CONFIDENCE_THRESHOLD:
            r["call"] = "no-call"
            r["no_call"] = True
            r["explanation"] += " [Confidence below threshold -- forced to no-call.]"
    return results
