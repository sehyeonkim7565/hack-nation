"""
Module 02 step A: deterministic molecular-target gate.

Klebsiella pneumoniae carries a chromosomal, species-intrinsic SHV class A
beta-lactamase. This confers resistance to aminopenicillins (ampicillin,
amoxicillin) regardless of any acquired resistance gene content, and this
fact is not something a statistical model should have to "learn" from a
tiny cohort -- it is textbook clinical microbiology. This gate is applied
BEFORE and independently of the logistic-regression predictor in
predict.py, and always wins over a model prediction for the antibiotics
it covers.

Scope note: none of the 3 antibiotics modeled statistically in this pilot
(trimethoprim/sulfamethoxazole, gentamicin, cefepime) are gated by this
rule -- ampicillin was excluded from the pilot's modeled set because its
label distribution is almost entirely resistant (see README). The gate is
still implemented and exercised here because it is a required, always-on
safety behavior of the overall decision system, independent of which 3
antibiotics happen to be statistically modeled this pilot.
"""

INTRINSIC_RESISTANCE_GATES = {
    "ampicillin": {
        "call": "likely to fail",
        "mechanism": "Chromosomal, species-intrinsic SHV class A beta-lactamase "
                     "(present in essentially all K. pneumoniae).",
        "confidence": 0.99,
    },
    "amoxicillin": {
        "call": "likely to fail",
        "mechanism": "Chromosomal, species-intrinsic SHV class A beta-lactamase "
                     "(present in essentially all K. pneumoniae).",
        "confidence": 0.99,
    },
}


def check_gate(antibiotic: str, organism: str = "Klebsiella pneumoniae") -> dict | None:
    """Return a gate override dict if this antibiotic/organism combination
    is covered by a deterministic rule, else None (fall through to the
    statistical predictor / evidence pipeline)."""
    if organism.strip().lower() != "klebsiella pneumoniae":
        return None
    key = antibiotic.strip().lower()
    if key in INTRINSIC_RESISTANCE_GATES:
        gate = INTRINSIC_RESISTANCE_GATES[key]
        return {
            "antibiotic": antibiotic,
            "call": gate["call"],
            "confidence": gate["confidence"],
            "evidence_type": "deterministic_gate",
            "explanation": gate["mechanism"],
        }
    return None
