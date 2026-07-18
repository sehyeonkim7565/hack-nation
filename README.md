# Genome Firewall — *K. pneumoniae* AMR Decision-Support Pilot

**Live app:** https://hack-nation-6vjg298pfnsih9rwfi9fov.streamlit.app/

> **This is a preliminary pilot built on n=25 lab-confirmed genomes due to time constraints of a 24-hour hackathon. Results should be treated as a proof-of-concept for the pipeline architecture, not as clinically meaningful predictions. Production deployment would require the full BV-BRC AMR Panel cohort (~6,700 *K. pneumoniae* genomes).**

**This is an extremely limited pilot cohort (n=25).** The original target was 1,000–3,000 genomes; what was actually available and lab-confirmed is ~1–2% of that. Every number in this document should be read with that in mind. This tool is a **defensive, decision-support prototype only**. It never designs, modifies, or synthesizes an organism. Every prediction carries a fixed disclaimer: **all results must be reconfirmed by standard laboratory antimicrobial susceptibility testing before any clinical or treatment decision is made.**

---

## 0. Data audit (read this first — it changed the project's scope)

Before building anything, the three input files were opened and inspected directly. Two findings fundamentally shaped what could honestly be built.

### 0.1 Corrupted / injected content in the raw files

Both `kpneu_amr_phenotypes.tsv` and `kpneu_genomes.fna` contain leftover **Windows terminal transcript noise** from however the data was exported via the PATRIC CLI (`p3-genome-fasta`):

- Command echoes: `C:\Program Files\PATRIC\cli\bin>p3-genome-fasta --contig <genome_id>`
- Repeated, EUC-KR-encoded batch-confirmation prompts: `일괄 작업을 끝내시겠습니까 (Y/N)?` ("Shall I finish the batch job (Y/N)?")

This is not valid FASTA or TSV content, and on first read it looked like it could be a prompt-injection attempt aimed at an AI agent processing the file. It was **not acted on** as an instruction. The parsing code (`split_genomes.py`) explicitly strips this noise and, where possible, uses the `--contig <genome_id>` markers to correctly attribute sequence to genome IDs — which is what revealed finding 0.2 below.

### 0.2 The FASTA file does not contain 25 genomes

`kpneu_genomes.fna` has 88 contigs total, not the expected ~92, and — critically — **real, usable sequence exists for only 2 distinct `genome_id`s**, not 25:

| genome_id | contigs | bp | in scope? |
|---|---|---|---|
| `1226680.5` | 2 (chromosome + plasmid) | 5.51 Mb | ✅ *K. pneumoniae* (Ecl8), matches metadata |
| `1389422.3` | 86 (accession `AYQE01...`) | 5.63 Mb | ❌ metadata lists this genome as ***Escherichia coli* LAU-EC1**, not *K. pneumoniae* — excluded per the species-scope rule |

Two further genome fetches are referenced in the leftover CLI transcript (`1250520.6`, `1304916.4`) but their sequence output is missing/truncated in the file — those genomes have **zero** sequence data available.

**Net result: of the 25 labeled *K. pneumoniae* genomes, exactly 1 (`1226680.5`) has real, species-matched nucleotide sequence.** The other 24 have phenotype labels (from `kpneu_amr_phenotypes.tsv`) but no corresponding FASTA sequence in this dataset. This is reported here in full rather than worked around silently, and it is the reason Module 02 below is not the genome-feature model originally scoped — see §2.2.

### 0.3 Label/FASTA matching summary

- `kpneu_amr_metadata.tsv`: 42,829 genome-level rows, 4 columns (`genome_id`, `genome_name`, AMR summary, evidence). Only used to check species; ~21% blank in the AMR summary column.
- `kpneu_amr_phenotypes.tsv`: 3,760 rows (1 trailing corrupt line dropped, see §0.1). Exactly **25 unique genome_ids** have any non-blank `resistant_phenotype` — matches the intended pilot size.
- All 25 labeled genome_ids are present in the metadata file (100% match). Only 1 of the 25 is present with real sequence in the FASTA file (4% match) — this is the actual, audited number, not the assumed one.
- Chosen 3 antibiotics, confirmed exact label counts: trimethoprim/sulfamethoxazole (R16/S9), gentamicin (R7/S17/I1), cefepime (R17/S3/I5).

---

## 1. Module 01 — Genome Reader

**Input:** genome FASTA. **Output:** detected acquired resistance genes + a genome × gene presence/absence matrix.

### 1.1 AMRFinderPlus could not be installed

The specified tool (AMRFinderPlus via conda, `--organism Klebsiella`) could not be set up in this environment:
- No `conda`/`mamba` is available (the specified install path).
- Direct network egress to `github.com` and `ftp.ncbi.nlm.nih.gov` — where AMRFinderPlus's binary and curated reference database are hosted — is blocked by the sandbox's egress policy (confirmed `403` policy denial on both).

### 1.2 Fallback used (as instructed): BLAST against a real curated reference DB

Per the project's stated fallback plan, gene detection instead uses `blastn` against the **ResFinder acquired-resistance-gene database** (`resfinder-db`, installed via `apt-get install resfinder-db` — the real upstream CGE ResFinder gene set and its gene→phenotype mapping table, not a placeholder or invented database). Thresholds: ≥90% nucleotide identity, ≥60% reference-gene coverage (ResFinder's own default calling thresholds).

**Known limitation vs. AMRFinderPlus:** ResFinder covers acquired (horizontally-transferred) genes at the allele-sequence level; it does not do AMRFinderPlus's point-mutation/functional curation. This produced a **real, observed false positive** during testing: `blaSHV-99` (a narrow-spectrum, non-ESBL chromosomal SHV allele) is broadly mapped to cefepime in ResFinder's phenotype table, but genome `1226680.5`'s real lab result is Susceptible to cefepime. BLAST-identity matching alone cannot distinguish ESBL from non-ESBL SHV alleles the way AMRFinderPlus's curation would. This is exactly why every prediction carries the mandatory lab-recheck disclaimer.

### 1.3 Pipeline

```
split_genomes.py        kpneu_genomes.fna -> data/genomes/<genome_id>.fna (noise-stripped, per-genome)
detect_amr_genes.py     blastn each genome vs. ResFinder DB, collapse to best hit per genomic locus
build_feature_matrix.py results/amr_gene_hits.tsv -> results/gene_presence_absence_matrix.tsv
```

Result on the only real, in-scope genome (`1226680.5`): 4 genes detected — `OqxA`, `OqxB` (intrinsic efflux), `blaSHV-99` (intrinsic chromosomal beta-lactamase), `fosA6` (intrinsic fosfomycin resistance). No acquired aminoglycoside or trimethoprim/sulfamethoxazole gene was found, consistent with this genome's real Susceptible phenotype for gentamicin and trimethoprim/sulfamethoxazole.

The pipeline is genome-count-agnostic: it will process any number of genomes given more FASTA data.

---

## 2. Module 02 — Predictor

### 2.1 Deterministic gate

`gate.py` always returns **likely to fail** for ampicillin/amoxicillin in *K. pneumoniae*, independent of any model, because of the chromosomal, species-intrinsic SHV class A beta-lactamase. This is domain knowledge, not something fit from 25 samples. (Ampicillin itself is not one of the 3 statistically modeled antibiotics — it was excluded from the pilot's modeled set per the class-imbalance rule — but the gate is still implemented and always-on, independent of which antibiotics happen to be statistically modeled.)

### 2.2 Honest scope pivot: why there is no genome-feature LOOCV model

The original plan was: AMRFinderPlus/BLAST gene features for all 25 genomes → L2 logistic regression, LOOCV. §0.2 makes this impossible with real data — there is exactly **1** genome with real sequence, not 25. Fabricating gene-presence features for the other 24 to make the numbers work would be exactly the kind of "inflate the data to look good" behavior this project explicitly forbids, so it was not done.

Instead, to still deliver a real, non-fabricated, LOOCV-evaluated statistical baseline across all n=25 labeled genomes, `train_loocv.py` trains **L2-regularized logistic regression (`class_weight="balanced"`)** using each genome's **other real, lab-observed antibiotic phenotypes** (all antibiotics tested on that genome except the target) as features — a "cross-antibiotic co-resistance" baseline. This uses only real, unmodified data. Missing values are imputed with the **training-fold mean only** (recomputed inside each LOOCV split, no test-set leakage). Features with fewer than 15/25 real observations across the cohort are dropped.

**This is why LOOCV is used:** with n≤25 per antibiotic, any held-out split leaves single-digit test sets with unstable metrics; LOOCV uses every genome as a test point exactly once.

**Important:** this baseline is a retrospective-cohort statistical exercise, **not a genome-derived model**. It cannot be applied to a freshly uploaded genome that has no prior lab susceptibility results for other drugs — which is exactly the situation Module 03 faces. Module 03 therefore does **not** call this model; it only uses the deterministic gate and direct gene evidence from Module 01 (see §3).

**Known overfitting risk:** for cefepime, 20 candidate features are fit against 20 eligible samples. Metrics below should be read as optimistic, not as a reliable estimate of true generalization — flagged here rather than hidden.

### 2.3 LOOCV results (real data, n as noted per antibiotic)

| Antibiotic | n (R/S) | Balanced acc. | F1 | R recall | S recall | AUROC | PR-AUC | Brier | No-call rate (conf<0.70) |
|---|---|---|---|---|---|---|---|---|---|
| trimethoprim/sulfamethoxazole | 25 | 0.684 | 0.788 | 0.812 | 0.556 | 0.688 | 0.816 | 0.214 | 0.560 |
| gentamicin | 24 (1 Intermediate excluded) | 0.664 | 0.545 | 0.857 | 0.471 | 0.739 | 0.501 | 0.198 | 0.792 |
| cefepime | 20 (5 Intermediate excluded) | 0.804 | 0.941 | 0.941 | 0.667 | 0.941 | 0.991 | 0.102 | 0.150 |

Reliability plot: `results/reliability_plot.png`. Full per-genome predictions: `results/loocv_per_genome_predictions.tsv`. Full metrics JSON: `results/loocv_metrics.json`.

---

## 3. Module 03 — Decision Report app

`streamlit run src/module03_app/app.py`. Upload a *K. pneumoniae* genome FASTA:

1. Module 01 runs BLAST gene detection on the upload.
2. Module 02's deterministic gate is checked first (always wins if it fires).
3. If a detected gene maps to the antibiotic → **likely to fail**, evidence type `known_resistance_gene`.
4. Otherwise → **no-call**, evidence type `no_signal`. There is deliberately **no** `statistical_association_only` path wired into the live upload flow, because the only statistical model this pilot could validate (§2.2) needs other real antibiotic results for the *same* genome as input — data a fresh upload does not have. Building that path in anyway would mean quietly guessing; the honest behavior is to no-call.
5. Any prediction with confidence < 0.70 is force-converted to no-call.
6. Every result carries the fixed disclaimer: *"This is a decision-support estimate only. Results MUST be confirmed by standard laboratory antimicrobial susceptibility testing before any clinical or treatment decision is made."*

Verified: the app boots (`streamlit run`, health check `200 OK`) and the exact gene-detection + decision code paths it calls were exercised directly against the one real in-scope genome (`1226680.5`), producing gate/gene/no-call outputs consistent with §1.3. A full interactive browser click-through was not possible in this headless remote environment.

---

## 4. Responsibility Requirements

1. **Defensive by construction.** The pipeline only reads and classifies existing sequence (BLAST lookups against a reference gene database). Nothing in it can design, edit, or synthesize a genome or gene. No step accepts or produces synthesis-ready output.
2. **Honest generalization.** §0 and §2.2 document, in detail and up front, exactly why the model cannot generalize beyond this pilot: n=25 total, 1 usable genomic-feature row, small-sample overfitting risk explicitly flagged rather than hidden behind a good-looking accuracy number.
3. **Calibrated confidence & no-call.** Every prediction carries a confidence score; predictions below a 0.70 threshold, or with no gate/gene evidence, are forced to no-call rather than guessed. The real, measured no-call rate (15–79% depending on antibiotic) is reported, not hidden.
4. **Honest explanations.** Every call states which evidence type produced it (deterministic gate / known resistance gene / no signal) and, where applicable, which specific gene and BLAST identity/coverage drove the call — including the documented `blaSHV-99`/cefepime false-positive as a concrete example of the fallback method's limits (§1.2).
5. **Human oversight.** Every output, in the app and in this document, carries the fixed disclaimer that results must be reconfirmed by standard laboratory susceptibility testing before any clinical use.

---

## 5. Repository layout

```
kpneu_amr_metadata.tsv, kpneu_amr_phenotypes.tsv, kpneu_genomes.fna   # raw inputs (unmodified)
reference_db/                    # ResFinder gene DB + gene->phenotype mapping (apt package, real data)
src/module01_genome_reader/      # FASTA cleaning, BLAST gene detection, feature matrix
src/module02_predictor/          # deterministic gate, LOOCV baseline, upload-time predictor, reliability plot
src/module03_app/app.py          # Streamlit decision-report app
data/genomes/                    # cleaned per-genome FASTA (output of split_genomes.py)
data/blast_db/                   # makeblastdb output
results/                         # gene hits, feature matrix, LOOCV metrics/predictions, reliability plot
```

## 6. How to run

```bash
pip install scikit-learn pandas numpy scipy matplotlib streamlit biopython

python3 src/module01_genome_reader/split_genomes.py
python3 src/module01_genome_reader/detect_amr_genes.py
python3 src/module01_genome_reader/build_feature_matrix.py

cd src/module02_predictor
python3 train_loocv.py
python3 make_reliability_plot.py
cd ../..

streamlit run src/module03_app/app.py
```

## 7. Deploying (Streamlit Community Cloud, free)

**Already deployed:** https://hack-nation-6vjg298pfnsih9rwfi9fov.streamlit.app/

`requirements.txt` (Python deps) and `packages.txt` (`ncbi-blast+`, an apt dependency) are included at the repo root, and `app.py` auto-builds the BLAST database on first run — so no manual setup step is needed on the host.

1. Go to [share.streamlit.io](https://share.streamlit.io), sign in with GitHub, click "New app".
2. Repo: `sehyeonkim7565/hack-nation`, branch: `claude/kpneu-amr-prediction-f3llev`, main file path: `src/module03_app/app.py`.
3. Deploy. First boot takes a bit longer while the BLAST DB is built; subsequent restarts are fast.
