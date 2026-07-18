"""
Module 01 step A: parse kpneu_genomes.fna and split it into one clean FASTA
file per genome_id.

Why this script exists (data-quality note, not hypothetical):
kpneu_genomes.fna is not a clean multi-FASTA. It contains leftover Windows
terminal transcript text from the PATRIC CLI export process
(`C:\\Program Files\\PATRIC\\cli\\bin>p3-genome-fasta --contig <genome_id>`
command echoes, and repeated EUC-KR "batch job confirm (Y/N)?" prompts).
This script uses those `--contig <genome_id>` markers to attribute each
FASTA record to its true genome_id, and silently drops any line that is
not a valid FASTA header or a valid nucleotide sequence line.

It writes one FASTA file per genome_id that has *actual* sequence data.
Genome IDs referenced by a `--contig` command but never followed by any
sequence (a failed/truncated fetch) are reported, not silently skipped.
"""
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
INPUT_FASTA = REPO_ROOT / "kpneu_genomes.fna"
OUTPUT_DIR = REPO_ROOT / "data" / "genomes"

SEQ_LINE_RE = re.compile(r"^[ACGTNacgtn]+$")
CONTIG_CMD_RE = re.compile(r"p3-genome-fasta --contig (\S+)")


def split_genomes(input_fasta: Path, output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    raw = input_fasta.read_bytes().decode("utf-8", errors="replace")
    lines = raw.split("\n")

    current_genome_id = None
    current_header = None
    buffers: dict[str, list[str]] = {}
    header_seen_for_genome: dict[str, list[str]] = {}
    genome_ids_commanded: list[str] = []

    for line in lines:
        cmd_match = CONTIG_CMD_RE.search(line)
        if cmd_match:
            current_genome_id = cmd_match.group(1)
            genome_ids_commanded.append(current_genome_id)
            continue

        if line.startswith(">"):
            current_header = line.strip()
            if current_genome_id is None:
                # No preceding command marker seen yet; skip (should not
                # happen given the observed file structure, but fail safe
                # rather than mis-attribute sequence to the wrong genome).
                continue
            buffers.setdefault(current_genome_id, []).append(current_header + "\n")
            header_seen_for_genome.setdefault(current_genome_id, []).append(current_header)
            continue

        if SEQ_LINE_RE.match(line):
            if current_genome_id is not None and current_genome_id in buffers:
                buffers[current_genome_id].append(line + "\n")
            continue

        # Anything else (blank lines, CLI transcript noise, garbled
        # confirmation prompts) is intentionally dropped.

    written = {}
    for genome_id, chunks in buffers.items():
        out_path = output_dir / f"{genome_id}.fna"
        out_path.write_text("".join(chunks))
        n_bases = sum(len(c.strip()) for c in chunks if not c.startswith(">"))
        n_contigs = len(header_seen_for_genome[genome_id])
        written[genome_id] = {"contigs": n_contigs, "bases": n_bases, "path": str(out_path)}

    empty_commands = [g for g in genome_ids_commanded if g not in buffers]
    return {
        "written": written,
        "commanded_but_no_sequence": sorted(set(empty_commands)),
    }


if __name__ == "__main__":
    result = split_genomes(INPUT_FASTA, OUTPUT_DIR)
    print("Genomes with real sequence recovered:")
    for gid, info in sorted(result["written"].items()):
        print(f"  {gid}: {info['contigs']} contigs, {info['bases']} bp -> {info['path']}")
    if result["commanded_but_no_sequence"]:
        print("\nWARNING: fetch was attempted but NO sequence was recovered for:")
        for gid in result["commanded_but_no_sequence"]:
            print(f"  {gid}  (source file transcript is truncated/corrupted here)")
    else:
        print("\nNo commanded-but-empty genome IDs found.")
