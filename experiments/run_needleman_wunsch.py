from itertools import combinations
from pathlib import Path

from src.pairwise.needleman_wunsch import needleman_wunsch
from src.utils.fasta import read_fasta


def short_name(header: str) -> str:
    parts = header.split("|")
    if len(parts) >= 3:
        return parts[2].split()[0]
    return header.split()[0]


def run_dataset_to_txt(
    fasta_path: str,
    output_path: str,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_penalty: int = -1,
) -> None:
    records = read_fasta(fasta_path)

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    print(f"\nReading FASTA file: {fasta_path}")
    print(f"Loaded {len(records)} sequences")
    print(f"Writing results to: {output_path}")

    with open(output_path, "w", encoding="utf-8") as out:
        out.write(f"Dataset: {fasta_path}\n")
        out.write(f"Sequences loaded: {len(records)}\n")
        out.write(f"Scoring: match={match_score}, mismatch={mismatch_score}, gap={gap_penalty}\n\n")

        for rec in records:
            out.write(f"{short_name(rec.header)} length = {len(rec.sequence)}\n")
        out.write("\n")

        pairs = list(combinations(records, 2))
        out.write(f"Total pairwise alignments: {len(pairs)}\n\n")

        for idx, (rec1, rec2) in enumerate(pairs, start=1):
            name1 = short_name(rec1.header)
            name2 = short_name(rec2.header)

            print(f"Running {idx}/{len(pairs)}: {name1} vs {name2}")

            result = needleman_wunsch(
                rec1.sequence,
                rec2.sequence,
                match_score=match_score,
                mismatch_score=mismatch_score,
                gap_penalty=gap_penalty,
            )

            out.write("=" * 80 + "\n")
            out.write(f"{name1}  vs  {name2}\n")
            out.write(f"Score: {result.score}\n")
            out.write(result.aligned_seq1 + "\n")
            out.write(result.aligned_seq2 + "\n")
            out.write(f"Matches: {result.matches}\n")
            out.write(f"Mismatches: {result.mismatches}\n")
            out.write(f"Gaps: {result.gaps}\n\n")

    print(f"Finished: {output_path}")


def run_all_datasets() -> None:
    datasets = [
        ("data/raw/globins.fasta", "results/globins_needleman_wunsch.txt"),
        ("data/raw/serine_proteases.fasta", "results/serine_proteases_needleman_wunsch.txt"),
        ("data/synthetic/synthetic_controls.fasta", "results/synthetic_controls_needleman_wunsch.txt"),
    ]

    available = []
    for fasta_path, output_path in datasets:
        if Path(fasta_path).exists():
            available.append((fasta_path, output_path))
        else:
            print(f"Skipping missing dataset: {fasta_path}")

    if not available:
        print("No dataset files found.")
        return

    for fasta_path, output_path in available:
        run_dataset_to_txt(
            fasta_path=fasta_path,
            output_path=output_path,
            match_score=1,
            mismatch_score=-1,
            gap_penalty=-1,
        )