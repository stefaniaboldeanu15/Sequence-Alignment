from itertools import combinations

from src.pairwise.needleman_wunsch import needleman_wunsch
from src.utils.fasta import read_fasta


def short_name(header: str) -> str:
    parts = header.split("|")
    if len(parts) >= 3:
        return parts[2].split()[0]
    return header.split()[0]


def run_dataset(fasta_path: str, match_score=1, mismatch_score=-1, gap_penalty=-1):
    records = read_fasta(fasta_path)

    print(f"\nLoaded {len(records)} sequences from {fasta_path}\n")

    for rec1, rec2 in combinations(records, 2):
        result = needleman_wunsch(
            rec1.sequence,
            rec2.sequence,
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_penalty=gap_penalty,
        )

        print("=" * 80)
        print(f"{short_name(rec1.header)}  vs  {short_name(rec2.header)}")
        print(f"Score: {result.score}")
        print(result.aligned_seq1)
        print(result.aligned_seq2)
        print(f"Matches: {result.matches}")
        print(f"Mismatches: {result.mismatches}")
        print(f"Gaps: {result.gaps}")


if __name__ == "__main__":
    run_dataset("data/raw/globins.fasta")