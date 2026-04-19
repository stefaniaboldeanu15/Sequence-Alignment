from dataclasses import dataclass


@dataclass
class FastaRecord:
    header: str
    sequence: str


def read_fasta(path: str) -> list[FastaRecord]:
    records = []
    header = None
    seq_lines = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append(FastaRecord(header=header, sequence="".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)

    if header is not None:
        records.append(FastaRecord(header=header, sequence="".join(seq_lines)))

    return records