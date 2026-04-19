# Sequence Alignment

This repository contains Python implementations of exact pairwise sequence alignment methods developed for the Sequence Alignment assignment.

## Implemented methods

- Needleman–Wunsch (global alignment)
- Smith–Waterman (local alignment)
- Gotoh (global alignment with affine gap penalties)
- Hirschberg (space-efficient global alignment)

## Project structure

```text
src/
  pairwise/      # alignment algorithm implementations
  utils/         # FASTA parser and helper functions
tests/           # unit tests
experiments/     # dataset runners
data/            # input FASTA files
results/         # exported alignment outputs
main.py          # entry point
```

## Requirements

- Python 3.11+
- pytest

```

## Run experiments

Edit `main.py` so it calls the method you want to run, then execute:

```bash
python main.py
```

The output files are written to the `results/` folder.

## Datasets

Input FASTA files are stored in:

```text
data/raw/
data/synthetic/
```

Current datasets:
- globins
- serine proteases
- synthetic controls

## Scoring

### Needleman–Wunsch / Smith–Waterman / Hirschberg
- match = +1
- mismatch = -1
- gap = -1

### Gotoh
- match = +1
- mismatch = -1
- gap opening = -2
- gap extension = -1

## Notes

These implementations were created for correctness testing, comparison of alignment methods, and reproducible experiments on biological and synthetic datasets.