#!/usr/bin/env python3

from pathlib import Path
import argparse


def read_lengths(fasta_path: Path) -> list[int]:
    lengths: list[int] = []
    cur_len = 0

    with fasta_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_len:
                    lengths.append(cur_len)
                cur_len = 0
            else:
                cur_len += len(line)

    if cur_len:
        lengths.append(cur_len)

    return lengths


def compute_n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    cutoff = sum(ordered) / 2.0
    running = 0
    for seq_len in ordered:
        running += seq_len
        if running >= cutoff:
            return seq_len
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute FASTA contiguity stats.")
    parser.add_argument("--fasta", required=True, help="Input FASTA file.")
    parser.add_argument("--n50-out", default="n50.txt", help="Output path for N50 value.")
    parser.add_argument(
        "--num-large-out",
        default="num_large_contigs.txt",
        help="Output path for >=1Mb contig count.",
    )
    parser.add_argument(
        "--large-threshold",
        type=int,
        default=1_000_000,
        help="Minimum contig length to count as large.",
    )
    args = parser.parse_args()

    fasta = Path(args.fasta)
    lengths = read_lengths(fasta)
    n50 = compute_n50(lengths)
    num_large = sum(1 for seq_len in lengths if seq_len >= args.large_threshold)

    Path(args.n50_out).write_text(f"{n50}\n")
    Path(args.num_large_out).write_text(f"{num_large}\n")


if __name__ == "__main__":
    main()
