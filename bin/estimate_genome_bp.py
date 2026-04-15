#!/usr/bin/env python3

import argparse
from pathlib import Path


def fasta_len(path_str: str) -> int | None:
    if not path_str:
        return None
    path = Path(path_str)
    if not path.exists():
        return None
    total = 0
    with path.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total if total > 0 else None


def parse_user_size(raw: str) -> int | None:
    if not raw:
        return None
    raw = raw.strip().lower()
    if not raw:
        return None

    multipliers = {"g": 1_000_000_000, "m": 1_000_000, "k": 1_000}
    try:
        if raw[-1] in multipliers:
            return int(float(raw[:-1]) * multipliers[raw[-1]])

        value = float(raw)
        # Historical ONTeater docs treated bare values as gigabases.
        if value < 1000:
            return int(value * 1_000_000_000)
        return int(value)
    except ValueError:
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Estimate genome size in base pairs.")
    parser.add_argument("--user-size", default="", help="User-provided genome size (e.g. 3g, 250m, 1.2).")
    parser.add_argument("--flye-asm", default="", help="Path to Flye assembly FASTA.")
    parser.add_argument("--assembly", required=True, help="Path to assembly FASTA under assessment.")
    args = parser.parse_args()

    user_size = parse_user_size(args.user_size)
    flye_size = fasta_len(args.flye_asm)
    asm_size = fasta_len(args.assembly)

    estimated = user_size or flye_size or asm_size or 0
    print(estimated)


if __name__ == "__main__":
    main()
