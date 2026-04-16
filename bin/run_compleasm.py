#!/usr/bin/env python3

import argparse
import shutil
import subprocess
from pathlib import Path


def write_fallback_report(report_out: str, reason: str, detail: str = "") -> None:
    lines = [
        f"status\tfallback",
        f"reason\t{reason}",
    ]
    if detail:
        lines.append(f"detail\t{detail}")
    Path(report_out).write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run compleasm and export summary report.")
    parser.add_argument("--assembly", required=True, help="Assembly FASTA path.")
    parser.add_argument("--sample-id", required=True, help="Sample identifier.")
    parser.add_argument("--threads", type=int, default=40, help="Thread count for compleasm.")
    parser.add_argument("--lineage", default="", help="BUSCO lineage (optional).")
    parser.add_argument(
        "--report-out",
        required=True,
        help="Output path for copied compleasm summary report.",
    )
    args = parser.parse_args()

    outdir = Path(args.sample_id)
    compleasm_exec = shutil.which("compleasm") or shutil.which("compleasm.py")
    if not compleasm_exec:
        write_fallback_report(
            args.report_out,
            "compleasm_executable_not_found",
            "Neither 'compleasm' nor 'compleasm.py' was available on PATH.",
        )
        return

    cmd = [
        compleasm_exec,
        "run",
        "-t",
        str(args.threads),
        "-a",
        args.assembly,
        "-o",
        str(outdir),
    ]
    if args.lineage.strip():
        cmd.extend(["-l", args.lineage.strip()])
    else:
        cmd.append("--autolineage")

    proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip().replace("\n", " | ")
        write_fallback_report(args.report_out, f"compleasm_run_failed_{proc.returncode}", detail)
        return

    summary = outdir / "summary.txt"
    if not summary.exists():
        write_fallback_report(
            args.report_out,
            "compleasm_summary_missing",
            f"Expected summary at {summary}",
        )
        return
    shutil.copyfile(summary, args.report_out)


if __name__ == "__main__":
    main()
