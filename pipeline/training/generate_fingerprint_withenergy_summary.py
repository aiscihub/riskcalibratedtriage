#!/usr/bin/env python3
"""Build a fingerprint summary CSV for energy-folder trajectory outputs."""

from pathlib import Path

try:
    from .generate_fingerprint_summary import build_fingerprint_summary, infer_ligand_name
except ImportError:
    from generate_fingerprint_summary import build_fingerprint_summary, infer_ligand_name


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Summarize fingerprints from fingerprints_ca_energy_4 folders."
    )
    parser.add_argument(
        "--top-dir",
        type=Path,
        required=True,
        help="Simulation root containing protein/simulation_explicit/pocket*/replica_1.",
    )
    parser.add_argument(
        "--labels",
        type=Path,
        required=True,
        help="Endpoint label CSV.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/fingerprint_summary_with_components_even_ac_4d_drift.csv"),
        help="Output summary CSV.",
    )
    parser.add_argument(
        "--ligand-name",
        default=None,
        help="Ligand name used in labels and trajectory filenames. Defaults from --top-dir.",
    )
    parser.add_argument(
        "--fingerprint-dir-name",
        default="fingerprints_ca_energy_4",
        help="Fingerprint directory under each replica folder.",
    )
    args = parser.parse_args()

    ligand_name = args.ligand_name or infer_ligand_name(args.top_dir)
    df = build_fingerprint_summary(
        top_dir=args.top_dir,
        label_file=args.labels,
        output=args.output,
        ligand_name=ligand_name,
        fingerprint_dir_name=args.fingerprint_dir_name,
    )
    print(f"Saved {args.output} ({len(df)} rows)")


if __name__ == "__main__":
    main()
