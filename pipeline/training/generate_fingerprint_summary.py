#!/usr/bin/env python3
"""Build a model-ready fingerprint summary CSV from local MD outputs."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


FEATURE_NAMES = [
    "slope",
    "rmsd_var",
    "mean_disp",
    "var_disp",
]


def load_fingerprints(fp_dir: Path) -> dict[str, np.ndarray] | None:
    paths = sorted(fp_dir.glob("fingerprint_*ns.npy"))
    if not paths:
        return None

    return {
        path.name.split("_")[-1].replace("ns.npy", ""): np.load(path)
        for path in paths
    }


def infer_ligand_name(root: Path) -> str:
    root_text = str(root).lower()
    if "beau" in root_text:
        return "Beauvericin"
    if "vera" in root_text:
        return "Verapamil"
    return "Milbemycin"


def build_fingerprint_summary(
    top_dir: Path,
    label_file: Path,
    output: Path,
    ligand_name: str,
    fingerprint_dir_name: str = "fingerprints_ca_4",
) -> pd.DataFrame:
    labels_df = pd.read_csv(label_file)

    labels_key = {}
    for _, r in labels_df.iterrows():
        key = (r["protein"], r["ligand_name"], r["pocket_id"], r["replica"])
        labels_key[key] = {
            "rmsd_20ns": r["rmsd_late_20ns"],
            "f_contact_20ns": r.get("f_contact_20ns", np.nan),
            "ligand_drift": r["drift"],
            "label_unstable": r["label_unstable"],
        }

    rows = []
    for protein_path in sorted(path for path in top_dir.iterdir() if path.is_dir()):
        protein = protein_path.name
        pocket_dirs = sorted((protein_path / "simulation_explicit").glob("pocket*"))

        for pocket_path in pocket_dirs:
            pocket = pocket_path.name
            replica = "replica_1"
            fp_dir = pocket_path / replica / fingerprint_dir_name

            if not fp_dir.is_dir():
                continue

            fps = load_fingerprints(fp_dir)
            if fps is None:
                continue

            label_info = labels_key.get(
                (protein, ligand_name, pocket, replica),
                {
                    "rmsd_20ns": np.nan,
                    "ligand_drift": np.nan,
                    "f_contact_20ns": np.nan,
                    "label_unstable": np.nan,
                },
            )

            for t in sorted(fps.keys(), key=float):
                vec = np.asarray(fps[t], dtype=float)
                row = {
                    "protein": protein,
                    "ligand_name": ligand_name,
                    "pocket": pocket,
                    "replica": replica,
                    "time_ns": float(t),
                    "mean": float(vec.mean()) if vec.size else np.nan,
                    "std": float(vec.std()) if vec.size else np.nan,
                    "n_features": len(FEATURE_NAMES),
                    "rmsd_20ns": label_info["rmsd_20ns"],
                    "ligand_drift": label_info["ligand_drift"],
                    "f_contact_20ns": label_info["f_contact_20ns"],
                    "label_unstable": label_info["label_unstable"],
                }

                for i, name in enumerate(FEATURE_NAMES):
                    row[name] = float(vec[i]) if i < vec.size else np.nan

                rows.append(row)

    df = pd.DataFrame(rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    return df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize per-horizon fingerprint .npy files into a CSV."
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
        default=Path("label_drift_20ns.csv"),
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
        default="fingerprints_ca_4",
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
