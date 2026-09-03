from __future__ import annotations

import argparse
import pandas as pd


DOMAIN_COLUMNS = [
    "Genetic association support (0-3)",
    "Colocalization/transcriptomic support (0-3)",
    "Protein-level support (0-2)",
    "Inflammatory/oxidative-stress context (0-1)",
    "Regulatory-context support (0-1)",
    "Penalty",
]


def evidence_tier(entity: str, total_score: float, penalty: float) -> str:
    if entity == "KDM1A":
        return "Exploratory regulator"
    if entity == "PRDX2":
        return "Exploratory"
    if entity == "APEX1" and total_score >= 6 and penalty == 0:
        return "Tier 1"
    if total_score >= 4:
        return "Tier 2"
    return "Exploratory"


def recompute_scores(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in DOMAIN_COLUMNS:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    df["TMPS total score"] = (
        df["Genetic association support (0-3)"]
        + df["Colocalization/transcriptomic support (0-3)"]
        + df["Protein-level support (0-2)"]
        + df["Inflammatory/oxidative-stress context (0-1)"]
        + df["Regulatory-context support (0-1)"]
        - df["Penalty"]
    ).round(3)
    df["Evidence tier"] = [
        evidence_tier(entity, score, penalty)
        for entity, score, penalty in zip(df["Entity"], df["TMPS total score"], df["Penalty"])
    ]
    return df.sort_values(["TMPS total score", "Entity"], ascending=[False, True])


def main() -> None:
    parser = argparse.ArgumentParser(description="Recompute TMPS candidate ranking.")
    parser.add_argument("--input", default="data/tmps_input_evidence.csv")
    parser.add_argument("--out", default="outputs/TMPS_recomputed_scores.csv")
    args = parser.parse_args()
    df = pd.read_csv(args.input)
    recompute_scores(df).to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
