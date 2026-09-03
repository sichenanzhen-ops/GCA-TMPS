from __future__ import annotations

import argparse
import pandas as pd


DOMAIN_COLS = {
    "genetic": "Genetic association support (0-3)",
    "coloc_expr": "Colocalization/transcriptomic support (0-3)",
    "protein": "Protein-level support (0-2)",
    "inflammatory": "Inflammatory/oxidative-stress context (0-1)",
    "regulatory": "Regulatory-context support (0-1)",
    "penalty": "Penalty",
    "tmps": "TMPS total score",
}


def prepare_scores(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in DOMAIN_COLS.values():
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    if DOMAIN_COLS["tmps"] not in df.columns:
        df[DOMAIN_COLS["tmps"]] = (
            df[DOMAIN_COLS["genetic"]]
            + df[DOMAIN_COLS["coloc_expr"]]
            + df[DOMAIN_COLS["protein"]]
            + df[DOMAIN_COLS["inflammatory"]]
            + df[DOMAIN_COLS["regulatory"]]
            - df[DOMAIN_COLS["penalty"]]
        )
    return df


def add_rank(rows: list[dict], df: pd.DataFrame, scenario: str, score_col: str, interpretation: str) -> None:
    tmp = df[["Entity", score_col]].copy()
    tmp["rank"] = tmp[score_col].rank(method="min", ascending=False).astype(int)
    tmp = tmp.sort_values([score_col, "Entity"], ascending=[False, True])
    for _, row in tmp.iterrows():
        rows.append(
            {
                "Scenario": scenario,
                "Entity": row["Entity"],
                "Score used for ranking": round(float(row[score_col]), 3),
                "Rank": int(row["rank"]),
                "Interpretation boundary": interpretation,
            }
        )


def benchmark_rankings(df: pd.DataFrame) -> pd.DataFrame:
    df = prepare_scores(df)
    df["Simple unweighted evidence sum"] = (
        df[DOMAIN_COLS["genetic"]]
        + df[DOMAIN_COLS["coloc_expr"]]
        + df[DOMAIN_COLS["protein"]]
        + df[DOMAIN_COLS["inflammatory"]]
        + df[DOMAIN_COLS["regulatory"]]
    )
    df["Positive evidence-domain count"] = (
        (df[DOMAIN_COLS["genetic"]] > 0).astype(int)
        + (df[DOMAIN_COLS["coloc_expr"]] > 0).astype(int)
        + (df[DOMAIN_COLS["protein"]] > 0).astype(int)
        + (df[DOMAIN_COLS["inflammatory"]] > 0).astype(int)
        + (df[DOMAIN_COLS["regulatory"]] > 0).astype(int)
    )

    rows: list[dict] = []
    add_rank(
        rows,
        df,
        "Full TMPS",
        DOMAIN_COLS["tmps"],
        "Primary transparent multi-domain score with traceability penalties.",
    )
    add_rank(
        rows,
        df,
        "Simple unweighted evidence sum",
        "Simple unweighted evidence sum",
        "Naive evidence aggregation without traceability penalties.",
    )
    add_rank(
        rows,
        df,
        "Positive evidence-domain count",
        "Positive evidence-domain count",
        "Counts whether each evidence domain is present, regardless of strength.",
    )
    add_rank(
        rows,
        df,
        "Genetic-evidence-only ranking",
        DOMAIN_COLS["genetic"],
        "Approximates a single-domain MR-first ranking and ignores contextual evidence.",
    )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare TMPS ranking with simpler baseline ranking rules.")
    parser.add_argument("--input", default="data/tmps_input_evidence.csv")
    parser.add_argument("--out", default="outputs/TMPS_benchmark_rankings.csv")
    args = parser.parse_args()
    df = pd.read_csv(args.input)
    benchmark_rankings(df).to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
