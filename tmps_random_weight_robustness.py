from __future__ import annotations

import argparse
import numpy as np
import pandas as pd


DOMAIN_COLS = {
    "genetic": "Genetic association support (0-3)",
    "coloc_expr": "Colocalization/transcriptomic support (0-3)",
    "protein": "Protein-level support (0-2)",
    "inflammatory": "Inflammatory/oxidative-stress context (0-1)",
    "regulatory": "Regulatory-context support (0-1)",
    "penalty": "Penalty",
}


def score_with_weights(df: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
    return (
        df[DOMAIN_COLS["genetic"]] * weights["genetic"]
        + df[DOMAIN_COLS["coloc_expr"]] * weights["coloc_expr"]
        + df[DOMAIN_COLS["protein"]] * weights["protein"]
        + df[DOMAIN_COLS["inflammatory"]] * weights["inflammatory"]
        + df[DOMAIN_COLS["regulatory"]] * weights["regulatory"]
        - df[DOMAIN_COLS["penalty"]] * weights["penalty"]
    )


def run_random_weight_robustness(
    df: pd.DataFrame,
    iterations: int,
    low: float,
    high: float,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = df.copy()
    for col in DOMAIN_COLS.values():
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    rng = np.random.default_rng(seed)
    rows = []
    for i in range(1, iterations + 1):
        weights = {
            key: float(rng.uniform(low, high))
            for key in ["genetic", "coloc_expr", "protein", "inflammatory", "regulatory", "penalty"]
        }
        scores = score_with_weights(df, weights)
        ranks = scores.rank(method="min", ascending=False).astype(int)
        top_score = scores.max()
        top_entities = sorted(df.loc[np.isclose(scores, top_score), "Entity"].tolist())
        for idx, entity in zip(df.index, df["Entity"]):
            rows.append(
                {
                    "Iteration": i,
                    "Entity": entity,
                    "Random-weight score": round(float(scores.loc[idx]), 4),
                    "Rank": int(ranks.loc[idx]),
                    "Top-ranked entities": "; ".join(top_entities),
                    "Sole top-ranked APEX1": entity == "APEX1" and ranks.loc[idx] == 1 and top_entities == ["APEX1"],
                    **{f"weight_{key}": round(value, 4) for key, value in weights.items()},
                }
            )

    detail = pd.DataFrame(rows)
    summary = (
        detail.groupby("Entity")
        .agg(
            iterations=("Iteration", "count"),
            median_rank=("Rank", "median"),
            best_rank=("Rank", "min"),
            worst_rank=("Rank", "max"),
            rank_1_frequency=("Rank", lambda x: float((x == 1).mean())),
            sole_APEX1_top_frequency=("Sole top-ranked APEX1", "mean"),
            median_score=("Random-weight score", "median"),
        )
        .reset_index()
    )
    for col in ["rank_1_frequency", "sole_APEX1_top_frequency"]:
        summary[col] = (summary[col] * 100).round(1)
    summary["median_score"] = summary["median_score"].round(3)
    return detail, summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Monte Carlo random-weight robustness analysis for TMPS.")
    parser.add_argument("--input", default="data/tmps_input_evidence.csv")
    parser.add_argument("--detail-out", default="outputs/TMPS_random_weight_detail.csv")
    parser.add_argument("--summary-out", default="outputs/TMPS_random_weight_summary.csv")
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--low", type=float, default=0.75)
    parser.add_argument("--high", type=float, default=1.25)
    parser.add_argument("--seed", type=int, default=20260903)
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    detail, summary = run_random_weight_robustness(df, args.iterations, args.low, args.high, args.seed)
    detail.to_csv(args.detail_out, index=False)
    summary.to_csv(args.summary_out, index=False)


if __name__ == "__main__":
    main()
