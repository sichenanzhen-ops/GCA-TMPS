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
}


def tier_from_score(entity: str, score: float, penalty: float) -> str:
    if entity == "KDM1A":
        return "Exploratory regulator"
    if entity == "PRDX2":
        return "Exploratory"
    if entity == "APEX1" and score >= 6 and penalty == 0:
        return "Tier 1"
    if score >= 4:
        return "Tier 2"
    return "Exploratory"


def rank_summary(scores: pd.DataFrame, score_col: str) -> dict:
    scores = scores.copy()
    scores["rank"] = scores[score_col].rank(method="min", ascending=False).astype(int)
    top = scores.sort_values([score_col, "Entity"], ascending=[False, True]).iloc[0]
    apex = scores[scores["Entity"] == "APEX1"].iloc[0]
    prdx = scores[scores["Entity"] == "PRDX2"].iloc[0]
    kdm = scores[scores["Entity"] == "KDM1A"].iloc[0]
    return {
        "Top-ranked entity": top["Entity"],
        "APEX1 score": round(float(apex[score_col]), 3),
        "APEX1 rank": int(apex["rank"]),
        "APEX1 tier": tier_from_score("APEX1", float(apex[score_col]), float(apex[DOMAIN_COLS["penalty"]])),
        "PRDX2 score": round(float(prdx[score_col]), 3),
        "PRDX2 rank": int(prdx["rank"]),
        "PRDX2 tier": tier_from_score("PRDX2", float(prdx[score_col]), float(prdx[DOMAIN_COLS["penalty"]])),
        "KDM1A score": round(float(kdm[score_col]), 3),
        "KDM1A rank": int(kdm["rank"]),
        "KDM1A tier": tier_from_score("KDM1A", float(kdm[score_col]), float(kdm[DOMAIN_COLS["penalty"]])),
    }


def ablation(df: pd.DataFrame) -> pd.DataFrame:
    scenarios = {
        "Full TMPS": {},
        "No protein-level evidence": {"protein": 0},
        "No SMR/expression-context evidence": {"coloc_expr": 0},
        "No inflammatory/redox context": {"inflammatory": 0},
        "No regulatory-context evidence": {"regulatory": 0},
        "No traceability penalty": {"penalty": 0},
    }
    rows = []
    for scenario, overrides in scenarios.items():
        temp = df.copy()
        for key, value in overrides.items():
            temp[DOMAIN_COLS[key]] = value
        temp["Scenario score"] = (
            temp[DOMAIN_COLS["genetic"]]
            + temp[DOMAIN_COLS["coloc_expr"]]
            + temp[DOMAIN_COLS["protein"]]
            + temp[DOMAIN_COLS["inflammatory"]]
            + temp[DOMAIN_COLS["regulatory"]]
            - temp[DOMAIN_COLS["penalty"]]
        )
        summary = rank_summary(temp, "Scenario score")
        summary.update(
            {
                "Analysis block": "Ablation",
                "Scenario": scenario,
                "Weighting/removed domain": "; ".join(overrides) if overrides else "none",
                "Interpretation": "APEX1 remained top-ranked" if summary["APEX1 rank"] == 1 else "APEX1 rank changed",
            }
        )
        rows.append(summary)
    return pd.DataFrame(rows)


def parameter_robustness(df: pd.DataFrame) -> pd.DataFrame:
    schemes = {
        "Primary weights": {"genetic": 1.0, "coloc_expr": 1.0, "protein": 1.0, "inflammatory": 1.0, "regulatory": 1.0, "penalty": 1.0},
        "Genetic-heavy weights": {"genetic": 1.25, "coloc_expr": 0.85, "protein": 0.85, "inflammatory": 0.75, "regulatory": 0.75, "penalty": 1.0},
        "Colocalization/expression-heavy weights": {"genetic": 0.85, "coloc_expr": 1.25, "protein": 0.85, "inflammatory": 0.75, "regulatory": 0.75, "penalty": 1.0},
        "Protein-heavy weights": {"genetic": 0.90, "coloc_expr": 0.90, "protein": 1.50, "inflammatory": 0.75, "regulatory": 0.75, "penalty": 1.0},
        "Penalty-heavy weights": {"genetic": 1.0, "coloc_expr": 1.0, "protein": 1.0, "inflammatory": 1.0, "regulatory": 1.0, "penalty": 1.5},
        "Equal-normalized domains": {"genetic": 1 / 3, "coloc_expr": 1 / 3, "protein": 1 / 2, "inflammatory": 1.0, "regulatory": 1.0, "penalty": 1.0},
    }
    rows = []
    for scenario, weights in schemes.items():
        temp = df.copy()
        temp["Scenario score"] = (
            temp[DOMAIN_COLS["genetic"]] * weights["genetic"]
            + temp[DOMAIN_COLS["coloc_expr"]] * weights["coloc_expr"]
            + temp[DOMAIN_COLS["protein"]] * weights["protein"]
            + temp[DOMAIN_COLS["inflammatory"]] * weights["inflammatory"]
            + temp[DOMAIN_COLS["regulatory"]] * weights["regulatory"]
            - temp[DOMAIN_COLS["penalty"]] * weights["penalty"]
        )
        summary = rank_summary(temp, "Scenario score")
        summary.update(
            {
                "Analysis block": "Parameter robustness",
                "Scenario": scenario,
                "Weighting/removed domain": "; ".join(f"{key}={round(value, 3)}" for key, value in weights.items()),
                "Interpretation": "APEX1 remained top-ranked; PRDX2 remained exploratory" if summary["APEX1 rank"] == 1 and summary["PRDX2 tier"] == "Exploratory" else "Rank/tier changed; interpret cautiously",
            }
        )
        rows.append(summary)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run TMPS ablation and parameter-robustness analyses.")
    parser.add_argument("--input", default="outputs/TMPS_recomputed_scores.csv")
    parser.add_argument("--out", default="outputs/TMPS_recomputed_ablation_robustness.csv")
    args = parser.parse_args()
    df = pd.read_csv(args.input)
    for col in DOMAIN_COLS.values():
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    result = pd.concat([ablation(df), parameter_robustness(df)], ignore_index=True)
    result.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
