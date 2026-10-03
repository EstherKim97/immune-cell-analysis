"""Part 3: compare immune-cell frequencies by miraclib response.

The pooled sample comparison answers the assessment question. Subject means and
baseline samples provide independent-subject sensitivity checks.
"""

import math
import sqlite3
from contextlib import closing
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu

from frequency_table import DB_PATH, POPULATIONS, get_frequency_table


OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"
POPULATION_ORDER = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")


def load_cohort(connection: sqlite3.Connection) -> pd.DataFrame:
    """Join Part 2 frequencies to the exact melanoma/miraclib/PBMC cohort."""
    metadata = pd.read_sql_query(
        """
        SELECT s.sample_id AS sample, s.subject_id, s.sample_type,
               s.time_from_treatment_start, u.project_id, u.condition,
               u.treatment, u.response
        FROM samples AS s
        JOIN subjects AS u ON u.subject_id = s.subject_id
        WHERE u.condition = 'melanoma' AND u.treatment = 'miraclib'
          AND s.sample_type = 'PBMC' AND u.response IN ('yes', 'no')
        ORDER BY s.sample_id
        """,
        connection,
    )
    if metadata.empty:
        raise ValueError("No samples match the Part 3 cohort")
    if metadata["sample"].duplicated().any():
        raise ValueError("Duplicate sample IDs in cohort metadata")
    if (metadata.groupby("subject_id")["response"].nunique() != 1).any():
        raise ValueError("Response changes within a subject")

    frequencies = pd.DataFrame.from_records(
        get_frequency_table(connection),
        columns=("sample", "total_count", "population", "count", "percentage"),
    )
    cohort = metadata.merge(frequencies, on="sample", how="left", validate="one_to_many")
    if len(cohort) != len(metadata) * len(POPULATIONS):
        raise ValueError("Cohort is missing frequency rows")
    if cohort["percentage"].isna().any():
        raise ValueError("Cohort contains missing percentages")
    if not (cohort.condition.eq("melanoma").all() and cohort.treatment.eq("miraclib").all()
            and cohort.sample_type.eq("PBMC").all() and cohort.response.isin(("yes", "no")).all()):
        raise ValueError("Cohort filter includes ineligible observations")
    if (cohort.groupby("sample")["population"].nunique() != len(POPULATIONS)).any():
        raise ValueError("A cohort sample is missing a population")
    return cohort


def summarize_cohort(cohort: pd.DataFrame) -> dict:
    """Return cohort counts and repeated-measures structure for review."""
    samples = cohort.drop_duplicates("sample")
    subjects = samples.drop_duplicates("subject_id")
    return {
        "subjects": len(subjects),
        "samples": len(samples),
        "responder_subjects": int((subjects.response == "yes").sum()),
        "nonresponder_subjects": int((subjects.response == "no").sum()),
        "responder_samples": int((samples.response == "yes").sum()),
        "nonresponder_samples": int((samples.response == "no").sum()),
        "samples_per_subject": samples.groupby("subject_id").size().value_counts().sort_index().to_dict(),
        "timepoints": samples.time_from_treatment_start.value_counts().sort_index().to_dict(),
        "population_observations": cohort.population.value_counts().to_dict(),
    }


def fdr_bh(p_values: list[float]) -> list[float]:
    """Benjamini–Hochberg adjusted p-values across the five populations."""
    adjusted = [0.0] * len(p_values)
    running_minimum = 1.0
    for rank, index in reversed(list(enumerate(sorted(range(len(p_values)), key=p_values.__getitem__), 1))):
        running_minimum = min(running_minimum, p_values[index] * len(p_values) / rank)
        adjusted[index] = running_minimum
    return adjusted


def compare_groups(data: pd.DataFrame) -> pd.DataFrame:
    """Distribution-free two-group tests; positive effect means responder dominance."""
    results = []
    for population in POPULATION_ORDER:
        subset = data.loc[data.population == population]
        responders = subset.loc[subset.response == "yes", "percentage"].to_numpy(dtype=float)
        nonresponders = subset.loc[subset.response == "no", "percentage"].to_numpy(dtype=float)
        if len(responders) == 0 or len(nonresponders) == 0:
            raise ValueError(f"Both response groups are required for {population}")
        test = mannwhitneyu(responders, nonresponders, alternative="two-sided", method="asymptotic")
        if not math.isfinite(test.pvalue) or not math.isfinite(test.statistic):
            raise ValueError(f"Invalid Mann–Whitney result for {population}")
        responder_median = float(pd.Series(responders).median())
        nonresponder_median = float(pd.Series(nonresponders).median())
        results.append({
            "population": population,
            "n_responder_observations": len(responders),
            "n_nonresponder_observations": len(nonresponders),
            "n_responder_subjects": subset.loc[subset.response == "yes", "subject_id"].nunique(),
            "n_nonresponder_subjects": subset.loc[subset.response == "no", "subject_id"].nunique(),
            "responder_mean": float(responders.mean()),
            "nonresponder_mean": float(nonresponders.mean()),
            "responder_median": responder_median,
            "nonresponder_median": nonresponder_median,
            "responder_std": float(responders.std(ddof=1)),
            "nonresponder_std": float(nonresponders.std(ddof=1)),
            "responder_iqr": float(pd.Series(responders).quantile(.75) - pd.Series(responders).quantile(.25)),
            "nonresponder_iqr": float(pd.Series(nonresponders).quantile(.75) - pd.Series(nonresponders).quantile(.25)),
            "effect_size": float(2 * test.statistic / (len(responders) * len(nonresponders)) - 1),
            "raw_p_value": float(test.pvalue),
            "direction": "responder median higher" if responder_median > nonresponder_median else
                         "responder median lower" if responder_median < nonresponder_median else "equal medians",
        })
    table = pd.DataFrame(results)
    table["fdr_q_value"] = fdr_bh(table.raw_p_value.tolist())
    table["nominal_significant"] = table.raw_p_value < .05
    table["fdr_significant"] = table.fdr_q_value < .05
    return table


def subject_mean_data(cohort: pd.DataFrame) -> pd.DataFrame:
    """Give each subject one mean percentage across available timepoints."""
    return cohort.groupby(["subject_id", "response", "population"], as_index=False)["percentage"].mean()


def baseline_data(cohort: pd.DataFrame) -> pd.DataFrame:
    """Select day-zero samples, with one observation per subject/population."""
    baseline = cohort.loc[cohort.time_from_treatment_start == 0].copy()
    if baseline.duplicated(["subject_id", "population"]).any():
        raise ValueError("Multiple baseline samples for a subject and population")
    if baseline.subject_id.nunique() != cohort.subject_id.nunique():
        raise ValueError("Some subjects have no baseline sample")
    return baseline


def subject_change_data(cohort: pd.DataFrame) -> pd.DataFrame:
    """Compare each subject's day-14 minus day-0 percentage-point change."""
    wide = cohort.pivot(index=["subject_id", "response", "population"],
                        columns="time_from_treatment_start", values="percentage")
    if 0 not in wide or 14 not in wide or wide[[0, 14]].isna().any().any():
        raise ValueError("Day-zero and day-14 values are required for every subject")
    return (wide[14] - wide[0]).rename("percentage").reset_index()


def make_boxplots(cohort: pd.DataFrame, results: pd.DataFrame, path: Path) -> None:
    """Plot pooled sample distributions with their calculated statistical results."""
    import numpy as np

    labels = {"b_cell": "B cells", "cd8_t_cell": "CD8+ T cells",
              "cd4_t_cell": "CD4+ T cells", "nk_cell": "NK cells",
              "monocyte": "Monocytes"}
    colors = ("#a7b0b8", "#477a9b")
    statistics = results.set_index("population")
    random = np.random.default_rng(2026)
    upper_limit = max(50, float(cohort.percentage.max()) * 1.15)

    fig = plt.figure(figsize=(13.5, 7.4), facecolor="white")
    grid = fig.add_gridspec(2, 6, left=.07, right=.99, top=.81, bottom=.13,
                            hspace=.6, wspace=.48)
    positions = ((0, slice(0, 2)), (0, slice(2, 4)), (0, slice(4, 6)),
                 (1, slice(1, 3)), (1, slice(3, 5)))
    first_axis = None
    for population, (row, columns) in zip(POPULATION_ORDER, positions):
        axis = fig.add_subplot(grid[row, columns], sharey=first_axis)
        if first_axis is None:
            first_axis = axis
        subset = cohort.loc[cohort.population == population]
        groups = [subset.loc[subset.response == response, "percentage"].to_numpy()
                  for response in ("no", "yes")]
        for position, (values, color) in enumerate(zip(groups, colors), start=1):
            jitter = random.uniform(-.19, .19, len(values))
            axis.scatter(position + jitter, values, s=5, color=color,
                         alpha=.16, edgecolors="none", zorder=1)
        boxes = axis.boxplot(groups, widths=.38, patch_artist=True, showfliers=False,
                             medianprops={"color": "#202934", "linewidth": 1.5},
                             boxprops={"color": "#34404b", "linewidth": 1},
                             whiskerprops={"color": "#34404b"},
                             capprops={"color": "#34404b"})
        for patch, color in zip(boxes["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(.85)
        result = statistics.loc[population]
        axis.text(.5, .98,
                  f"p = {result.raw_p_value:.3g}  |  FDR q = {result.fdr_q_value:.3g}",
                  transform=axis.transAxes, ha="center", va="top", fontsize=8.5,
                  color="#34404b")
        axis.set_title(labels[population], fontsize=11, pad=8)
        axis.set_xticks((1, 2),
                        (f"Nonresponders\n(n = {len(groups[0])} samples)",
                         f"Responders\n(n = {len(groups[1])} samples)"))
        axis.set_ylim(0, upper_limit)
        axis.set_axisbelow(True)
        axis.grid(axis="y", color="#e5e9ed", linewidth=.6)
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(axis="both", labelsize=8.5)

    fig.suptitle("Immune cell frequencies by response to miraclib", y=.965, fontsize=13)
    fig.text(.5, .915, "Melanoma PBMC samples", ha="center", fontsize=10, color="#52606d")
    fig.supylabel("Relative frequency (%)", x=.025, fontsize=10)
    fig.text(.5, .035, "Pooled samples; repeated measures may occur within subjects.",
             ha="center", fontsize=8.5, color="#52606d")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def conclusions(table: pd.DataFrame, label: str) -> list[str]:
    """Describe calculated associations without treating them as causal."""
    lines = []
    for result in table.itertuples():
        evidence = ("survived FDR correction" if result.fdr_significant else
                    "nominal only; did not survive FDR" if result.nominal_significant else
                    "no nominal evidence")
        lines.append(f"{label}: {result.population} — {evidence}; {result.direction}; "
                     f"rank-biserial={result.effect_size:.3f}, "
                     f"p={result.raw_p_value:.4g}, q={result.fdr_q_value:.4g}.")
    return lines


def run_analysis(connection: sqlite3.Connection) -> dict:
    """Return dashboard-ready cohort, summaries, and all three comparisons."""
    cohort = load_cohort(connection)
    summary = summarize_cohort(cohort)
    sample_names = {"n_responder_observations": "n_responder_samples",
                    "n_nonresponder_observations": "n_nonresponder_samples"}
    pooled = compare_groups(cohort).rename(columns=sample_names)
    subject_mean = compare_groups(subject_mean_data(cohort))
    baseline = compare_groups(baseline_data(cohort)).rename(columns=sample_names)
    subject_change = compare_groups(subject_change_data(cohort))
    return {"cohort": cohort, "summary": summary, "pooled": pooled,
            "subject_mean": subject_mean, "baseline": baseline,
            "subject_change": subject_change}


def main() -> None:
    if not DB_PATH.is_file():
        raise FileNotFoundError(f"Database not found: {DB_PATH}. Run load_data.py first.")
    with closing(sqlite3.connect(DB_PATH)) as connection:
        analysis = run_analysis(connection)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    analysis["pooled"].to_csv(OUTPUT_DIR / "response_analysis_results.csv", index=False)
    analysis["subject_mean"].to_csv(OUTPUT_DIR / "subject_mean_sensitivity_results.csv", index=False)
    analysis["baseline"].to_csv(OUTPUT_DIR / "baseline_sensitivity_results.csv", index=False)
    analysis["subject_change"].to_csv(OUTPUT_DIR / "subject_change_sensitivity_results.csv", index=False)
    make_boxplots(analysis["cohort"], analysis["pooled"],
                  OUTPUT_DIR / "figures" / "response_frequency_boxplots.png")

    summary = analysis["summary"]
    print(f"Cohort: {summary['subjects']} subjects ({summary['responder_subjects']} responders, "
          f"{summary['nonresponder_subjects']} nonresponders); {summary['samples']} samples "
          f"({summary['responder_samples']} responders, {summary['nonresponder_samples']} nonresponders).")
    print(f"Samples per subject: {summary['samples_per_subject']}; timepoints: {summary['timepoints']}")
    print(f"Population observations: {summary['population_observations']}")
    print("\nRequested pooled sample-level comparison (exploratory; repeated samples are correlated):")
    print(analysis["pooled"][['population', 'responder_median', 'nonresponder_median',
                               'effect_size', 'raw_p_value', 'fdr_q_value']].to_string(index=False))
    for label, key in (("Pooled", "pooled"), ("Subject mean", "subject_mean"),
                       ("Baseline", "baseline"), ("Day 14 − day 0 change", "subject_change")):
        print(f"\n{label} conclusions:")
        for line in conclusions(analysis[key], label):
            print("  " + line)
    print("\nPercentages describe relative composition, so changes in one population affect the others. "
          "Subject-level mean and change comparisons give each patient one observation. "
          "Direction refers to medians; rank-biserial sign refers to pairwise dominance and can differ. "
          "Day zero is treated as baseline, but the CSV does not confirm pretreatment collection. "
          "Baseline associations alone do not establish predictive utility.")


if __name__ == "__main__":
    main()
