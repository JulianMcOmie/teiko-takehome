"""Part 3 - do cell population frequencies differ between miraclib responders and non-responders?

Cohort: melanoma patients treated with miraclib, PBMC samples only.

The headline analysis is the one Part 3 asks for directly: every sample in that cohort,
responders against non-responders, one test per population. Three design choices sit
underneath it, and each is reported rather than assumed.

* Multiple testing. Five populations means five tests, so p-values are adjusted with
  Benjamini-Hochberg and only the adjusted values are used to call significance.
* Non-independence. Each subject contributes three samples (days 0, 7 and 14), so the
  pooled comparison treats one patient as three observations. The same comparison run
  on subject-level averages is reported next to it.
* Timepoint. Pooling mixes the pre-treatment draw into a comparison about a drug's
  effect. Splitting by timepoint turns out to matter a great deal here, so the
  timecourse analyses below are reported as part of the answer rather than as an
  appendix: responders and non-responders are indistinguishable at baseline and
  separate only once miraclib is underway.

The test throughout is Mann-Whitney U. Relative frequencies are bounded percentages
with no particular reason to be normally distributed, and a rank test avoids assuming
they are.
"""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from src.config import ALPHA, OUTPUT_DIR, POPULATIONS
from src.db import query

COHORT_SQL = """
SELECT
    f.sample,
    f.population,
    f.percentage,
    d.subject_id,
    d.response,
    d.time_from_treatment_start
FROM sample_frequencies f
JOIN sample_details d ON d.sample_id = f.sample
WHERE d.condition   = ?
  AND d.treatment   = ?
  AND d.sample_type = ?
  AND d.response IS NOT NULL
ORDER BY f.sample, f.population
"""

# The cohort Part 3 specifies. Parameterising the query lets the dashboard point the
# same analysis at another cohort without a second copy of the SQL.
PART3_COHORT = ("melanoma", "miraclib", "PBMC")

BASELINE, ON_TREATMENT = 0, (7, 14)

STATS_PATH = OUTPUT_DIR / "part3_responder_stats.csv"
TIMECOURSE_PATH = OUTPUT_DIR / "part3_timecourse_stats.csv"
BOXPLOT_PATH = OUTPUT_DIR / "part3_boxplots.png"
TIMECOURSE_PLOT_PATH = OUTPUT_DIR / "part3_boxplots_by_timepoint.png"


def load_cohort(condition=None, treatment=None, sample_type=None):
    """Relative frequencies for every sample in a cohort that has a response recorded.

    Defaults to the Part 3 cohort: melanoma patients on miraclib, PBMC samples only.
    """
    default_condition, default_treatment, default_sample_type = PART3_COHORT
    return query(
        COHORT_SQL,
        (
            condition or default_condition,
            treatment or default_treatment,
            sample_type or default_sample_type,
        ),
    )


def by_subject(cohort):
    """One value per subject per population: the mean across their three timepoints."""
    return cohort.groupby(
        ["subject_id", "population", "response"], as_index=False
    )["percentage"].mean()


def subject_timecourse(cohort):
    """One row per subject and population, with each timepoint in its own column.

    Also derives the two summaries the timecourse analyses need: the mean of the
    on-treatment draws, and the change from the patient's own baseline. Differencing
    against a patient's own day 0 removes between-patient variation, which makes it
    the most sensitive view of what the drug actually did.
    """
    wide = cohort.pivot_table(
        index=["subject_id", "response", "population"],
        columns="time_from_treatment_start",
        values="percentage",
    ).reset_index()
    wide["on_treatment"] = wide[list(ON_TREATMENT)].mean(axis=1)
    wide["change_from_baseline"] = wide[ON_TREATMENT[-1]] - wide[BASELINE]
    return wide


def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg adjusted p-values, controlling the false discovery rate."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    order = np.argsort(p)
    scaled = p[order] * n / np.arange(1, n + 1)
    # Adjusted values must be non-decreasing in p, so sweep back from the largest.
    scaled = np.minimum.accumulate(scaled[::-1])[::-1]
    adjusted = np.empty(n)
    adjusted[order] = np.clip(scaled, 0, 1)
    return adjusted


def compare_groups(frame, value="percentage"):
    """Compare responders against non-responders for each population.

    Returns one row per population with group sizes, medians, the Mann-Whitney U
    statistic, raw and FDR-adjusted p-values, and a rank-biserial effect size. The
    effect size is positive when responders carry the higher value.
    """
    rows = []
    for population in POPULATIONS:
        values = frame[frame["population"] == population]
        responder = values.loc[values["response"] == "yes", value]
        non_responder = values.loc[values["response"] == "no", value]

        statistic, p_value = mannwhitneyu(
            responder, non_responder, alternative="two-sided"
        )

        rows.append(
            {
                "population": population,
                "n_responder": len(responder),
                "n_non_responder": len(non_responder),
                "median_responder": responder.median(),
                "median_non_responder": non_responder.median(),
                "median_difference": responder.median() - non_responder.median(),
                "u_statistic": statistic,
                "p_value": p_value,
                "rank_biserial": 2 * statistic / (len(responder) * len(non_responder)) - 1,
            }
        )

    results = pd.DataFrame(rows)
    results["p_adjusted"] = benjamini_hochberg(results["p_value"])
    results["significant"] = results["p_adjusted"] < ALPHA
    return results.sort_values("p_value").reset_index(drop=True)


def timecourse_analyses(cohort, wide):
    """Run the supplementary comparisons and stack them into one tidy table.

    Each analysis is corrected for multiple testing within itself, because each answers
    a separate question about five populations.
    """
    analyses = {
        "subject_average": (by_subject(cohort), "percentage"),
        "baseline_only": (
            cohort[cohort["time_from_treatment_start"] == BASELINE],
            "percentage",
        ),
        "on_treatment": (wide, "on_treatment"),
        "change_from_baseline": (wide, "change_from_baseline"),
    }
    for day in sorted(cohort["time_from_treatment_start"].unique()):
        analyses[f"day_{day}"] = (
            cohort[cohort["time_from_treatment_start"] == day],
            "percentage",
        )

    stacked = []
    for name, (frame, value) in analyses.items():
        result = compare_groups(frame, value)
        result.insert(0, "analysis", name)
        stacked.append(result)
    return pd.concat(stacked, ignore_index=True)


def _style_boxes(boxes, colours):
    for patch, colour in zip(boxes["boxes"], colours):
        patch.set_facecolor(colour)
        patch.set_alpha(0.65)


RESPONDER_COLOUR, NON_RESPONDER_COLOUR = "#2A9D8F", "#B25A7E"


def _figure():
    import matplotlib

    matplotlib.use("Agg")  # Render without a display; the pipeline runs headless.
    import matplotlib.pyplot as plt

    return plt


def make_boxplot(cohort, results, path=BOXPLOT_PATH):
    """Part 3 as asked: one boxplot per population, responders beside non-responders."""
    plt = _figure()
    lookup = results.set_index("population")
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(16, 4.8))

    for axis, (population, label) in zip(axes, POPULATIONS.items()):
        values = cohort[cohort["population"] == population]
        boxes = axis.boxplot(
            [
                values.loc[values["response"] == "yes", "percentage"],
                values.loc[values["response"] == "no", "percentage"],
            ],
            labels=["Responder", "Non-responder"],
            patch_artist=True,
            widths=0.55,
            medianprops={"color": "#1D3B36", "linewidth": 1.6},
            flierprops={"marker": ".", "markersize": 3, "alpha": 0.35},
        )
        _style_boxes(boxes, [RESPONDER_COLOUR, NON_RESPONDER_COLOUR])

        row = lookup.loc[population]
        axis.set_title(
            f"{label}{' *' if row['significant'] else ''}\nq = {row['p_adjusted']:.3g}",
            fontsize=10,
        )
        axis.set_ylabel("Relative frequency (%)" if population == "b_cell" else "")
        axis.tick_params(axis="x", labelsize=9)
        axis.grid(axis="y", alpha=0.25, linewidth=0.6)
        axis.set_axisbelow(True)

    fig.suptitle(
        "Melanoma patients on miraclib (PBMC), all timepoints pooled", fontsize=13
    )
    fig.text(
        0.5,
        0.015,
        f"Mann-Whitney U, Benjamini-Hochberg adjusted. * marks q < {ALPHA}.",
        ha="center",
        fontsize=9,
        color="#4A5C60",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.94])
    path.parent.mkdir(exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def make_timecourse_plot(cohort, timecourse, path=TIMECOURSE_PLOT_PATH):
    """The same populations split by timepoint, which is where the story actually is."""
    plt = _figure()
    days = sorted(cohort["time_from_treatment_start"].unique())
    lookup = timecourse.set_index(["analysis", "population"])
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(17, 4.8))

    for axis, (population, label) in zip(axes, POPULATIONS.items()):
        for offset, response, colour in [
            (-0.19, "yes", RESPONDER_COLOUR),
            (0.19, "no", NON_RESPONDER_COLOUR),
        ]:
            groups = [
                cohort.loc[
                    (cohort["population"] == population)
                    & (cohort["time_from_treatment_start"] == day)
                    & (cohort["response"] == response),
                    "percentage",
                ]
                for day in days
            ]
            boxes = axis.boxplot(
                groups,
                positions=[i + offset for i in range(len(days))],
                widths=0.32,
                patch_artist=True,
                medianprops={"color": "#1D3B36", "linewidth": 1.4},
                flierprops={"marker": ".", "markersize": 2.5, "alpha": 0.3},
            )
            _style_boxes(boxes, [colour] * len(days))

        marks = " ".join(
            "*" if lookup.loc[(f"day_{day}", population), "significant"] else ""
            for day in days
        )
        axis.set_xticks(range(len(days)))
        axis.set_xticklabels([f"Day {day}" for day in days], fontsize=9)
        axis.set_title(f"{label} {marks}".strip(), fontsize=10)
        axis.set_ylabel("Relative frequency (%)" if population == "b_cell" else "")
        axis.grid(axis="y", alpha=0.25, linewidth=0.6)
        axis.set_axisbelow(True)

    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=RESPONDER_COLOUR, alpha=0.65),
        plt.Rectangle((0, 0), 1, 1, facecolor=NON_RESPONDER_COLOUR, alpha=0.65),
    ]
    fig.legend(
        handles,
        ["Responder", "Non-responder"],
        loc="upper right",
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    fig.suptitle(
        "Responders and non-responders separate only after miraclib begins",
        fontsize=13,
    )
    fig.text(
        0.5,
        0.015,
        "* marks q < %.2f within that timepoint (Mann-Whitney U, Benjamini-Hochberg)."
        % ALPHA,
        ha="center",
        fontsize=9,
        color="#4A5C60",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.92])
    path.parent.mkdir(exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


DISPLAY_COLUMNS = [
    "population",
    "n_responder",
    "n_non_responder",
    "median_responder",
    "median_non_responder",
    "median_difference",
    "p_value",
    "p_adjusted",
    "rank_biserial",
    "significant",
]

ROUNDING = {
    "median_responder": 2,
    "median_non_responder": 2,
    "median_difference": 2,
    "p_value": 4,
    "p_adjusted": 4,
    "rank_biserial": 3,
}


def describe(results, label):
    print(f"\n{label}")
    print(results[DISPLAY_COLUMNS].round(ROUNDING).to_string(index=False))


def report_conclusions(pooled, timecourse):
    """Print the findings in the order Bob would want to hear them."""
    at = lambda analysis: timecourse[timecourse["analysis"] == analysis].set_index(
        "population"
    )
    baseline, on_treatment, change = (
        at("baseline_only"),
        at("on_treatment"),
        at("change_from_baseline"),
    )

    print("\n" + "=" * 78)
    print("FINDINGS")
    print("=" * 78)

    if not baseline["significant"].any():
        print(
            "\n1. Nothing at baseline separates the two groups. Before miraclib is given,\n"
            "   no population differs between future responders and non-responders\n"
            f"   (smallest q = {baseline['p_adjusted'].min():.2f}). On this cohort, a\n"
            "   pre-treatment frequency on its own does not predict who will respond."
        )

    movers = on_treatment[on_treatment["significant"]].sort_values("p_adjusted")
    if not movers.empty:
        print("\n2. Once on treatment, the groups do separate:")
        for population, row in movers.iterrows():
            direction = "higher" if row["median_difference"] > 0 else "lower"
            print(
                f"   {POPULATIONS[population]}: {direction} in responders "
                f"({row['median_responder']:.2f}% vs {row['median_non_responder']:.2f}%), "
                f"q = {row['p_adjusted']:.4f}"
            )

    shifts = change[change["significant"]].sort_values("p_adjusted")
    if not shifts.empty:
        print("\n3. Measured as change from each patient's own baseline:")
        for population, row in shifts.iterrows():
            print(
                f"   {POPULATIONS[population]}: responders {row['median_responder']:+.2f} "
                f"points vs {row['median_non_responder']:+.2f} in non-responders, "
                f"q = {row['p_adjusted']:.4f}"
            )

    pooled_significant = pooled[pooled["significant"]]
    print(
        "\n4. The pooled comparison across all timepoints - the one Part 3 asks for\n"
        "   directly - finds "
        + (
            "no population significant after correction"
            if pooled_significant.empty
            else ", ".join(pooled_significant["population"])
        )
        + f" (smallest q = {pooled['p_adjusted'].min():.3f}).\n"
        "   That is expected rather than contradictory: pooling averages the\n"
        "   pre-treatment draw, where there is no difference, into the on-treatment\n"
        "   draws, where there is. Splitting by timepoint recovers the signal."
    )
    print("=" * 78)


def main():
    cohort = load_cohort()
    wide = subject_timecourse(cohort)

    pooled = compare_groups(cohort)
    timecourse = timecourse_analyses(cohort, wide)

    OUTPUT_DIR.mkdir(exist_ok=True)
    pooled.to_csv(STATS_PATH, index=False)
    timecourse.to_csv(TIMECOURSE_PATH, index=False)
    make_boxplot(cohort, pooled)
    make_timecourse_plot(cohort, timecourse)

    print(
        f"Part 3 - melanoma / miraclib / PBMC: {cohort['sample'].nunique():,} samples "
        f"from {cohort['subject_id'].nunique():,} subjects"
    )
    describe(pooled, "Pooled across timepoints, per sample (the comparison Part 3 asks for)")
    for analysis, label in [
        ("subject_average", "Per subject, timepoints averaged (non-independence check)"),
        ("baseline_only", "Baseline only (day 0)"),
        ("on_treatment", "On treatment (days 7 and 14 averaged, per subject)"),
        ("change_from_baseline", "Change from own baseline (day 14 minus day 0)"),
    ]:
        describe(timecourse[timecourse["analysis"] == analysis], label)

    report_conclusions(pooled, timecourse)
    print(
        f"\nWritten to {STATS_PATH.name}, {TIMECOURSE_PATH.name}, "
        f"{BOXPLOT_PATH.name}, {TIMECOURSE_PLOT_PATH.name}"
    )


if __name__ == "__main__":
    main()
