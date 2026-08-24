"""Interactive dashboard for the Teiko cell-count analysis.

Start it with `make dashboard` (or `streamlit run dashboard.py`). The app reads the
SQLite database the pipeline builds, so the database is the interface between the two:
nothing here recomputes a number that the pipeline defines differently.

The statistics tab reuses the functions in src/stats.py rather than reimplementing
them, which is why a cohort chosen in the sidebar produces the same answer the
pipeline would print for that cohort.
"""

import pandas as pd
import plotly.express as px
import streamlit as st

from src import stats as analysis
from src.config import ALPHA, CSV_PATH, DB_PATH, POPULATIONS
from src.db import query

RESPONDER_COLOUR, NON_RESPONDER_COLOUR = "#2A9D8F", "#B25A7E"
RESPONSE_LABELS = {"yes": "Responder", "no": "Non-responder"}

# One colour per population, in panel order, so a population keeps its colour
# wherever it appears in the app.
POPULATION_COLOURS = ["#2A9D8F", "#3D6E9C", "#7D4180", "#C4813C", "#5E8C4E"]

st.set_page_config(
    page_title="Teiko | Immune cell populations", page_icon="🧬", layout="wide"
)


@st.cache_resource(show_spinner=False)
def ensure_database():
    """Build the database if it is not there yet.

    The .db file is a build artifact and is not committed, so a fresh deployment has
    the CSV but no database and nothing has run `make pipeline`. Building it on first
    load takes a couple of seconds and makes the app self-sufficient, whether it is
    running on a hosted service or on a checkout where only `make dashboard` was run.
    """
    if DB_PATH.exists():
        return DB_PATH

    import load_data

    if not CSV_PATH.exists():
        return None
    load_data.main()
    return DB_PATH


@st.cache_data(show_spinner=False)
def run_query(sql, params=()):
    return query(sql, params)


@st.cache_data(show_spinner=False)
def cohort_frequencies(condition, treatment, sample_type):
    return analysis.load_cohort(condition, treatment, sample_type)


def metadata():
    return run_query("SELECT * FROM sample_details ORDER BY sample_id")


def distinct(column):
    rows = run_query(f"SELECT DISTINCT {column} AS v FROM sample_details ORDER BY v")
    return [v for v in rows["v"].tolist() if v is not None]


def population_label(series):
    return series.map(POPULATIONS)


def significance_note(results):
    """One sentence summarising which populations survived correction."""
    hits = results[results["significant"]]
    if hits.empty:
        return (
            f"No population differs significantly at FDR < {ALPHA} "
            f"(smallest adjusted p = {results['p_adjusted'].min():.3f})."
        )
    parts = [
        f"**{POPULATIONS[row.population]}** "
        f"({'higher' if row.median_difference > 0 else 'lower'} in responders, "
        f"q = {row.p_adjusted:.4f})"
        for row in hits.itertuples()
    ]
    return "Significant after Benjamini-Hochberg correction: " + ", ".join(parts) + "."


def style_results(results):
    """Render the statistics table with the columns Bob actually reads."""
    table = results.rename(
        columns={
            "population": "Population",
            "n_responder": "n responder",
            "n_non_responder": "n non-responder",
            "median_responder": "Median responder",
            "median_non_responder": "Median non-responder",
            "median_difference": "Difference",
            "p_value": "p",
            "p_adjusted": "q (FDR)",
            "rank_biserial": "Effect size",
            "significant": "Significant",
        }
    )
    table["Population"] = table["Population"].map(POPULATIONS)
    numeric = [
        "Median responder",
        "Median non-responder",
        "Difference",
        "Effect size",
    ]
    return table[
        ["Population", "n responder", "n non-responder", *numeric, "p", "q (FDR)", "Significant"]
    ].style.format(
        {
            **{column: "{:.2f}" for column in numeric},
            "p": "{:.4f}",
            "q (FDR)": "{:.4f}",
        }
    ).apply(
        lambda row: [
            "background-color: rgba(42,157,143,0.16)" if row["Significant"] else ""
        ]
        * len(row),
        axis=1,
    )


with st.spinner("Building the database from cell-count.csv…"):
    if ensure_database() is None:
        st.title("Immune cell populations")
        st.error(
            f"Neither `{DB_PATH.name}` nor `{CSV_PATH.name}` is present, so there is "
            "nothing to show. Check out the repository with its data file and run "
            "`make pipeline` from the root."
        )
        st.stop()

meta = metadata()

st.title("Immune cell populations in the miraclib trial")
st.caption(
    "Relative frequencies of five immune cell populations across "
    f"{meta['sample_id'].nunique():,} samples from {meta['subject_id'].nunique():,} "
    f"subjects in {meta['project_id'].nunique()} projects."
)

overview_tab, frequency_tab, response_tab, subset_tab = st.tabs(
    [
        "Overview",
        "Part 2 · Cell frequencies",
        "Part 3 · Response",
        "Part 4 · Baseline subset",
    ]
)


with overview_tab:
    st.subheader("What is in the database")

    columns = st.columns(4)
    columns[0].metric("Samples", f"{meta['sample_id'].nunique():,}")
    columns[1].metric("Subjects", f"{meta['subject_id'].nunique():,}")
    columns[2].metric("Projects", f"{meta['project_id'].nunique()}")
    columns[3].metric("Populations measured", len(POPULATIONS))

    left, right = st.columns(2)
    for container, column, title in [
        (left, "condition", "Samples by condition"),
        (right, "treatment", "Samples by treatment"),
    ]:
        counts = meta[column].value_counts().reset_index()
        counts.columns = [column, "samples"]
        figure = px.bar(
            counts, x=column, y="samples", title=title, color_discrete_sequence=[RESPONDER_COLOUR]
        )
        figure.update_layout(showlegend=False, xaxis_title="", height=320)
        container.plotly_chart(figure, use_container_width=True)

    left, right = st.columns(2)
    breakdown = (
        meta.groupby(["project_id", "sample_type"]).size().reset_index(name="samples")
    )
    figure = px.bar(
        breakdown,
        x="project_id",
        y="samples",
        color="sample_type",
        title="Samples by project and specimen type",
        color_discrete_sequence=[RESPONDER_COLOUR, NON_RESPONDER_COLOUR],
    )
    figure.update_layout(xaxis_title="", height=320)
    left.plotly_chart(figure, use_container_width=True)

    response_counts = (
        meta.assign(response=meta["response"].fillna("not assessed"))["response"]
        .value_counts()
        .reset_index()
    )
    response_counts.columns = ["response", "samples"]
    figure = px.bar(
        response_counts,
        x="response",
        y="samples",
        title="Samples by recorded response",
        color_discrete_sequence=[NON_RESPONDER_COLOUR],
    )
    figure.update_layout(showlegend=False, xaxis_title="", height=320)
    right.plotly_chart(figure, use_container_width=True)

    st.info(
        "Healthy subjects are untreated and have no response assessment, so their "
        "`response` is NULL rather than an empty string. SQL comparisons exclude them "
        "automatically instead of silently counting them as non-responders."
    )


with frequency_tab:
    st.subheader("Relative frequency of each population in each sample")
    st.caption(
        "For every sample the five counts are summed, and each population is shown as "
        "a percentage of that total. This is the Part 2 summary table."
    )

    filters = st.columns(4)
    condition = filters[0].multiselect(
        "Condition", distinct("condition"), default=distinct("condition")
    )
    treatment = filters[1].multiselect(
        "Treatment", distinct("treatment"), default=distinct("treatment")
    )
    sample_type = filters[2].multiselect(
        "Specimen type", distinct("sample_type"), default=distinct("sample_type")
    )
    timepoint = filters[3].multiselect(
        "Timepoint (days)",
        distinct("time_from_treatment_start"),
        default=distinct("time_from_treatment_start"),
    )

    summary = run_query(
        "SELECT sample, total_count, population, count, percentage FROM sample_frequencies"
    )
    selected = meta[
        meta["condition"].isin(condition)
        & meta["treatment"].isin(treatment)
        & meta["sample_type"].isin(sample_type)
        & meta["time_from_treatment_start"].isin(timepoint)
    ]
    table = summary[summary["sample"].isin(selected["sample_id"])]

    if table.empty:
        st.warning("No samples match the current filters.")
    else:
        st.caption(
            f"{table['sample'].nunique():,} samples · {len(table):,} rows "
            f"({table['sample'].nunique():,} samples x {len(POPULATIONS)} populations)"
        )

        figure = px.box(
            table.assign(Population=population_label(table["population"])),
            x="Population",
            y="percentage",
            color="Population",
            title="Distribution of relative frequency by population",
            category_orders={"Population": list(POPULATIONS.values())},
            color_discrete_sequence=POPULATION_COLOURS,
        )
        figure.update_layout(
            showlegend=False, yaxis_title="Relative frequency (%)", xaxis_title="", height=420
        )
        st.plotly_chart(figure, use_container_width=True)

        display = table.assign(percentage=table["percentage"].round(2))
        st.dataframe(display, use_container_width=True, hide_index=True, height=380)
        st.download_button(
            "Download this table as CSV",
            display.to_csv(index=False).encode(),
            file_name="summary_table.csv",
            mime="text/csv",
        )


with response_tab:
    st.subheader("Responders versus non-responders")

    selectors = st.columns(3)
    condition = selectors[0].selectbox(
        "Condition", distinct("condition"), index=distinct("condition").index("melanoma")
    )
    treatments = [t for t in distinct("treatment") if t != "none"]
    treatment = selectors[1].selectbox(
        "Treatment", treatments, index=treatments.index("miraclib")
    )
    sample_type = selectors[2].selectbox(
        "Specimen type", distinct("sample_type"), index=distinct("sample_type").index("PBMC")
    )

    default_cohort = (condition, treatment, sample_type) == analysis.PART3_COHORT
    if default_cohort:
        st.caption(
            "This is the cohort Part 3 specifies: melanoma patients on miraclib, "
            "PBMC samples only."
        )

    cohort = cohort_frequencies(condition, treatment, sample_type)

    if cohort.empty or cohort["response"].nunique() < 2:
        st.warning(
            "This cohort does not contain both responders and non-responders, so there "
            "is nothing to compare. Try melanoma / miraclib / PBMC."
        )
        st.stop()

    view = st.radio(
        "How to treat each patient's three timepoints",
        [
            "All timepoints pooled (per sample)",
            "Per subject (timepoints averaged)",
            "Baseline only (day 0)",
            "On treatment (days 7 and 14)",
            "Change from own baseline (day 14 − day 0)",
        ],
        horizontal=True,
        help=(
            "Each patient gives three samples, so pooling counts one patient three "
            "times. The other views address that, and separate the pre-treatment draw "
            "from the on-treatment ones."
        ),
    )

    wide = analysis.subject_timecourse(cohort)
    frames = {
        "All timepoints pooled (per sample)": (cohort, "percentage"),
        "Per subject (timepoints averaged)": (analysis.by_subject(cohort), "percentage"),
        "Baseline only (day 0)": (
            cohort[cohort["time_from_treatment_start"] == analysis.BASELINE],
            "percentage",
        ),
        "On treatment (days 7 and 14)": (wide, "on_treatment"),
        "Change from own baseline (day 14 − day 0)": (wide, "change_from_baseline"),
    }
    frame, value = frames[view]
    results = analysis.compare_groups(frame, value)

    plot_data = frame.assign(
        Population=population_label(frame["population"]),
        Response=frame["response"].map(RESPONSE_LABELS),
    )
    axis_title = (
        "Change in relative frequency (percentage points)"
        if value == "change_from_baseline"
        else "Relative frequency (%)"
    )
    figure = px.box(
        plot_data,
        x="Population",
        y=value,
        color="Response",
        color_discrete_map={
            "Responder": RESPONDER_COLOUR,
            "Non-responder": NON_RESPONDER_COLOUR,
        },
        category_orders={
            "Population": list(POPULATIONS.values()),
            "Response": ["Responder", "Non-responder"],
        },
    )
    figure.update_layout(
        boxmode="group", yaxis_title=axis_title, xaxis_title="", height=460
    )
    if value == "change_from_baseline":
        figure.add_hline(y=0, line_dash="dot", line_color="#7A8B88")
    st.plotly_chart(figure, use_container_width=True)

    st.markdown(significance_note(results))
    st.dataframe(style_results(results), use_container_width=True, hide_index=True)
    st.caption(
        "Mann-Whitney U, two-sided, with Benjamini-Hochberg correction across the five "
        "populations. Effect size is the rank-biserial correlation; positive means "
        "higher in responders."
    )

    if default_cohort:
        with st.expander("What the timecourse shows, and why the pooled test looks quiet"):
            st.markdown(
                """
Switch the view above between **Baseline only** and **On treatment** to see the
shape of the result.

At baseline nothing separates the two groups: before miraclib is given, no population
distinguishes future responders from non-responders. Once treatment is underway they
do separate — the CD4+ T cell fraction runs higher in responders, and the B cell
fraction falls away from each responder's own starting point while non-responders stay
flat.

Pooling all three timepoints, which is the comparison Part 3 asks for directly, mixes
the pre-treatment draw into the on-treatment ones and dilutes both effects. That is why
the pooled view shows nothing significant after correction while the timepoint-aware
views do. The pooled result is not evidence of no effect; it is the wrong average for
this question.
                """
            )


with subset_tab:
    st.subheader("Melanoma PBMC samples at baseline, from miraclib-treated patients")
    st.caption(
        "time_from_treatment_start = 0. Samples are counted per project; responders "
        "and sex are counted per subject, following the wording of Part 4."
    )

    from src import subsets

    baseline = run_query(subsets.BASELINE_SAMPLES_SQL)
    breakdowns = subsets.breakdowns()

    columns = st.columns(3)
    columns[0].metric("Baseline samples", f"{len(baseline):,}")
    columns[1].metric("Subjects", f"{baseline['subject_id'].nunique():,}")
    columns[2].metric("Projects represented", baseline["project_id"].nunique())

    charts = st.columns(3)
    specs = [
        ("samples_per_project", "Samples per project", "Project"),
        ("subjects_by_response", "Subjects by response", "Response"),
        ("subjects_by_sex", "Subjects by sex", "Sex"),
    ]
    for container, (key, title, axis_label) in zip(charts, specs):
        data = breakdowns[breakdowns["breakdown"] == key].copy()
        if key == "subjects_by_response":
            data["value"] = data["value"].map(RESPONSE_LABELS)
        figure = px.bar(
            data,
            x="value",
            y="n",
            title=title,
            text="n",
            color_discrete_sequence=[RESPONDER_COLOUR],
        )
        figure.update_layout(
            showlegend=False,
            xaxis_title=axis_label,
            yaxis_title="",
            height=330,
            # Headroom so the value labels above each bar are not clipped.
            yaxis_range=[0, data["n"].max() * 1.18],
        )
        figure.update_traces(textposition="outside")
        container.plotly_chart(figure, use_container_width=True)

    st.dataframe(baseline, use_container_width=True, hide_index=True, height=340)
    st.download_button(
        "Download the baseline subset as CSV",
        baseline.to_csv(index=False).encode(),
        file_name="part4_baseline_samples.csv",
        mime="text/csv",
    )

    focus = subsets.melanoma_male_baseline_bcell()
    st.metric(
        "Mean B cell count · melanoma males, responders, day 0, any specimen and treatment",
        f"{focus['mean_b_cell']:,.2f}",
        help=f"Averaged over {int(focus['n_samples']):,} samples.",
    )
