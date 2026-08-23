# Immune cell populations in the miraclib trial

A small analytics pipeline over `cell-count.csv`: it models the data in SQLite, computes
the relative frequency of five immune cell populations in every sample, tests whether
those frequencies differ between miraclib responders and non-responders, and serves the
whole thing through an interactive dashboard.

## Quick start

```bash
make setup       # install dependencies from requirements.txt
make pipeline    # build the database, then write every table and figure to outputs/
make dashboard   # serve the dashboard at http://localhost:8501
```

`make pipeline` runs start to finish with no manual intervention and takes a few seconds.
It is idempotent: the schema is dropped and rebuilt each run, so it can be run repeatedly.

Part 1 on its own is just:

```bash
python load_data.py     # creates cell_counts.db in the repository root
```

Requires Python 3.9+. `cell_counts.db` is deliberately not committed — it is a build
artifact that `make pipeline` regenerates from the CSV in seconds.

**Dashboard:** <!-- DASHBOARD_LINK -->_deployment link to be added_<!-- /DASHBOARD_LINK -->
(runs locally with `make dashboard` in the meantime).

## What the pipeline produces

| File in `outputs/` | Part | Contents |
| --- | --- | --- |
| `summary_table.csv` | 2 | One row per sample-population pair: `sample, total_count, population, count, percentage` (52,500 rows) |
| `part3_responder_stats.csv` | 3 | Responder vs non-responder test per population, all timepoints pooled |
| `part3_timecourse_stats.csv` | 3 | The same comparison run six more ways: per subject, per timepoint, on-treatment, and change from baseline |
| `part3_boxplots.png` | 3 | Boxplot per population, responders beside non-responders |
| `part3_boxplots_by_timepoint.png` | 3 | The same populations split by day 0 / 7 / 14 |
| `part4_baseline_samples.csv` | 4 | The melanoma / PBMC / miraclib / baseline subset |
| `part4_breakdowns.csv` | 4 | That subset counted by project, by response, and by sex |

## Database schema

```
projects ──< subjects ──< samples ──< cell_counts >── populations
```

| Table | Grain | Columns |
| --- | --- | --- |
| `projects` | one project | `project_id` |
| `subjects` | one patient | `subject_id`, `project_id`, `condition`, `age`, `sex`, `treatment`, `response` |
| `samples` | one specimen at one timepoint | `sample_id`, `subject_id`, `sample_type`, `time_from_treatment_start` |
| `populations` | one measurable cell type | `population`, `display_name` |
| `cell_counts` | one measurement | `sample_id`, `population`, `count` |

Two views sit on top: `sample_details` flattens sample → subject → project so that a
cohort filter does not have to repeat the same joins, and `sample_frequencies` is the
Part 2 calculation itself.

### Why this shape

**Clinical attributes belong to the patient, not the specimen.** The CSV repeats
`condition`, `age`, `sex`, `treatment` and `response` on all three of a subject's rows.
Storing them once on `subjects` removes the update anomaly where one sample could come to
disagree with another about the same patient's response. The loader verifies that the
repetition really was consistent before collapsing it, so this normalization can never
quietly discard a disagreement.

**The cytometry panel is data, not schema.** This is the one design decision worth
arguing about. Keeping five count columns on the sample row would mirror the CSV and make
loading marginally simpler. Instead `cell_counts` is long — one row per
(sample, population) — because:

- Adding a population becomes an `INSERT`, not an `ALTER TABLE` plus a migration of every
  downstream query. Real panels carry dozens to hundreds of populations and grow over time.
- Analytics stay uniform. The Part 2 percentage is one window function over one table
  rather than an unpivot in pandas, and it is expressible in SQL, which is what lets it
  live in a view that every consumer shares.
- Sparsity is free. A panel that measures a population on some samples but not others
  costs nothing, where a wide table would fill with NULLs.

The cost is that a human reading `cell_counts` directly sees five rows where they expected
one. `sample_frequencies` gives back the readable per-sample view.

**`response` is NULL, not `''`.** Healthy subjects are untreated and have no response
assessment. Storing NULL means `WHERE response = 'no'` cannot accidentally sweep them in,
and `COUNT(DISTINCT ...)` grouped by response reports them separately rather than as a
third category that looks like data.

### How this scales

At hundreds of projects and thousands of samples this schema is still comfortable — that
is roughly 10⁵–10⁶ rows in `cell_counts`, which SQLite handles without complaint given the
indexes in `schema.sql` (foreign keys, plus a composite index on the
`condition/treatment/response` cohort filter and on `sample_type/time_from_treatment_start`).

What changes as it grows:

- **Analytics shape.** `cell_counts` is already a fact table and `subjects`/`samples`/
  `populations` are already dimensions, so the model is a star schema. That is the shape
  warehouses want, so moving to DuckDB, BigQuery or Redshift is a load, not a redesign.
- **Concurrency, not volume, is what forces the move off SQLite.** One writer at a time is
  the real ceiling. A pipeline that ingests several projects concurrently, or a dashboard
  serving many analysts, wants Postgres; nothing in the schema or the queries has to change,
  because everything here is standard SQL.
- **Repeated aggregations get materialized.** `sample_frequencies` recomputes on every
  read, which is correct and cheap now. At scale it becomes a materialized table refreshed
  by the pipeline, or a column on `cell_counts` maintained at load time. Because consumers
  already go through the view rather than reimplementing the arithmetic, that swap is
  invisible to them.
- **New measurement types slot in.** A second panel, or a different assay entirely, adds
  rows to `populations` and `cell_counts` rather than columns to `samples`. If assays
  diverge enough to need their own metadata, `cell_counts` gains an `assay_id` dimension.
- **Provenance is the obvious next dimension.** Batch, instrument and operator are what
  you reach for first when a result looks odd, and they attach to `samples` as foreign keys
  to their own small tables.

## Code structure

```
load_data.py        Part 1: build the database from the CSV (root, no arguments)
dashboard.py        Streamlit app
src/
  schema.sql        the DDL, kept as SQL rather than embedded in Python
  config.py         paths, the panel, the significance threshold
  db.py             connection and query helpers
  summary.py        Part 2
  stats.py          Part 3
  subsets.py        Part 4
outputs/            generated tables and figures
```

The organizing idea is that **the database is the interface between stages**. Each stage
reads from SQLite and writes to `outputs/`; none of them passes DataFrames to another or
depends on having been run in the same process. That is what makes `make pipeline` a plain
sequence of commands, and it is why the dashboard can be started independently — it reads
the same database rather than recomputing anything.

A few consequences worth naming:

- **Definitions live in one place.** The relative-frequency calculation is in
  `schema.sql`, not in Python, so the exported CSV, the statistics and the dashboard cannot
  drift apart. Anything more than one stage must agree on — the panel, the paths, the alpha
  level — is in `config.py`.
- **The dashboard calls the analysis code, it does not reimplement it.** The Part 3 tab
  imports `src/stats.py` and calls the same `compare_groups` the pipeline calls, which is
  why choosing a cohort in the UI gives exactly the answer the pipeline would print for
  that cohort. The cohort query is parameterized for this reason.
- **The loader verifies itself.** A partial load is the expensive failure, because every
  downstream number would still look plausible. `load_data.py` checks row counts against
  the source, confirms the subject-level attributes really were constant before collapsing
  them, and confirms the total cell count survived the wide-to-long reshape.
- **SQL stays SQL.** Filters are expressed as SQL in the module that uses them rather than
  built out of dictionaries, because the queries are the substance of Parts 3 and 4 and are
  easier to check when you can read them whole.

## Results

### Part 2

Every sample's five counts are summed and each population expressed as a percentage of
that total, giving 52,500 rows over 10,500 samples. Percentages are rounded to two decimals
on export only; the statistics read unrounded values from the view.

### Part 3

Cohort: melanoma patients on miraclib, PBMC samples only — 1,968 samples from 656 subjects.
Comparisons use a two-sided Mann-Whitney U test (relative frequencies are bounded
percentages with no reason to be normal), with Benjamini-Hochberg correction across the
five populations. Reported effect size is the rank-biserial correlation.

**The pooled comparison finds nothing significant** after correction — the smallest
adjusted p-value is 0.067, for CD4+ T cells. That result is real but it is the wrong
average for the question, and the timecourse says why:

| | Baseline (day 0) | On treatment (days 7 and 14) | Change from own baseline |
| --- | --- | --- | --- |
| CD4+ T cell | 29.63% vs 29.53%, q = 0.89 | **30.88% vs 29.80%, q = 0.008** | +1.16 vs +0.08, q = 0.27 |
| B cell | 9.79% vs 9.76%, q = 0.89 | 9.59% vs 9.88%, q = 0.064 | **−1.00 vs +0.15, q = 0.031** |

(responder value first; bold marks significance at FDR < 0.05)

Three findings follow, and they are consistent with each other:

1. **Nothing at baseline separates the groups.** Before miraclib is given, no population
   distinguishes future responders from non-responders — the smallest adjusted p-value
   across all five is 0.89. On this cohort, a pre-treatment frequency on its own does not
   predict who will respond.
2. **Once treatment is underway, the groups separate.** Responders carry a higher CD4+ T
   cell fraction (30.88% vs 29.80%, q = 0.008).
3. **The B cell effect is a within-patient shift.** Responders' B cell fraction falls a
   full percentage point from their own baseline while non-responders drift slightly up
   (−1.00 vs +0.15, q = 0.031). Differencing against each patient's own day 0 removes
   between-patient variation, which is why the effect is clearest measured that way.

Pooling all three timepoints averages the pre-treatment draw, where there is no
difference, into the on-treatment draws, where there is — which dilutes both effects
below significance. The pooled result is not evidence of no effect; it is the wrong
average. The dashboard lets you switch between these views directly.

Two caveats stated plainly: each subject contributes three samples, so the pooled test
treats one patient as three observations, which is why every finding above is reported at
subject level. And these are associations within one trial arm — the design cannot
separate "responders' immune profile shifts" from "the shift is what response consists of".

### Part 4

Melanoma PBMC samples at baseline from miraclib-treated patients: **656 samples** from 656
subjects (each subject contributes exactly one baseline sample).

| Breakdown | Unit | Counts |
| --- | --- | --- |
| Samples per project | samples | prj1: 384, prj3: 272 (prj2 contributes none) |
| Responders vs non-responders | subjects | yes: 331, no: 325 |
| Males vs females | subjects | M: 344, F: 312 |

Samples are counted per project while response and sex are counted per subject, following
the wording of the question. Those two units coincide here, but they are different
questions and the queries keep them distinct.

A related question asked alongside the assignment — melanoma males of all sample and
treatment types, mean B cell count for responders at time 0 — is computed by the same
module: **10206.15**, over 485 samples.

## Dashboard

Four tabs, reading directly from the database:

- **Overview** — what is in the database, and how samples break down by condition,
  treatment, project and response.
- **Part 2** — the summary table with filters on condition, treatment, specimen type and
  timepoint, plus a distribution plot and a CSV download.
- **Part 3** — cohort selectors and the five analysis views described above. Boxplots and
  the statistics table recompute live, so the timecourse argument is something you can
  check rather than take on faith.
- **Part 4** — the baseline subset with its three breakdowns and the underlying rows.
