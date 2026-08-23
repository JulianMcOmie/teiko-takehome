"""Part 4 - baseline melanoma samples from miraclib-treated patients.

The base subset is every melanoma PBMC sample taken at baseline
(time_from_treatment_start = 0) from a patient treated with miraclib. The breakdowns
then follow the wording of the question: samples are counted per project, while
responders/non-responders and males/females are counted per *subject*. Those two units
coincide in this dataset, because each subject contributes exactly one baseline
sample, but they are different questions and the queries keep them distinct.
"""

import pandas as pd

from src.config import OUTPUT_DIR
from src.db import query

BASELINE_FILTER = """
    WHERE d.condition                 = 'melanoma'
      AND d.sample_type               = 'PBMC'
      AND d.treatment                 = 'miraclib'
      AND d.time_from_treatment_start = 0
"""

BASELINE_SAMPLES_SQL = f"""
SELECT
    d.sample_id,
    d.project_id,
    d.subject_id,
    d.condition,
    d.age,
    d.sex,
    d.treatment,
    d.response,
    d.sample_type,
    d.time_from_treatment_start
FROM sample_details d
{BASELINE_FILTER}
ORDER BY d.sample_id
"""

SAMPLES_PER_PROJECT_SQL = f"""
SELECT d.project_id, COUNT(*) AS n_samples
FROM sample_details d
{BASELINE_FILTER}
GROUP BY d.project_id
ORDER BY d.project_id
"""

SUBJECTS_BY_RESPONSE_SQL = f"""
SELECT d.response, COUNT(DISTINCT d.subject_id) AS n_subjects
FROM sample_details d
{BASELINE_FILTER}
GROUP BY d.response
ORDER BY d.response
"""

SUBJECTS_BY_SEX_SQL = f"""
SELECT d.sex, COUNT(DISTINCT d.subject_id) AS n_subjects
FROM sample_details d
{BASELINE_FILTER}
GROUP BY d.sex
ORDER BY d.sex
"""

# A separate question, deliberately looser than the subset above: melanoma males of
# every sample type and every treatment, responders only, at baseline.
MELANOMA_MALE_BASELINE_BCELL_SQL = """
SELECT
    COUNT(*)       AS n_samples,
    AVG(c.count)   AS mean_b_cell
FROM cell_counts c
JOIN sample_details d ON d.sample_id = c.sample_id
WHERE c.population                = 'b_cell'
  AND d.condition                 = 'melanoma'
  AND d.sex                       = 'M'
  AND d.response                  = 'yes'
  AND d.time_from_treatment_start = 0
"""

BASELINE_PATH = OUTPUT_DIR / "part4_baseline_samples.csv"
BREAKDOWN_PATH = OUTPUT_DIR / "part4_breakdowns.csv"


def baseline_samples():
    """The melanoma / PBMC / miraclib / baseline subset."""
    return query(BASELINE_SAMPLES_SQL)


def breakdowns():
    """The three requested counts, stacked into one tidy long table."""
    per_project = query(SAMPLES_PER_PROJECT_SQL).rename(
        columns={"project_id": "value", "n_samples": "n"}
    )
    per_project.insert(0, "breakdown", "samples_per_project")
    per_project.insert(2, "unit", "samples")

    by_response = query(SUBJECTS_BY_RESPONSE_SQL).rename(
        columns={"response": "value", "n_subjects": "n"}
    )
    by_response.insert(0, "breakdown", "subjects_by_response")
    by_response.insert(2, "unit", "subjects")

    by_sex = query(SUBJECTS_BY_SEX_SQL).rename(
        columns={"sex": "value", "n_subjects": "n"}
    )
    by_sex.insert(0, "breakdown", "subjects_by_sex")
    by_sex.insert(2, "unit", "subjects")

    return pd.concat([per_project, by_response, by_sex], ignore_index=True)


def melanoma_male_baseline_bcell():
    """Mean B cell count for male melanoma responders at baseline, any type/treatment."""
    return query(MELANOMA_MALE_BASELINE_BCELL_SQL).iloc[0]


def main():
    samples = baseline_samples()
    table = breakdowns()

    OUTPUT_DIR.mkdir(exist_ok=True)
    samples.to_csv(BASELINE_PATH, index=False)
    table.to_csv(BREAKDOWN_PATH, index=False)

    print(
        f"Part 4 - melanoma PBMC baseline samples on miraclib: {len(samples):,} samples "
        f"from {samples['subject_id'].nunique():,} subjects"
    )
    for name, group in table.groupby("breakdown", sort=False):
        unit = group["unit"].iloc[0]
        counts = ", ".join(f"{row.value}: {row.n:,}" for row in group.itertuples())
        print(f"  {name} ({unit}) -> {counts}")

    focus = melanoma_male_baseline_bcell()
    print(
        f"\nMelanoma males, responders, baseline (all sample and treatment types): "
        f"mean B cell count = {focus['mean_b_cell']:.2f} over {int(focus['n_samples']):,} samples"
    )
    print(f"\nWritten to {BASELINE_PATH.name}, {BREAKDOWN_PATH.name}")


if __name__ == "__main__":
    main()
