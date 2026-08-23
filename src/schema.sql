-- Schema for the Teiko cell-count database.
--
-- Shape: projects -> subjects -> samples -> cell_counts, with the measurements held
-- in a long (tall) fact table rather than five fixed count columns. The rationale is
-- in README.md; the short version is that a cytometry panel is data, not schema, so
-- adding a population must be an INSERT and never an ALTER TABLE.
--
-- Re-running this file drops and rebuilds everything, which keeps `make pipeline`
-- idempotent.

DROP VIEW  IF EXISTS sample_frequencies;
DROP VIEW  IF EXISTS sample_details;
DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS populations;
DROP TABLE IF EXISTS projects;

-- A project is one study or collection effort that contributes subjects.
CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

-- Clinical attributes describe the *person*, not an individual blood draw, and in
-- this dataset they are constant across a subject's three samples. Storing them here
-- rather than on every sample row removes the update anomaly where one sample could
-- disagree with another about the same patient's response.
CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects (project_id),
    condition  TEXT NOT NULL,
    age        INTEGER,
    sex        TEXT CHECK (sex IN ('M', 'F')),
    treatment  TEXT,
    -- NULL rather than '' for subjects with no response assessment (the healthy,
    -- untreated arm), so that SQL's three-valued logic excludes them automatically.
    response   TEXT CHECK (response IN ('yes', 'no'))
);

-- One row per biological sample: a single specimen from one subject at one timepoint.
CREATE TABLE samples (
    sample_id                 TEXT PRIMARY KEY,
    subject_id                TEXT NOT NULL REFERENCES subjects (subject_id),
    sample_type               TEXT NOT NULL,
    time_from_treatment_start INTEGER
);

-- The panel itself, so that population names are a foreign key rather than free text.
CREATE TABLE populations (
    population   TEXT PRIMARY KEY,
    display_name TEXT NOT NULL
);

-- The measurement fact table: one row per (sample, population).
CREATE TABLE cell_counts (
    sample_id  TEXT    NOT NULL REFERENCES samples (sample_id),
    population TEXT    NOT NULL REFERENCES populations (population),
    count      INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population)
);

-- Indexes chosen for the access patterns the analyses actually use: rolling samples
-- up to subjects and projects, and slicing subjects by cohort.
CREATE INDEX idx_subjects_project  ON subjects (project_id);
CREATE INDEX idx_subjects_cohort   ON subjects (condition, treatment, response);
CREATE INDEX idx_samples_subject   ON samples (subject_id);
CREATE INDEX idx_samples_slice     ON samples (sample_type, time_from_treatment_start);
CREATE INDEX idx_counts_population ON cell_counts (population);

-- Flattens the sample -> subject -> project chain so analytical queries express a
-- cohort filter without repeating the same two joins.
CREATE VIEW sample_details AS
SELECT
    s.sample_id,
    s.sample_type,
    s.time_from_treatment_start,
    sub.subject_id,
    sub.project_id,
    sub.condition,
    sub.age,
    sub.sex,
    sub.treatment,
    sub.response
FROM samples s
JOIN subjects sub ON sub.subject_id = s.subject_id;

-- Part 2 as a view, so the exported CSV, the dashboard and any ad-hoc SQL all share
-- one definition of "relative frequency" instead of reimplementing the arithmetic.
CREATE VIEW sample_frequencies AS
SELECT
    c.sample_id AS sample,
    SUM(c.count) OVER (PARTITION BY c.sample_id) AS total_count,
    c.population,
    c.count,
    100.0 * c.count / SUM(c.count) OVER (PARTITION BY c.sample_id) AS percentage
FROM cell_counts c;
