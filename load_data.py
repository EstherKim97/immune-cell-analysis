"""Load cell-count.csv into a normalized SQLite database beside this script."""

import csv
import os
import sqlite3
from contextlib import closing
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "loblaw.db"
CELL_TYPES = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")
FIELDS = (
    "project", "subject", "condition", "age", "sex", "treatment", "response",
    "sample", "sample_type", "time_from_treatment_start", *CELL_TYPES,
)

SCHEMA = """
CREATE TABLE projects (
    project_id TEXT PRIMARY KEY NOT NULL CHECK (length(project_id) > 0)
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY NOT NULL CHECK (length(subject_id) > 0),
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    condition TEXT NOT NULL CHECK (length(condition) > 0),
    age INTEGER NOT NULL CHECK (age BETWEEN 0 AND 120),
    sex TEXT NOT NULL CHECK (length(sex) > 0),
    treatment TEXT NOT NULL CHECK (length(treatment) > 0),
    response TEXT CHECK (response IN ('yes', 'no'))
);

CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY NOT NULL CHECK (length(sample_id) > 0),
    subject_id TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type TEXT NOT NULL CHECK (length(sample_type) > 0),
    time_from_treatment_start INTEGER NOT NULL,
    UNIQUE (subject_id, sample_type, time_from_treatment_start)
);

CREATE TABLE cell_counts (
    sample_id TEXT NOT NULL REFERENCES samples(sample_id),
    cell_type TEXT NOT NULL CHECK (cell_type IN
        ('b_cell', 'cd8_t_cell', 'cd4_t_cell', 'nk_cell', 'monocyte')),
    cell_count INTEGER NOT NULL CHECK (cell_count >= 0),
    PRIMARY KEY (sample_id, cell_type)
) WITHOUT ROWID;

CREATE INDEX idx_subjects_project ON subjects(project_id);
CREATE INDEX idx_subjects_cohort ON subjects(condition, treatment, response);
CREATE INDEX idx_samples_subject ON samples(subject_id);
CREATE INDEX idx_samples_type_time ON samples(sample_type, time_from_treatment_start);
"""


def required_text(row, field, line_number):
    value = row[field].strip()
    if not value:
        raise ValueError(f"Line {line_number}: {field} is empty")
    return value


def integer(row, field, line_number):
    value = required_text(row, field, line_number)
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(f"Line {line_number}: invalid integer for {field}: {value!r}") from exc


def load(connection):
    projects = set()
    subjects = {}
    samples = []
    counts = []
    sample_ids = set()

    with CSV_PATH.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames != list(FIELDS):
            raise ValueError(f"Unexpected CSV columns: {reader.fieldnames!r}")

        for line_number, row in enumerate(reader, start=2):
            project_id = required_text(row, "project", line_number)
            subject_id = required_text(row, "subject", line_number)
            sample_id = required_text(row, "sample", line_number)
            response = row["response"].strip() or None
            if response not in (None, "yes", "no"):
                raise ValueError(f"Line {line_number}: invalid response: {response!r}")

            subject = (
                subject_id, project_id,
                required_text(row, "condition", line_number),
                integer(row, "age", line_number),
                required_text(row, "sex", line_number),
                required_text(row, "treatment", line_number), response,
            )
            if subject_id in subjects and subjects[subject_id] != subject:
                raise ValueError(f"Line {line_number}: conflicting metadata for {subject_id}")
            if sample_id in sample_ids:
                raise ValueError(f"Line {line_number}: duplicate sample {sample_id}")

            projects.add(project_id)
            subjects[subject_id] = subject
            sample_ids.add(sample_id)
            samples.append((
                sample_id, subject_id,
                required_text(row, "sample_type", line_number),
                integer(row, "time_from_treatment_start", line_number),
            ))
            for cell_type in CELL_TYPES:
                count = integer(row, cell_type, line_number)
                if count < 0:
                    raise ValueError(f"Line {line_number}: negative {cell_type} count")
                counts.append((sample_id, cell_type, count))

    if not samples:
        raise ValueError("The CSV contains no samples")

    connection.executemany("INSERT INTO projects VALUES (?)", [(value,) for value in sorted(projects)])
    connection.executemany("INSERT INTO subjects VALUES (?, ?, ?, ?, ?, ?, ?)",
                           [subjects[key] for key in sorted(subjects)])
    connection.executemany("INSERT INTO samples VALUES (?, ?, ?, ?)", samples)
    connection.executemany("INSERT INTO cell_counts VALUES (?, ?, ?)", counts)
    return len(projects), len(subjects), len(samples), len(counts)


def main():
    temporary_path = DB_PATH.with_suffix(".db.tmp")
    temporary_path.unlink(missing_ok=True)
    try:
        with closing(sqlite3.connect(temporary_path)) as connection:
            with connection:
                connection.execute("PRAGMA foreign_keys = ON")
                connection.executescript(SCHEMA)
                totals = load(connection)
                if connection.execute("PRAGMA foreign_key_check").fetchall():
                    raise ValueError("Foreign key validation failed")
        os.replace(temporary_path, DB_PATH)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise
    print(f"Loaded {totals[0]} projects, {totals[1]} subjects, "
          f"{totals[2]} samples, and {totals[3]} cell counts into {DB_PATH}")


if __name__ == "__main__":
    main()
