"""Part 4 database queries and subset summaries."""

import csv
import sqlite3
from collections import Counter
from contextlib import closing
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "loblaw.db"
OUTPUT_DIR = ROOT / "outputs"

BASELINE_MALE_B_CELL_FROM = """
    FROM subjects AS u
    JOIN samples AS s ON s.subject_id = u.subject_id
    JOIN cell_counts AS c ON c.sample_id = s.sample_id
    WHERE u.condition = 'melanoma'
      AND u.sex = 'M'
      AND u.response = 'yes'
      AND s.time_from_treatment_start = 0
      AND c.cell_type = 'b_cell'
"""


def get_baseline_cohort(connection: sqlite3.Connection) -> list[dict]:
    """Return samples matching all four Part 4 filters."""
    cursor = connection.execute(
        """
        SELECT s.sample_id AS sample, s.subject_id AS subject,
               u.project_id AS project, u.response, u.sex,
               u.condition, u.treatment, s.sample_type,
               s.time_from_treatment_start
        FROM samples AS s
        JOIN subjects AS u ON u.subject_id = s.subject_id
        WHERE u.condition = 'melanoma'
          AND s.sample_type = 'PBMC'
          AND u.treatment = 'miraclib'
          AND s.time_from_treatment_start = 0
        ORDER BY s.sample_id
        """
    )
    fields = [column[0] for column in cursor.description]
    return [dict(zip(fields, row)) for row in cursor]


def get_summary_tables(cohort: list[dict]) -> dict[str, list[dict]]:
    """Count samples by project and distinct subjects by response and sex."""
    if not cohort:
        raise ValueError("No samples match the Part 4 filters")

    sample_ids = set()
    subject_details = {}
    project_counts = Counter()
    for row in cohort:
        if (row["condition"], row["sample_type"], row["treatment"],
                row["time_from_treatment_start"]) != ("melanoma", "PBMC", "miraclib", 0):
            raise ValueError(f"Sample {row['sample']} does not match all Part 4 filters")
        if row["sample"] in sample_ids:
            raise ValueError(f"Duplicate sample {row['sample']}")
        if row["response"] not in ("yes", "no") or row["sex"] not in ("M", "F"):
            raise ValueError(f"Unexpected response or sex for subject {row['subject']}")
        details = (row["project"], row["response"], row["sex"])
        if row["subject"] in subject_details and subject_details[row["subject"]] != details:
            raise ValueError(f"Conflicting metadata for subject {row['subject']}")
        sample_ids.add(row["sample"])
        subject_details[row["subject"]] = details
        project_counts[row["project"]] += 1

    response_counts = Counter(details[1] for details in subject_details.values())
    sex_counts = Counter(details[2] for details in subject_details.values())
    if (sum(project_counts.values()) != len(sample_ids)
            or sum(response_counts.values()) != len(subject_details)
            or sum(sex_counts.values()) != len(subject_details)):
        raise ValueError("Part 4 counts do not reconcile with the filtered cohort")

    return {
        "samples_by_project": [
            {"project": project, "sample_count": project_counts[project]}
            for project in sorted(project_counts)
        ],
        "subjects_by_response": [
            {"response": response, "subject_count": response_counts[response]}
            for response in ("yes", "no")
        ],
        "subjects_by_sex": [
            {"sex": sex, "subject_count": sex_counts[sex]}
            for sex in ("M", "F")
        ],
    }


def get_baseline_male_responder_b_cell_average(connection: sqlite3.Connection) -> dict:
    """Average raw B-cell counts across all sample and treatment types."""
    selected = connection.execute(
        """SELECT s.sample_id, s.subject_id, u.condition, u.sex, u.response,
                  s.time_from_treatment_start, s.sample_type, u.treatment,
                  c.cell_type, c.cell_count
        """ + BASELINE_MALE_B_CELL_FROM
    ).fetchall()
    if not selected:
        raise ValueError("No B-cell observations match the requested cohort")

    sample_ids = set()
    subject_ids = set()
    sample_types = set()
    treatments = set()
    for (sample, subject, condition, sex, response, timepoint,
         sample_type, treatment, cell_type, count) in selected:
        if (condition, sex, response, timepoint, cell_type) != ("melanoma", "M", "yes", 0, "b_cell"):
            raise ValueError(f"Sample {sample} does not match the requested filters")
        if sample in sample_ids:
            raise ValueError(f"Duplicate B-cell observation for sample {sample}")
        if count is None or count < 0:
            raise ValueError(f"Invalid raw B-cell count for sample {sample}")
        sample_ids.add(sample)
        subject_ids.add(subject)
        sample_types.add(sample_type)
        treatments.add(treatment)

    aggregate = connection.execute(
        "SELECT AVG(c.cell_count), COUNT(DISTINCT s.sample_id), "
        "COUNT(DISTINCT s.subject_id) " + BASELINE_MALE_B_CELL_FROM
    ).fetchone()
    eligible_samples = connection.execute(
        """SELECT COUNT(DISTINCT s.sample_id)
           FROM subjects AS u JOIN samples AS s ON s.subject_id = u.subject_id
           WHERE u.condition = 'melanoma' AND u.sex = 'M'
             AND u.response = 'yes' AND s.time_from_treatment_start = 0"""
    ).fetchone()[0]
    if aggregate[1] != len(sample_ids) or aggregate[2] != len(subject_ids):
        raise ValueError("B-cell observation counts do not match distinct database IDs")
    if eligible_samples != len(sample_ids):
        raise ValueError("A qualifying sample is missing its B-cell measurement")

    return {
        "average_b_cell_count": aggregate[0],
        "qualifying_samples": len(sample_ids),
        "qualifying_subjects": len(subject_ids),
        "sample_types": sorted(sample_types),
        "treatments": sorted(treatments),
    }


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if not DB_PATH.is_file():
        raise FileNotFoundError(f"Database not found: {DB_PATH}. Run load_data.py first.")
    with closing(sqlite3.connect(DB_PATH.as_uri() + "?mode=ro", uri=True)) as connection:
        cohort = get_baseline_cohort(connection)
        b_cell_average = get_baseline_male_responder_b_cell_average(connection)
    tables = get_summary_tables(cohort)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUTPUT_DIR / "part4_cohort.csv", tuple(cohort[0]), cohort)
    write_csv(OUTPUT_DIR / "part4_samples_by_project.csv",
              ("project", "sample_count"), tables["samples_by_project"])
    write_csv(OUTPUT_DIR / "part4_subjects_by_response.csv",
              ("response", "subject_count"), tables["subjects_by_response"])
    write_csv(OUTPUT_DIR / "part4_subjects_by_sex.csv",
              ("sex", "subject_count"), tables["subjects_by_sex"])

    print(f"Filtered cohort: {len(cohort)} samples, "
          f"{len({row['subject'] for row in cohort})} subjects")
    for name, rows in tables.items():
        print(f"{name}: {rows}")
    print("\nMelanoma Male Responders at Baseline")
    print("All sample types and all treatment types")
    print(f"Qualifying samples: {b_cell_average['qualifying_samples']}")
    print(f"Qualifying subjects: {b_cell_average['qualifying_subjects']}")
    print(f"Average B-cell count: {b_cell_average['average_b_cell_count']:.2f}")
    print(f"Sample types represented: {', '.join(b_cell_average['sample_types'])}")
    print(f"Treatment types represented: {', '.join(b_cell_average['treatments'])}")


if __name__ == "__main__":
    main()
