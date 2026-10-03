"""Calculate sample-level immune-cell frequencies from loblaw.db."""

import csv
import math
import os
import sqlite3
from collections import defaultdict
from contextlib import closing
from pathlib import Path
from typing import NamedTuple


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "loblaw.db"
OUTPUT_PATH = ROOT / "outputs" / "frequency_table.csv"
POPULATIONS = frozenset(("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"))


class FrequencyRow(NamedTuple):
    sample: str
    total_count: int
    population: str
    count: int
    percentage: float


def get_frequency_table(connection: sqlite3.Connection) -> list[FrequencyRow]:
    """Return one validated row per sample and immune-cell population."""
    sample_count = connection.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    if sample_count == 0:
        raise ValueError("The database contains no samples")

    query = """
        SELECT s.sample_id, c.cell_type, c.cell_count,
               SUM(c.cell_count) OVER (PARTITION BY s.sample_id) AS total_count
        FROM samples AS s
        JOIN cell_counts AS c ON c.sample_id = s.sample_id
        ORDER BY s.sample_id, c.cell_type
    """
    rows = []
    populations_by_sample = defaultdict(set)
    counts_by_sample = defaultdict(int)
    percentages_by_sample = defaultdict(float)

    for sample, population, count, total_count in connection.execute(query):
        if population not in POPULATIONS:
            raise ValueError(f"Unexpected population {population!r} for {sample}")
        if population in populations_by_sample[sample]:
            raise ValueError(f"Duplicate {sample} × {population} measurement")
        if count < 0:
            raise ValueError(f"Negative count for {sample} × {population}")
        if total_count <= 0:
            raise ValueError(f"Sample {sample} has zero total cell count")

        percentage = count / total_count * 100
        if not 0 <= percentage <= 100:
            raise ValueError(f"Invalid percentage for {sample} × {population}")
        rows.append(FrequencyRow(sample, total_count, population, count, percentage))
        populations_by_sample[sample].add(population)
        counts_by_sample[sample] += count
        percentages_by_sample[sample] += percentage

    if len(populations_by_sample) != sample_count or len(rows) != sample_count * len(POPULATIONS):
        raise ValueError("The frequency table does not cover every sample and population")
    for sample, populations in populations_by_sample.items():
        if populations != POPULATIONS:
            raise ValueError(f"Sample {sample} has an incomplete population set: {populations}")
    for row in rows:
        if row.total_count != counts_by_sample[row.sample]:
            raise ValueError(f"Incorrect total count for {row.sample}")
    for sample, percentage_sum in percentages_by_sample.items():
        if not math.isclose(percentage_sum, 100.0, rel_tol=0, abs_tol=1e-9):
            raise ValueError(f"Percentages do not sum to 100 for {sample}")
    return rows


def main():
    if not DB_PATH.is_file():
        raise FileNotFoundError(f"Database not found: {DB_PATH}. Run load_data.py first.")
    with closing(sqlite3.connect(DB_PATH)) as connection:
        rows = get_frequency_table(connection)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = OUTPUT_PATH.with_suffix(".csv.tmp")
    try:
        with temporary_path.open("w", newline="", encoding="utf-8") as target:
            writer = csv.writer(target)
            writer.writerow(FrequencyRow._fields)
            writer.writerows(rows)
        os.replace(temporary_path, OUTPUT_PATH)
    finally:
        temporary_path.unlink(missing_ok=True)
    print(f"Wrote {len(rows):,} population rows for "
          f"{len(rows) // len(POPULATIONS):,} samples to {OUTPUT_PATH}")
    print("\nSample frequency table (first 10 rows):")
    print(f"{'sample':<15} {'total_count':>11} {'population':<12} {'count':>9} {'percentage':>11}")
    for row in rows[:10]:
        print(f"{row.sample:<15} {row.total_count:>11} {row.population:<12} "
              f"{row.count:>9} {row.percentage:>10.2f}%")


if __name__ == "__main__":
    main()
