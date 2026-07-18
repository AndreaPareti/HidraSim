#!/usr/bin/env python3
"""Map DREMTubes tower IDs to front-face tower labels.

The front-face label convention is T<column><row>, as in T311:
  - column is 1..5 from left to right as seen from the beam;
  - row is 01..20 from bottom to top.

The active map below mirrors the TB26 `modflag` geometry in
include/DREMTubesGeoPar.hh. The simulation grid column is mirrored when
converted to front-face labels because the figure is drawn as seen from
the incoming beam.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import asdict, dataclass


NOF_MODULES_X = 5
NOF_MODULES_Y = 20

MODFLAG = [
    -1, -1, 0, -1, -1,
    -1, 1, 2, 3, -1,
    -1, 4, 5, 6, -1,
    -1, 7, 8, 9, -1,
    10, 11, 12, 13, 14,
    15, 16, 17, 18, 19,
    20, 21, 22, 23, 24,
    25, 26, 27, 28, 29,
    30, 31, 32, 33, 34,
    35, 36, 37, 38, 39,
    40, 41, 42, 43, 44,
    45, 46, 47, 48, 49,
    50, 51, 52, 53, 54,
    55, 56, 57, 58, 59,
    60, 61, 62, 63, 64,
    65, 66, 67, 68, 69,
    -1, 70, 71, 72, -1,
    -1, 73, 74, 75, -1,
    -1, 76, 77, 78, -1,
    -1, -1, 79, -1, -1,
]


TUBE_RADIUS_MM = 1.0
SQRT3_APPROX = 1.733
NOF_FIBERS_COLUMN = 64
NOF_FIBERS_ROW = 16


@dataclass(frozen=True)
class TowerMapEntry:
    tower_id: int
    mapped_tower: str
    front_column: int
    front_row: int
    x_mm: float
    y_mm: float
    simulation_column: int
    simulation_row: int


def build_mapping() -> list[TowerMapEntry]:
    mapping: list[TowerMapEntry] = []
    pitch_x = 2.0 * TUBE_RADIUS_MM * NOF_FIBERS_COLUMN
    pitch_y = SQRT3_APPROX * TUBE_RADIUS_MM * NOF_FIBERS_ROW

    for index, tower_id in enumerate(MODFLAG):
        if tower_id < 0:
            continue

        simulation_row = index // NOF_MODULES_X
        simulation_column = index % NOF_MODULES_X

        front_column = NOF_MODULES_X - simulation_column
        front_row = simulation_row + 1
        mapped_tower = f"T{front_column}{front_row:02d}"
        x_mm = -pitch_x * simulation_column + pitch_x * (NOF_MODULES_X - 1) / 2.0
        y_mm = pitch_y * simulation_row - pitch_y * (NOF_MODULES_Y - 1) / 2.0

        mapping.append(
            TowerMapEntry(
                tower_id=tower_id,
                mapped_tower=mapped_tower,
                front_column=front_column,
                front_row=front_row,
                x_mm=x_mm,
                y_mm=y_mm,
                simulation_column=simulation_column,
                simulation_row=simulation_row,
            )
        )

    return sorted(mapping, key=lambda entry: entry.tower_id)


def mapping_by_tower_id() -> dict[int, TowerMapEntry]:
    return {entry.tower_id: entry for entry in build_mapping()}


def mapping_by_label() -> dict[str, TowerMapEntry]:
    return {entry.mapped_tower: entry for entry in build_mapping()}


def print_entry(entry: TowerMapEntry) -> None:
    print(
        f"{entry.tower_id} -> {entry.mapped_tower} "
        f"(front column {entry.front_column}, front row {entry.front_row}, "
        f"x={entry.x_mm:.3f} mm, y={entry.y_mm:.3f} mm)"
    )


def dump_table(entries: list[TowerMapEntry]) -> None:
    print("tower_id mapped_tower front_column front_row x_mm y_mm simulation_column simulation_row")
    for entry in entries:
        print(
            f"{entry.tower_id:8d} {entry.mapped_tower:12s} "
            f"{entry.front_column:12d} {entry.front_row:9d} "
            f"{entry.x_mm:8.3f} {entry.y_mm:8.3f} "
            f"{entry.simulation_column:17d} {entry.simulation_row:14d}"
        )


def dump_csv(entries: list[TowerMapEntry]) -> None:
    writer = csv.DictWriter(sys.stdout, fieldnames=list(asdict(entries[0]).keys()))
    writer.writeheader()
    for entry in entries:
        writer.writerow(asdict(entry))


def dump_json(entries: list[TowerMapEntry]) -> None:
    print(json.dumps([asdict(entry) for entry in entries], indent=2))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Map DREMTubes tower IDs to front-face labels such as T311."
    )
    parser.add_argument("--tower-id", type=int, help="Tower ID to map, for example 37.")
    parser.add_argument("--mapped-tower", help="Mapped tower label to map back, for example T311.")
    parser.add_argument(
        "--dump",
        choices=("table", "csv", "json"),
        help="Dump the full tower map.",
    )
    args = parser.parse_args()

    if args.tower_id is None and args.mapped_tower is None and args.dump is None:
        parser.print_help()
        return 0

    entries = build_mapping()

    if args.tower_id is not None:
        by_tower = {entry.tower_id: entry for entry in entries}
        if args.tower_id not in by_tower:
            print(f"Unknown tower ID: {args.tower_id}", file=sys.stderr)
            return 1
        print_entry(by_tower[args.tower_id])

    if args.mapped_tower is not None:
        mapped_tower = args.mapped_tower.upper()
        by_label = {entry.mapped_tower: entry for entry in entries}
        if mapped_tower not in by_label:
            print(f"Unknown mapped tower: {args.mapped_tower}", file=sys.stderr)
            return 1
        entry = by_label[mapped_tower]
        print(f"{entry.mapped_tower} -> {entry.tower_id}")

    if args.dump == "table":
        dump_table(entries)
    elif args.dump == "csv":
        dump_csv(entries)
    elif args.dump == "json":
        dump_json(entries)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
