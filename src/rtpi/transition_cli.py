"""Anonymized local RTPLAN CP-transition / VMAT-sector CSV exports."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from .dicom import read_rtplan
from .transitions import MachineLimits


def export_transition_tables(
    paths: list[Path], output_dir: Path, *, sector_width_deg: int = 30,
    machine_limits: MachineLimits | None = None,
) -> dict[str, int]:
    """No names, UIDs, original paths or patient metadata in resulting CSVs."""
    summaries, transitions, sectors = [], [], []
    for idx, path in enumerate(paths, 1):
        try:
            beams = read_rtplan(
                path, include_transitions=True, sector_width_deg=sector_width_deg,
                machine_limits=machine_limits,
            )
        except Exception as exc:
            # Avoid leaking exception messages that may contain source filenames.
            raise ValueError(
                f"RTPLAN #{idx}: parsing failed ({type(exc).__name__})"
            ) from None
        for beam in beams:
            keys={"plan_index":idx,"beam_number":beam["beam_number"]}
            summary=beam["transition_summary"]
            summaries.append({
                **keys,
                "beam_mu":beam["monitor_units"],
                "mcs_vmat_jaw_clipped":beam["mcs_vmat_jaw_clipped"],
                **summary,
            })
            transitions.extend({**keys,**r} for r in beam["transition_data"])
            sectors.extend({**keys,**r} for r in beam["sector_data"])
    output_dir.mkdir(parents=True,exist_ok=True)
    def write(name: str, rows: list[dict]):
        fields=["plan_index","beam_number"]
        fields.extend(sorted(set().union(*(set(row)-set(fields) for row in rows))))
        with (output_dir/name).open("w",encoding="utf-8",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    write("rtpi_beam_transition_summary.csv",summaries)
    write("rtpi_cp_transitions.csv",transitions)
    write("rtpi_vmat_sectors.csv",sectors)
    return {"beams":len(summaries),"transitions":len(transitions),"sectors":len(sectors)}


def main() -> None:
    parser=argparse.ArgumentParser(
        description="Anonymized technical CP/MLC research export, NOT clinical QA."
    )
    parser.add_argument("plans",nargs="+",type=Path,help="Local DICOM RTPLAN paths")
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--sector-width-deg",type=int,default=30)
    parser.add_argument("--max-leaf-speed-mm-s",type=float,default=None)
    parser.add_argument("--max-gantry-speed-deg-s",type=float,default=None)
    parser.add_argument("--max-dose-rate-mu-min",type=float,default=None)
    args=parser.parse_args()
    limits = None
    if any(x is not None for x in (
            args.max_leaf_speed_mm_s,args.max_gantry_speed_deg_s,
            args.max_dose_rate_mu_min)):
        limits=MachineLimits(
            max_leaf_speed_mm_s=args.max_leaf_speed_mm_s,
            max_gantry_speed_deg_s=args.max_gantry_speed_deg_s,
            max_dose_rate_mu_min=args.max_dose_rate_mu_min,
        )
    counts=export_transition_tables(
        args.plans,args.output_dir,sector_width_deg=args.sector_width_deg,
        machine_limits=limits,
    )
    print(
        f"Analyzed {counts['beams']} beams, {counts['transitions']} CP intervals, "
        f"{counts['sectors']} angular sectors. No clinical identifiers exported."
    )
    if limits is not None:
        print("Kinematic values are idealized lower bounds, NOT observed times or QA verdicts.")


if __name__=="__main__":
    main()
