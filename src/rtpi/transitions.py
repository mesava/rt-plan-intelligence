"""Research-only RTPLAN CP-to-CP geometry and optional kinematic lower bounds.

No log timestamps, actual leaf velocities, measured fluence, gamma or dose
agreement are inferred. Output is numerical and contains no patient identifiers.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from math import floor, isfinite
from typing import Sequence

import numpy as np

from .metrics import Aperture


@dataclass(frozen=True)
class MachineLimits:
    """Optional VERIFIED user-supplied maxima; not embedded machine specifications."""
    max_leaf_speed_mm_s: float | None = None
    max_gantry_speed_deg_s: float | None = None
    max_dose_rate_mu_min: float | None = None

    def __post_init__(self):
        for val in (self.max_leaf_speed_mm_s, self.max_gantry_speed_deg_s,
                    self.max_dose_rate_mu_min):
            if val is not None and (not isfinite(val) or val <= 0):
                raise ValueError("Machine limit must be strictly positive and finite")


@dataclass(frozen=True)
class TransitionAnalysis:
    summary: dict
    transitions: list[dict]
    sectors: list[dict]

    def to_dict(self) -> dict:
        return asdict(self)


def _signed_short_arc_deg(a: float | None, b: float | None) -> float | None:
    """Minimal signed angular change, not a DICOM gantry direction resolver."""
    if a is None or b is None:
        return None
    if not isfinite(a) or not isfinite(b):
        raise ValueError("Non-finite gantry angle")
    return ((b - a + 180.) % 360.) - 180.


def analyze_transitions(
    apertures: Sequence[Aperture],
    cumulative_meterset: Sequence[float],
    gantry_angles: Sequence[float | None],
    beam_meterset_mu: float | None,
    *,
    final_cmw: float | None = None,
    sector_width_deg: int = 30,
    limits: MachineLimits | None = None,
) -> TransitionAnalysis:
    """Inspect each CP interval, then aggregate rotating intervals into sectors.

    ΔMU = BeamMeterset × ΔCMW/FinalCumulativeMetersetWeight.
    Gantry sectors assign a whole interval to its angular midpoint: this is
    deliberately an approximation, not a dose influence map or fluence model.
    Max leaf displacement considers BOTH banks, including closed/off-jaw pairs.
    An additional active-leaf displacement uses jaw-exposed OPEN pairs.
    """
    n = len(apertures)
    if not (len(cumulative_meterset) == len(gantry_angles) == n):
        raise ValueError("CP positions, CMW and gantry angles must align")
    if sector_width_deg < 1 or 360 % sector_width_deg:
        raise ValueError("sector_width_deg must be a positive divisor of 360")
    if limits is not None and not isinstance(limits, MachineLimits):
        raise TypeError("limits must be MachineLimits or None")
    cmw = np.asarray(cumulative_meterset, dtype=float)
    if not np.all(np.isfinite(cmw)) or np.any(cmw < -1e-8) or np.any(np.diff(cmw) < -1e-8):
        raise ValueError("Invalid cumulative meterset weights")
    if beam_meterset_mu is not None and (
        not isfinite(beam_meterset_mu) or beam_meterset_mu < 0
    ):
        raise ValueError("BeamMeterset must be finite and nonnegative")
    end = float(cmw[-1]) if n and final_cmw is None else final_cmw
    if end is not None and (not isfinite(end) or end < 0):
        raise ValueError("Invalid FinalCumulativeMetersetWeight")
    if end is not None and n and abs(float(cmw[-1]) - end) > 1e-5:
        raise ValueError("Final CMW differs from last CP")
    mu_known = beam_meterset_mu is not None and end is not None and end > 0
    rows: list[dict] = []
    aggregated: dict[int, dict] = {}
    for i in range(1, n):
        a, b = apertures[i-1], apertures[i]
        if len(a.left) != len(b.left) or not np.array_equal(a.boundaries, b.boundaries):
            raise ValueError("Inconsistent MLC geometry across CP")
        dl = np.abs(b.left-a.left)
        dr = np.abs(b.right-a.right)
        max_step = float(max(np.max(dl), np.max(dr)))
        sum_step = float(np.sum(dl) + np.sum(dr))
        active = (((a.exposed_heights() > 0) & (a.gaps() > 0)) |
                  ((b.exposed_heights() > 0) & (b.gaps() > 0)))
        max_active = (float(max(np.max(dl[active]), np.max(dr[active])))
                      if np.any(active) else None)
        dcmw = max(0., float(cmw[i]-cmw[i-1]))
        delta_mu = (float(beam_meterset_mu*dcmw/end) if mu_known else None)
        signed_angle = _signed_short_arc_deg(gantry_angles[i-1], gantry_angles[i])
        sweep = abs(signed_angle) if signed_angle is not None else None
        midpoint = ((float(gantry_angles[i-1]) + signed_angle/2.) % 360.
                    if signed_angle is not None else None)
        mm_mu = max_step/delta_mu if delta_mu is not None and delta_mu > 0 else None
        mm_degree = max_step/sweep if sweep is not None and sweep > 0 else None
        flags = []
        if dcmw == 0 and max_step > 1e-8:
            flags.append("leaf_motion_without_incremental_mu")
        if dcmw == 0 and sweep is not None and sweep > 1e-8:
            flags.append("gantry_motion_without_incremental_mu")
        if signed_angle is None:
            flags.append("gantry_angle_unavailable")
        if signed_angle is not None and abs(signed_angle) >= 180.-1e-7:
            flags.append("gantry_shortest_path_ambiguous")
        if dcmw > 0 and not mu_known:
            flags.append("absolute_mu_unavailable")
        bound = None
        if limits is not None:
            components = []
            if limits.max_leaf_speed_mm_s is not None:
                components.append(max_step/limits.max_leaf_speed_mm_s)
            if limits.max_gantry_speed_deg_s is not None and sweep is not None:
                components.append(sweep/limits.max_gantry_speed_deg_s)
            if limits.max_dose_rate_mu_min is not None and delta_mu is not None:
                components.append(60.*delta_mu/limits.max_dose_rate_mu_min)
            if components:
                bound = float(max(components))
        row = {
            "cp_from":i-1, "cp_to":i,
            "delta_cmw":dcmw,
            "delta_mu":delta_mu,
            "gantry_angle_from_deg":(float(gantry_angles[i-1])
                                     if gantry_angles[i-1] is not None else None),
            "gantry_angle_to_deg":(float(gantry_angles[i])
                                   if gantry_angles[i] is not None else None),
            "gantry_sweep_shortest_deg":sweep,
            "gantry_midpoint_deg":midpoint,
            "max_single_leaf_step_mm":max_step,
            "max_active_leaf_step_mm":max_active,
            "total_leaf_step_mm":sum_step,
            "max_leaf_mm_per_mu":mm_mu,
            "max_leaf_mm_per_degree":mm_degree,
            "kinematic_min_seconds_from_limits":bound,
            "flags":";".join(flags),
        }
        rows.append(row)
        if midpoint is not None and sweep is not None and sweep > 1e-8:
            start = int(floor(midpoint / sector_width_deg)) * sector_width_deg
            if start not in aggregated:
                aggregated[start] = {
                    "sector_start_deg":start,
                    "sector_end_deg":start+sector_width_deg,
                    "transition_count":0,
                    "sum_delta_mu":0. if mu_known else None,
                    "sum_gantry_sweep_deg":0.,
                    "sum_leaf_step_mm":0.,
                    "max_single_leaf_step_mm":0.,
                    "zero_mu_leaf_motion_count":0,
                    "max_leaf_mm_per_mu":None,
                }
            sec=aggregated[start]
            sec["transition_count"] += 1
            if delta_mu is not None:
                sec["sum_delta_mu"] += delta_mu
            sec["sum_gantry_sweep_deg"] += sweep
            sec["sum_leaf_step_mm"] += sum_step
            sec["max_single_leaf_step_mm"] = max(sec["max_single_leaf_step_mm"],max_step)
            if dcmw == 0 and max_step > 1e-8:
                sec["zero_mu_leaf_motion_count"] += 1
            if mm_mu is not None:
                sec["max_leaf_mm_per_mu"] = max(
                    sec["max_leaf_mm_per_mu"] if sec["max_leaf_mm_per_mu"] is not None else 0.,
                    mm_mu,
                )
    moving = [r for r in rows if r["gantry_sweep_shortest_deg"] is not None
              and r["gantry_sweep_shortest_deg"] > 1e-8]
    all_mu=[r["delta_mu"] for r in rows]
    summary = {
        "transition_count":len(rows),
        "rotating_transition_count":len(moving),
        "max_single_leaf_step_mm":max((r["max_single_leaf_step_mm"] for r in rows),default=None),
        "max_active_leaf_step_mm":max((r["max_active_leaf_step_mm"] for r in rows
                                     if r["max_active_leaf_step_mm"] is not None),default=None),
        "largest_max_leaf_mm_per_mu":max((r["max_leaf_mm_per_mu"] for r in rows
                                        if r["max_leaf_mm_per_mu"] is not None),default=None),
        "zero_mu_leaf_motion_count":sum("leaf_motion_without_incremental_mu" in r["flags"]
                                         for r in rows),
        "total_delta_mu":sum(all_mu) if all(x is not None for x in all_mu) else None,
        "gantry_sweep_shortest_total_deg":sum(
            r["gantry_sweep_shortest_deg"] for r in moving),
        "kinematic_min_total_seconds":(
            sum(r["kinematic_min_seconds_from_limits"] for r in rows)
            if rows and all(r["kinematic_min_seconds_from_limits"] is not None
                            for r in rows) else None
        ),
        "sector_width_deg":sector_width_deg,
        "angular_assignment":"interval_midpoint_shortest_path",
        "time_axis_available":False,
        "measured_deliverability_or_gamma":False,
    }
    sectors=sorted(aggregated.values(),key=lambda x:x["sector_start_deg"])
    for sec in sectors:
        mu=sec["sum_delta_mu"]
        deg=sec["sum_gantry_sweep_deg"]
        sec["mu_per_swept_degree"] = mu/deg if mu is not None and deg>0 else None
        sec["leaf_step_mm_per_mu"] = (sec["sum_leaf_step_mm"]/mu
                                       if mu is not None and mu>0 else None)
    return TransitionAnalysis(summary, rows, sectors)
