"""Transparent, research-only RTPLAN MLCX metrics. See docs/METRICS.md."""
from __future__ import annotations
from dataclasses import dataclass, asdict
import numpy as np

@dataclass(frozen=True)
class Aperture:
    """MLCX coordinates in mm at isocentre; boundaries are leaf Y edges."""
    left: np.ndarray
    right: np.ndarray
    boundaries: np.ndarray
    xjaw: tuple[float,float]
    yjaw: tuple[float,float]

    def __post_init__(self):
        if len(self.boundaries) < 2 or len(self.left) != len(self.right) or len(self.left) != len(self.boundaries)-1:
            raise ValueError("Inconsistent MLC leaf geometry")
        if not np.all(np.diff(self.boundaries) > 0):
            raise ValueError("Leaf position boundaries must be strictly increasing")
        if any(not np.all(np.isfinite(v)) for v in (self.left,self.right,self.boundaries)):
            raise ValueError("Non-finite MLC positions")
        if self.xjaw[0] >= self.xjaw[1] or self.yjaw[0] >= self.yjaw[1]:
            raise ValueError("Jaw limits must be increasing")

    def gaps(self) -> np.ndarray:
        return np.maximum(0., np.minimum(self.right,self.xjaw[1])-np.maximum(self.left,self.xjaw[0]))

    def exposed_heights(self) -> np.ndarray:
        return np.maximum(0.,np.minimum(self.boundaries[1:],self.yjaw[1])-np.maximum(self.boundaries[:-1],self.yjaw[0]))

    def area_mm2(self) -> float:
        return float(np.dot(self.gaps(),self.exposed_heights()))

@dataclass(frozen=True)
class BeamMetrics:
    beam_number: int
    beam_name: str
    control_points: int
    monitor_units: float | None
    aperture_area_mean_mm2: float | None
    aperture_area_cv: float | None
    leaf_gap_mean_mm: float | None
    leaf_gap_median_mm: float | None
    sas5_mu_weighted: float | None
    sas10_mu_weighted: float | None
    # Explicit v0.2 jaw-aware / MU-interval metrics (NOT universal reference values).
    mcs_vmat_jaw_clipped: float | None
    lsv_cp_mu_weighted: float | None
    aav_beam_cp_mu_weighted: float | None
    leaf_travel_total_mm: float
    leaf_travel_mean_per_transition_mm: float
    leaf_travel_max_per_transition_mm: float
    mean_cp_aperture_area_change_fraction: float | None
    notes: str = "Exploratory RTPLAN metrics; no patient-specific QA verdict."

    def to_dict(self) -> dict:
        return asdict(self)

def _weighted_sas(apertures: list[Aperture], weights: np.ndarray, threshold: float) -> float | None:
    """Fraction of OPEN jaw-exposed pairs below threshold; ending CP ΔCMW proxy."""
    num=den=0.
    for ap,w in zip(apertures,weights):
        gap=ap.gaps()
        opened=(gap>0)&(ap.exposed_heights()>0)
        num+=float(w*np.count_nonzero(opened&(gap<threshold)))
        den+=float(w*np.count_nonzero(opened))
    return num/den if den>0 else None

def _changes(aps: list[Aperture]) -> np.ndarray:
    """Sum of movements of BOTH banks including off-field pairs, each transition."""
    return np.array([float(np.abs(b.left-a.left).sum()+np.abs(b.right-a.right).sum())
                     for a,b in zip(aps[:-1],aps[1:])],dtype=float)


def control_point_complexity(apertures: list[Aperture]) -> tuple[np.ndarray, np.ndarray]:
    """Return per-CP LSV and beam-reference AAV for one MLCX beam.

    Research definition: LSV uses only OPEN jaw-exposed leaf pairs and
    adjacent positions within each bank; a single (or no) active pair
    has LSV=1 by convention. BeamAAV denominator is the jaw-clipped
    virtual union envelope of all control points, with leaf widths.
    This definition is intentionally *not* aperture_area_cv.
    """
    if not apertures:
        return np.array([], dtype=float), np.array([], dtype=float)
    left_env = np.min(np.stack([
        np.maximum(ap.left, ap.xjaw[0]) for ap in apertures]), axis=0)
    right_env = np.max(np.stack([
        np.minimum(ap.right, ap.xjaw[1]) for ap in apertures]), axis=0)
    height_env = np.max(np.stack([
        ap.exposed_heights() for ap in apertures]), axis=0)
    envelope = float(np.dot(np.maximum(0., right_env-left_env),height_env))
    lsv = []
    aav = []
    for ap in apertures:
        selected = (ap.exposed_heights() > 0) & (ap.gaps() > 0)
        bank_scores = []
        for bank in (ap.left[selected], ap.right[selected]):
            if bank.size < 2:
                bank_scores.append(1.)
                continue
            spread = float(np.ptp(bank))
            bank_scores.append(
                1. if spread == 0 else
                float(np.clip(1. - np.mean(np.abs(np.diff(bank)))/spread, 0., 1.)))
        lsv.append(bank_scores[0] * bank_scores[1])
        aav.append(ap.area_mm2() / envelope if envelope > 0 else np.nan)
    return np.asarray(lsv, dtype=float), np.asarray(aav, dtype=float)


def modulation_scores(
    apertures: list[Aperture], cmw: np.ndarray
) -> tuple[float | None, float | None, float | None]:
    """Return (MCSv, LSV, beamAAV) using trapezoidal CP values and ΔCMW."""
    if len(apertures) < 2:
        return None, None, None
    dm = np.diff(cmw)
    positive = dm > 0
    if not np.any(positive):
        return None, None, None
    lsv, aav = control_point_complexity(apertures)
    if not np.all(np.isfinite(aav)):
        return None, None, None
    lsv_intervals = (lsv[:-1] + lsv[1:]) / 2
    aav_intervals = (aav[:-1] + aav[1:]) / 2
    weights = dm[positive] / dm[positive].sum()
    l_mean = float(np.dot(weights, lsv_intervals[positive]))
    a_mean = float(np.dot(weights, aav_intervals[positive]))
    mcs = float(np.dot(weights,
                       lsv_intervals[positive] * aav_intervals[positive]))
    return mcs, l_mean, a_mean

def calculate_beam(beam_number:int, beam_name:str, apertures:list[Aperture],
                   cumulative_meterset:list[float], monitor_units:float|None) -> BeamMetrics:
    if not apertures:
        raise ValueError("No control points")
    a0=apertures[0]
    if any(len(a.left)!=len(a0.left) or len(a.right)!=len(a0.right) or not np.array_equal(a.boundaries,a0.boundaries)
           for a in apertures):
        raise ValueError("Inconsistent MLC geometry across control points")
    if len(cumulative_meterset)!=len(apertures):
        raise ValueError("Meterset weights must align with control points")
    cmw=np.asarray(cumulative_meterset,dtype=float)
    if not np.all(np.isfinite(cmw)) or np.any(cmw < -1e-8) or np.any(np.diff(cmw) < -1e-8):
        raise ValueError("Invalid or decreasing cumulative meterset weights")
    if monitor_units is not None and (not np.isfinite(monitor_units) or monitor_units<0):
        raise ValueError("Invalid beam monitor units")
    weights=np.r_[0.,np.maximum(0.,np.diff(cmw))]
    if len(apertures)==1 and cmw[0]>0: weights[0]=cmw[0]
    mcs, lsv_mean, aav_mean = modulation_scores(apertures, cmw)
    areas=np.asarray([a.area_mm2() for a in apertures],dtype=float)
    area_mean=float(np.average(areas,weights=weights)) if weights.sum()>0 else None
    cv=float(areas.std()/areas.mean()) if areas.mean()>0 else None
    g=np.concatenate([a.gaps()[a.exposed_heights()>0] for a in apertures])
    opened=g[g>0]
    changes=_changes(apertures)
    denom=areas[:-1]+areas[1:]
    ratio=np.divide(np.abs(np.diff(areas)),denom,out=np.full(len(denom),np.nan),where=denom>0)
    valid=ratio[np.isfinite(ratio)]
    return BeamMetrics(
        beam_number=int(beam_number),beam_name=str(beam_name),control_points=len(apertures),
        monitor_units=float(monitor_units) if monitor_units is not None else None,
        aperture_area_mean_mm2=area_mean,aperture_area_cv=cv,
        leaf_gap_mean_mm=float(opened.mean()) if len(opened) else None,
        leaf_gap_median_mm=float(np.median(opened)) if len(opened) else None,
        sas5_mu_weighted=_weighted_sas(apertures,weights,5.),
        sas10_mu_weighted=_weighted_sas(apertures,weights,10.),
        mcs_vmat_jaw_clipped=mcs,
        lsv_cp_mu_weighted=lsv_mean,
        aav_beam_cp_mu_weighted=aav_mean,
        leaf_travel_total_mm=float(changes.sum()),
        leaf_travel_mean_per_transition_mm=float(changes.mean()) if len(changes) else 0.,
        leaf_travel_max_per_transition_mm=float(changes.max()) if len(changes) else 0.,
        mean_cp_aperture_area_change_fraction=float(valid.mean()) if len(valid) else None
    )

def control_point_records(apertures:list[Aperture], cmw:list[float], angles:list[float|None]) -> list[dict]:
    """Anonymous CP geometries for interactive visualization."""
    travels=_changes(apertures)
    lsv, aav = control_point_complexity(apertures)
    return [{
        "cp_index":i, "gantry_deg":angles[i],
        "cumulative_meterset_weight":float(cmw[i]),
        "aperture_area_mm2":ap.area_mm2(),
        "lsv_cp":float(lsv[i]),
        "aav_beam_cp":float(aav[i]) if np.isfinite(aav[i]) else None,
        "travel_from_previous_cp_mm":0. if i==0 else float(travels[i-1]),
        "left_mm":ap.left.tolist(),"right_mm":ap.right.tolist(),
        "leaf_boundaries_mm":ap.boundaries.tolist(),
        "xjaw_mm":list(ap.xjaw),"yjaw_mm":list(ap.yjaw)
    } for i,ap in enumerate(apertures)]
