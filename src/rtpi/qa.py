"""Private local PTW VeriSoft PDF/XCC QA extraction; no cloud or identifiers."""
from __future__ import annotations
import base64
from dataclasses import asdict, dataclass
from pathlib import Path
import re
import xml.etree.ElementTree as ET

def _num(value: str) -> float:
    return float(value.replace(",", "."))

def _integer(value: str) -> int:
    return int(re.sub(r"\s+", "", value))

@dataclass(frozen=True)
class GammaReport:
    gamma_pass_percent: float | None = None
    gamma_fail_percent: float | None = None
    evaluated_points: int | None = None
    passed_points: int | None = None
    failed_points: int | None = None
    dose_difference_percent: float | None = None
    dta_mm: float | None = None
    low_dose_threshold_percent: float | None = None
    low_dose_reference: str | None = None
    gamma_passing_limit: float | None = None
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict:
        return asdict(self)

def parse_verisoft_pdf_text(text: str) -> GammaReport:
    """Extract displayed gamma metrics; ambiguities are flagged, not resolved."""
    if not text or not text.strip():
        return GammaReport(warnings=("pdf_text_empty_needs_manual_review_or_ocr",))
    flags = []
    def percent_row(name: str) -> tuple[int | None, float | None]:
        expr = rf"(?im)^\s*{name}\s+([0-9][0-9 ]*)\s*\(\s*([0-9]+(?:[.,][0-9]+)?)\s*%\s*\)"
        values = {(_integer(a), _num(b)) for a,b in re.findall(expr,text)}
        if len(values)>1:
            flags.append("multiple_"+name.lower()+"_results")
        return next(iter(values)) if len(values)==1 else (None,None)
    passed,pass_pct=percent_row("Passed")
    failed,fail_pct=percent_row("Failed")
    evaluated={_integer(n) for n in re.findall(
        r"(?im)^\s*Evaluated\s+(?:Dose\s+Points|Voxels)\s+([0-9][0-9 ]*)\s*\(",text)}
    if len(evaluated)>1:
        flags.append("multiple_evaluated_point_totals")
    eval_count=next(iter(evaluated)) if len(evaluated)==1 else None
    def unique_num(expr,label):
        values={_num(s) for s in re.findall(expr,text,flags=re.I|re.M)}
        if len(values)>1: flags.append("multiple_"+label)
        return next(iter(values)) if len(values)==1 else None
    dta=unique_num(r"([0-9]+(?:[.,][0-9]+)?)\s*mm\s*Distance\s*[-–]?\s*To\s*[-–]?\s*Agreement","dta")
    dd=unique_num(r"([0-9]+(?:[.,][0-9]+)?)\s*%\s*Dose\s+difference","dose_difference")
    dose_lines=re.findall(r"(?im)^\s*Suppress\s+dose\s+below\s+([0-9]+(?:[.,][0-9]+)?)\s*%\s*([^\r\n]*)",text)
    threshold_set={_num(s) for s,_ in dose_lines}
    threshold=next(iter(threshold_set)) if len(threshold_set)==1 else None
    if len(threshold_set)>1: flags.append("multiple_low_dose_thresholds")
    modes=set()
    for _,s in dose_lines:
        ss=s.lower()
        modes.add("normalization_value" if "normalization value" in ss else
                  "calculated_volume_max" if "max. dose of calculated volume" in ss else
                  "other_or_unspecified")
    basis=next(iter(modes)) if len(modes)==1 else None
    if len(modes)>1: flags.append("multiple_low_dose_reference_modes")
    limit=unique_num(r"Gamma\s*(?:≤|<=|<)\s*([0-9]+(?:[.,][0-9]+)?)","gamma_passing_limit")
    if pass_pct is None: flags.append("gpr_not_extracted")
    if dta is None or dd is None or threshold is None:
        flags.append("gamma_parameters_incomplete")
    if eval_count is not None and passed is not None and failed is not None and passed+failed!=eval_count:
        flags.append("point_counts_inconsistent")
    if pass_pct is not None and fail_pct is not None and abs(pass_pct+fail_pct-100)>0.2:
        flags.append("percentage_totals_inconsistent")
    return GammaReport(pass_pct,fail_pct,eval_count,passed,failed,dd,dta,threshold,basis,limit,tuple(flags))

@dataclass(frozen=True)
class XccMetadata:
    format_version: str | None
    software_module: str | None
    detector_model: str | None
    storage_mode: str | None
    measurement_count: int
    detector_count: int | None
    values_per_measurement: int | None
    sample_interval_s: float | None
    last_sample_time_s: float | None
    warnings: tuple[str,...] = ()

    def to_dict(self) -> dict:
        return asdict(self)

def parse_xcc_metadata(source: str | Path) -> XccMetadata:
    """Parse XCC timing/shape; does not calculate gamma from detector readings."""
    source=Path(source)
    if source.stat().st_size>250_000_000:
        raise ValueError("XCC size exceeds local parser limit")
    with source.open("rb") as f:
        if b"<!DOCTYPE" in f.read(4096).upper():
            raise ValueError("DOCTYPE not accepted")
    root=ET.parse(source).getroot()
    if root.tag!="PTW":
        raise ValueError("Expected PTW XCC")
    count_text=root.findtext(".//DetectorNumbers")
    detector_count=None
    if count_text:
        x=base64.b64decode(count_text.strip(),validate=True)
        if len(x)%2: raise ValueError("Invalid DetectorNumbers array length")
        detector_count=len(x)//2
    measurements=root.findall(".//MeasData/Measurement")
    if not measurements: measurements=root.findall(".//Measurement")
    value_count=None
    last_time=None
    if measurements:
        first=measurements[0].findtext("Data")
        if first:
            blob=base64.b64decode(first.strip(),validate=True)
            if len(blob)%4: raise ValueError("Invalid float32 Data length")
            value_count=len(blob)//4
        s=measurements[-1].findtext("Time")
        last_time=_num(s) if s else None
    storage=root.findtext(".//DeviceStoreMode")
    flags=[]
    if not measurements: flags.append("xcc_has_no_measurements")
    if value_count is not None and detector_count is not None and value_count!=detector_count:
        flags.append("detector_value_count_mismatch")
    if storage!="DoseAccumulated": flags.append("xcc_storage_mode_not_dose_accumulated")
    interval=root.findtext(".//IntervalTime")
    return XccMetadata(
        root.findtext("Version"),root.findtext(".//Adminsitrative/Module"),
        root.findtext(".//Detector"),storage,len(measurements),detector_count,
        value_count,_num(interval) if interval else None,last_time,tuple(flags))
