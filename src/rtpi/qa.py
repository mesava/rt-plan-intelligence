"""Private local PTW VeriSoft PDF/XCC QA extraction; no cloud or identifiers."""
from __future__ import annotations
import base64
import csv
import statistics
from typing import TextIO
from dataclasses import asdict, dataclass
from pathlib import Path
import re
from defusedxml import ElementTree as ET

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


# Backward-compatible public research QA interface already used in the project.
@dataclass(frozen=True)
class QAReport:
    gamma_pass_percent: float | None
    dose_difference_percent: float | None
    dta_mm: float | None
    low_dose_threshold_percent: float | None
    passed_count: int | None
    failed_count: int | None
    analysis_kind: str | None
    notes: str = "Parsed QA text; independently verify criteria and normalization."

    def to_dict(self) -> dict:
        return asdict(self)


def _decimal(s: str | None) -> float | None:
    return float(s.replace(",", ".")) if s is not None else None


def _match_number(pattern: str, text: str) -> float | None:
    found = re.search(pattern, text, re.IGNORECASE | re.MULTILINE)
    return _decimal(found.group(1)) if found else None


def _count(label: str, text: str) -> int | None:
    found = re.search(
        rf"^\s*{label}\s+([\d ]+)\s*\(\s*\d+[.,]\d+\s*%\s*\)",
        text, re.IGNORECASE | re.MULTILINE
    )
    return int(found.group(1).replace(" ", "")) if found else None


def parse_qa_text(text: str) -> QAReport:
    """Parse human-readable PTW-style PDF report text without patient fields."""
    gamma = _match_number(r"\bResult\s+(\d+(?:[,.]\d+)?)\s*%", text)
    difference = _match_number(r"(\d+(?:[,.]\d+)?)\s*%\s*Dose\s+difference", text)
    distance = _match_number(r"(\d+(?:[,.]\d+)?)\s*mm\s*Distance\s*-\s*To\s*-\s*Agreement", text)
    threshold = _match_number(r"Suppress\s+dose\s+below\s+(\d+(?:[,.]\d+)?)\s*%", text)
    kind = "Gamma 3D" if re.search(r"Gamma\s*3D", text, re.IGNORECASE) else (
        "Volume Analysis" if re.search(r"Volume\s+Analysis", text, re.IGNORECASE) else None
    )
    if gamma is not None and not 0 <= gamma <= 100:
        raise ValueError("Gamma pass percent must be between 0 and 100")
    return QAReport(gamma, difference, distance, threshold,
                    _count("Passed", text), _count("Failed", text), kind)


def read_qa_pdf(source: str | Path) -> QAReport:
    """Extract text using pypdf, without OCR and without patient identifiers."""
    from pypdf import PdfReader
    reader = PdfReader(str(source))
    if len(reader.pages) > 200:
        raise ValueError("Unexpectedly long QA report PDF")
    text = "\n".join((page.extract_text() or "") for page in reader.pages)
    return parse_qa_text(text)


def cohort_summary_csv(source: str | Path | TextIO) -> dict:
    """Summarize existing QA CSV; never return patient ID/name or row details."""
    if isinstance(source, (str, Path)):
        with open(source, "r", newline="", encoding="utf-8-sig") as stream:
            return cohort_summary_csv(stream)
    reader = csv.DictReader(source)
    required = {"case", "gamma_result_percent"}
    if not required.issubset(set(reader.fieldnames or [])):
        raise ValueError("CSV must contain case and gamma_result_percent")
    rows = list(reader)
    case_ids = []
    gamma = []
    with_dta = with_diff = with_threshold = 0
    for row in rows:
        try:
            case = int(row["case"])
            g = _decimal(row.get("gamma_result_percent", "").strip() or None)
        except (TypeError, ValueError) as exc:
            raise ValueError("Invalid case or gamma value") from exc
        if g is not None and not 0 <= g <= 100:
            raise ValueError("Gamma percent must be between 0 and 100")
        case_ids.append(case)
        if g is not None:
            gamma.append(g)
        with_dta += bool(row.get("dta_mm", "").strip())
        with_diff += bool(row.get("dose_diff_percent", "").strip())
        with_threshold += bool(row.get("threshold_percent", "").strip())
    if len(case_ids) != len(set(case_ids)):
        raise ValueError("Duplicate QA case index")
    return {
        "plans": len(rows),
        "gamma_observed": len(gamma),
        "gamma_mean_percent": float(statistics.mean(gamma)) if gamma else None,
        "gamma_median_percent": float(statistics.median(gamma)) if gamma else None,
        "gamma_below_95_count": sum(v < 95 for v in gamma),
        "explicit_dta_count": with_dta,
        "explicit_dose_difference_count": with_diff,
        "explicit_low_dose_threshold_count": with_threshold,
        "warning": (
            "Existing CSV values are not independently verified against PDFs. "
            "Mixed/missing gamma criteria and repeated patients may prevent "
            "valid pooled correlation analysis."
        ),
    }
