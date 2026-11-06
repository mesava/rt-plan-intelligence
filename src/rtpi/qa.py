"""Privacy-minimizing QA ingestion. PDF text extraction is not OCR.

Only explicit Gamma QA findings are extracted; all patient identifiers are
ignored. The source CSV may still contain PHI: process it on an authorized
local workstation, never as part of the public demonstration.
"""
from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path
import re
import statistics
from typing import TextIO


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
