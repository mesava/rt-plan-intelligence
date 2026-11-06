"""Read PTW XCC measurement structure without exposing patient identity or detector arrays.

Research-only: XCC frames are accumulated detector readings, NOT gamma, film
fluence, a dose-to-water map, or a validated reconstructed dose distribution.
"""
from __future__ import annotations

import base64
from dataclasses import asdict, dataclass
import math
from pathlib import Path
import statistics
from typing import BinaryIO
from defusedxml import ElementTree as SafeET


@dataclass(frozen=True)
class XCCMetadata:
    frame_count: int
    detectors_per_frame: int
    start_time_s: float
    last_time_s: float
    median_interval_s: float | None
    first_angle_deg: float
    last_angle_deg: float
    minimum_angle_deg: float
    maximum_angle_deg: float
    measurement_preset_s: float | None
    device_store_mode: str | None
    matrix_resolution_lr_mm: float | None
    matrix_resolution_gt_mm: float | None
    notes: str = (
        "PTW XCC structural metadata only; no gamma, clinical QA verdict "
        "or validated dose/fluence reconstruction."
    )

    def to_dict(self) -> dict:
        return asdict(self)


def _optional_float(text: str | None) -> float | None:
    if not text or not text.strip():
        return None
    v = float(text.strip().replace(",", "."))
    if not math.isfinite(v):
        raise ValueError("Non-finite numeric value in XCC")
    return v


def read_xcc_metadata(source: str | Path | BinaryIO) -> XCCMetadata:
    """Parse PTW XML/XCC with DTD/entity-safe XML parsing.

    Frame payloads are checked as base64 float32-little-endian blobs and
    their detector counts must match. Only structural metadata is returned.
    No patient fields, raw detector values, filenames or paths are returned.
    """
    if isinstance(source, (str, Path)) and Path(source).stat().st_size > 500_000_000:
        raise ValueError("XCC file exceeds configured 500 MB safety limit")
    root = SafeET.parse(source).getroot()
    if root.tag != "PTW":
        raise ValueError("Expected PTW root element in XCC")
    content = root.find("Content")
    if content is None:
        raise ValueError("Missing PTW Content node")
    meas = content.find("MeasData")
    if meas is None:
        raise ValueError("Missing XCC MeasData")
    frames = meas.findall("Measurement")
    if not frames:
        raise ValueError("XCC has no measurement frames")
    if len(frames) > 20000:
        raise ValueError("XCC contains too many measurement frames")
    times: list[float] = []
    angles: list[float] = []
    detectors: int | None = None
    for frame in frames:
        t = _optional_float(frame.findtext("Time"))
        angle = _optional_float(frame.findtext("Angle"))
        if t is None or angle is None:
            raise ValueError("Measurement frame missing Time/Angle")
        encoded = (frame.findtext("Data") or "").strip()
        if len(encoded) > 10_000_000:
            raise ValueError("Oversized XCC detector payload")
        try:
            blob = base64.b64decode(encoded, validate=True)
        except (ValueError, base64.binascii.Error) as exc:
            raise ValueError("Invalid XCC base64 detector payload") from exc
        if not blob or len(blob) % 4:
            raise ValueError("XCC detector payload is not float32-compatible")
        n = len(blob) // 4
        if detectors is None:
            detectors = n
        elif n != detectors:
            raise ValueError("Inconsistent XCC detector counts")
        times.append(t)
        angles.append(angle)
    if any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError("XCC time stamps must increase strictly")
    intervals = [b-a for a,b in zip(times,times[1:])]
    preset = _optional_float(content.findtext("MeasurementPreset/MeasTime"))
    return XCCMetadata(
        frame_count=len(frames),
        detectors_per_frame=detectors or 0,
        start_time_s=times[0],
        last_time_s=times[-1],
        median_interval_s=float(statistics.median(intervals)) if intervals else None,
        first_angle_deg=angles[0],
        last_angle_deg=angles[-1],
        minimum_angle_deg=min(angles),
        maximum_angle_deg=max(angles),
        measurement_preset_s=preset,
        device_store_mode=content.findtext("DetectorArray/DeviceStoreMode"),
        matrix_resolution_lr_mm=_optional_float(content.findtext("DetectorArray/MatrixResolutionLR")),
        matrix_resolution_gt_mm=_optional_float(content.findtext("DetectorArray/MatrixResolutionGT")),
    )
