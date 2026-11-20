"""Private local RTSTRUCT index and cautious plan-to-structure matching.

Never publish input RTPLAN/RTSTRUCT. Uses identifiers only in process memory.
Matching patient/study does not demonstrate geometric registration.
"""
from __future__ import annotations
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
import re
from zipfile import ZipFile, BadZipFile


@dataclass(frozen=True)
class StructureIndex:
    # Sensitive raw metadata: NEVER serialize this object to public outputs.
    patient_id: str
    study_uid: str
    sop_uid: str
    frame_uid: str
    roi_names: tuple[str, ...]


def read_rtstruct_archive(source: str | Path) -> list[StructureIndex]:
    """Read a ZIP in memory, with decompression limits and SOP checks."""
    import pydicom
    from pydicom.uid import RTStructureSetStorage
    from pydicom.errors import InvalidDicomError
    found = []
    with ZipFile(source) as z:
        files = [item for item in z.infolist() if not item.is_dir()]
        if len(files) > 128 or sum(f.file_size for f in files) > 550_000_000:
            raise ValueError("RTSTRUCT archive exceeds size or file count limit")
        for info in files:
            if info.file_size > 60_000_000 or info.compress_size > 60_000_000:
                raise ValueError("Individual RTSTRUCT exceeds permitted size")
            if not info.filename.lower().endswith(".dcm"):
                raise ValueError("Only DICOM files are allowed in RTSTRUCT archive")
            b = z.read(info)
            if len(b) != info.file_size:
                raise ValueError("Incomplete RTSTRUCT archive entry")
            bio = BytesIO(b)
            try:
                ds = pydicom.dcmread(bio, stop_before_pixels=True)
            except InvalidDicomError:
                bio.seek(0)
                ds = pydicom.dcmread(bio, force=True, stop_before_pixels=True)
            if str(getattr(ds, "Modality", "")) != "RTSTRUCT" or (
                str(getattr(ds, "SOPClassUID", "")) != str(RTStructureSetStorage)
            ):
                raise ValueError("Archive contains an invalid RTSTRUCT SOP Class")
            roi_names = tuple(
                str(getattr(item, "ROIName", ""))[:256]
                for item in getattr(ds, "StructureSetROISequence", [])
            )
            found.append(StructureIndex(
                patient_id=str(getattr(ds, "PatientID", "")).strip(),
                study_uid=str(getattr(ds, "StudyInstanceUID", "")).strip(),
                sop_uid=str(getattr(ds, "SOPInstanceUID", "")).strip(),
                frame_uid=str(getattr(ds, "FrameOfReferenceUID", "")).strip(),
                roi_names=roi_names,
            ))
    if not found:
        raise ValueError("No RTSTRUCT found")
    return found


def infer_site_candidate(plan_label: str, roi_names: tuple[str, ...]) -> tuple[str, str, str]:
    """Nonvalidated plan-specific anatomical *hint*, never verified diagnosis."""
    s = " ".join(roi_names).lower()
    label = plan_label.lower()
    def has(*patterns: str) -> bool:
        return any(re.search(p, s, flags=re.IGNORECASE) for p in patterns)
    if "esoph" in label and has(r"ptv.?esoph", r"gtvp.?esoph", r"ctvp.?esoph"):
        return ("thorax_esophagus", "explicit_plan_and_target_roi", "high")
    if "larynx" in label and has(r"ptv.?larynx", r"gtvp.?larynx", r"ctvp.?larynx"):
        return ("head_neck", "explicit_plan_and_target_roi", "high")
    if has(r"ctv.?mesorectum"):
        return ("pelvis_rectum", "mesorectum_target_roi", "high")
    if has("prostate") and has(r"ptv"):
        return ("pelvis_prostate", "prostate_roi_with_ptv", "high")
    if has("breast", "mammary", r"ptv.?mam", r"ptv.?breast") and has(r"ptv", r"ctv"):
        return ("breast", "breast_roi_with_targets", "high")
    if has(r"ptv.?esoph", r"gtvp.?esoph"):
        return ("thorax_esophagus", "esophageal_target_roi", "high")
    if has(r"ptv.?larynx"):
        return ("head_neck", "larynx_target_roi", "high")
    if has("parotid", r"bone_mandib", "mandible", r"ptv.?70", r"ptv.?63") and has("ptv", "gtv"):
        return ("head_neck", "head_neck_roi_pattern", "moderate")
    if has("brain") and has("optic", "chiasm", "lens") and has("ptv"):
        return ("brain", "cranial_oar_pattern", "moderate")
    if has("cochlea", "eye", "lens") and has("ptv"):
        return ("cranial_or_orbit", "cranial_roi_pattern", "low")
    return ("unknown", "insufficient_roi_evidence", "none")


def link_rtplan_rtstruct(plan, structures: list[StructureIndex]) -> dict:
    """Return only non-identifying status and site candidate, no raw metadata."""
    pid = str(getattr(plan, "PatientID", "")).strip()
    study = str(getattr(plan, "StudyInstanceUID", "")).strip()
    reference_uid = {
        str(getattr(r, "ReferencedSOPInstanceUID", "")).strip()
        for r in getattr(plan, "ReferencedStructureSetSequence", [])
    }
    matches = [s for s in structures if pid and study and
               s.patient_id == pid and s.study_uid == study]
    if len(matches) != 1:
        return {
            "rtstruct_match": "missing" if not matches else "ambiguous",
            "rtstruct_sop_reference": "unverified",
            "rtstruct_frame_of_reference": "unverified",
            "rtstruct_roi_count": None,
            "site_candidate": "unknown",
            "site_evidence": "no_unique_matching_structure",
            "site_candidate_confidence": "none",
            "site_verified": False,
            "geometry_safe_to_combine": False,
        }
    s = matches[0]
    sop_status = "match" if s.sop_uid and s.sop_uid in reference_uid else (
        "unverified" if not reference_uid or not s.sop_uid else "mismatch"
    )
    plan_frame = str(getattr(plan, "FrameOfReferenceUID", "")).strip()
    frame_status = "match" if s.frame_uid and s.frame_uid == plan_frame else (
        "unverified" if not s.frame_uid or not plan_frame else "mismatch"
    )
    candidate, evidence, confidence = infer_site_candidate(
        str(getattr(plan, "RTPlanLabel", "")), s.roi_names
    )
    # A same-frame/same-SOP check is necessary but not sufficient for spatial QA.
    return {
        "rtstruct_match": "patient_study_match",
        "rtstruct_sop_reference": sop_status,
        "rtstruct_frame_of_reference": frame_status,
        "rtstruct_roi_count": len(s.roi_names),
        "site_candidate": candidate,
        "site_evidence": evidence,
        "site_candidate_confidence": confidence,
        "site_verified": False,
        "geometry_safe_to_combine": False,
    }
