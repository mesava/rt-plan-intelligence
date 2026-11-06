"""Private cohort stratification from RTPLAN/RTDOSE, without patient identifiers.

Research-only. Metadata labels describe DICOM, not independently verified
machine settings or clinical delivery. Original IDs never leave memory.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from collections import Counter
from typing import Any


def read_dicom(path: str | Path, expected_modality: str):
    """Accept Part-10 and validated headerless TPS exports."""
    import pydicom
    from pydicom.errors import InvalidDicomError
    from pydicom.uid import RTPlanStorage, RTDoseStorage
    try:
        ds=pydicom.dcmread(path, stop_before_pixels=True)
    except InvalidDicomError:
        ds=pydicom.dcmread(path, stop_before_pixels=True, force=True)
        expected_uid=str(RTPlanStorage if expected_modality=="RTPLAN" else RTDoseStorage)
        if str(getattr(ds,"SOPClassUID",""))!=expected_uid:
            raise ValueError("Headerless source has unexpected SOP Class")
    if str(getattr(ds,"Modality",""))!=expected_modality:
        raise ValueError("Unexpected DICOM modality")
    return ds


def beam_stratum(beam: Any) -> dict:
    """Beam-level classification, using CP inheritance and gantry geometry."""
    cps=list(getattr(beam,"ControlPointSequence",[]) or [])
    if not cps:
        raise ValueError("Beam has no control points")
    gantry=[];last_angle=None;energy=[];last_energy=None
    for cp in cps:
        if hasattr(cp,"GantryAngle"):
            last_angle=float(cp.GantryAngle)
        gantry.append(last_angle)
        if hasattr(cp,"NominalBeamEnergy"):
            last_energy=float(cp.NominalBeamEnergy)
        if last_energy is not None:
            energy.append(last_energy)
    if any(a is None for a in gantry):
        rotates=None
    else:
        rotates=any(
            abs(((b-a+180.) % 360.)-180.) > 1e-5
            for a,b in zip(gantry[:-1],gantry[1:]))
    beam_type=str(getattr(beam,"BeamType","UNKNOWN")).strip().upper()
    if rotates is None:
        technique="UNDETERMINED"
    elif rotates:
        technique="ROTATING_ARC"
    elif beam_type=="DYNAMIC":
        technique="DYNAMIC_MLC_FIXED_GANTRY"
    elif beam_type=="STATIC":
        technique="STATIC_MLC_FIXED_GANTRY"
    else:
        technique="FIXED_GANTRY_OTHER"
    fluence=[]
    for f in getattr(beam,"PrimaryFluenceModeSequence",[]) or []:
        x=str(getattr(f,"FluenceMode","")).strip()
        if x: fluence.append(x)
    fluence_modes=sorted(set(fluence))
    energies=sorted(set(energy))
    return {
        "beam_number":int(beam.BeamNumber),
        "technique":technique,
        "beam_type":beam_type,
        "nominal_energy_mv":float(energies[0]) if len(energies)==1 else None,
        "nominal_energy_values_mv":"|".join(str(x) for x in energies),
        "fluence_mode":"|".join(fluence_modes) if fluence_modes else "UNSPECIFIED",
        "control_points":len(cps),
        "gantry_rotates":rotates,
    }


def summarize_rtplan(ds: Any) -> dict:
    """Metadata only; no PlanLabel, BeamName, name, UID, IDs or free text."""
    if str(getattr(ds,"Modality",""))!="RTPLAN":
        raise ValueError("Expected RTPLAN")
    beams=[beam_stratum(b) for b in getattr(ds,"BeamSequence",[]) or []]
    if not beams:
        raise ValueError("RTPLAN has no beams")
    modes=sorted(set(b["technique"] for b in beams))
    energies=sorted(set(x for b in beams for x in b["nominal_energy_values_mv"].split("|") if x))
    fluence=sorted(set(b["fluence_mode"] for b in beams))
    kind=next(iter(modes)) if len(modes)==1 else "MIXED"
    return {
        "beam_count":len(beams),
        "cp_count":sum(b["control_points"] for b in beams),
        "rotating_beams":sum(b["gantry_rotates"] is True for b in beams),
        "delivery_group":kind,
        "delivery_modes":"|".join(modes),
        "nominal_energy_values_mv":"|".join(energies),
        "fluence_modes":"|".join(fluence),
    }


def dose_link_status(plan: Any,dose: Any) -> dict:
    """Structural DICOM linkage only; never asserts measured QA correspondence."""
    refs=getattr(dose,"ReferencedRTPlanSequence",None) or []
    plan_uid=str(getattr(plan,"SOPInstanceUID","")).strip()
    uid_match=(len(refs)==1 and bool(plan_uid) and
               str(getattr(refs[0],"ReferencedSOPInstanceUID","")).strip()==plan_uid)
    if len(refs)!=1 or not plan_uid:
        uid_status="unverified"
    else:
        uid_status="match" if uid_match else "mismatch"
    a=str(getattr(plan,"PatientID","")).strip()
    b=str(getattr(dose,"PatientID","")).strip()
    pid_status=("unverified" if not a or not b else "match" if a==b else "mismatch")
    return {"dose_rtplan_uid_link":uid_status,"dose_patient_id_link":pid_status}


def classify_private_folders(root: Path,site_labels: dict[int,str]|None=None) -> list[dict]:
    """Only numeric case indices and non-PHI enumerated groups are exported."""
    if not root.is_dir():
        raise ValueError("Cohort root is not a directory")
    files=sorted((p for p in root.iterdir() if p.is_dir() and p.name.isdecimal()),
                 key=lambda p:int(p.name))
    if not files: raise ValueError("No numeric case folders")
    cases=[]
    patient_keys={}
    machine_keys={}
    for folder in files:
        plans=[p for p in folder.iterdir() if p.is_file() and
               p.suffix.lower()==".dcm" and not p.name.lower().endswith("_dose.dcm")]
        doses=[p for p in folder.iterdir() if p.is_file() and p.name.lower().endswith("_dose.dcm")]
        if len(plans)!=1:
            raise ValueError("Expected exactly one RTPLAN per numbered case")
        plan=read_dicom(plans[0],"RTPLAN")
        row={"case_index":int(folder.name),**summarize_rtplan(plan)}
        raw_patient_id=str(getattr(plan,"PatientID","")).strip()
        if raw_patient_id:
            key=("id",raw_patient_id)
            if key not in patient_keys:
                patient_keys[key]=f"PG{len(patient_keys)+1:03d}"
            row["patient_group"]=patient_keys[key]
            row["patient_id_present"]=1
        else:
            row["patient_group"]=f"UNKNOWN_CASE_{int(folder.name)}"
            row["patient_id_present"]=0
        machines={str(getattr(b,"TreatmentMachineName","")).strip() for b in plan.BeamSequence}
        machines.discard("")
        if len(machines)==1:
            key=next(iter(machines))
            if key not in machine_keys:
                machine_keys[key]=f"M{len(machine_keys)+1:02d}"
            row["machine_group"]=machine_keys[key]
        else:
            row["machine_group"]="unknown_or_mixed"
        row["site_label"]=(site_labels or {}).get(int(folder.name),"unknown")
        row["site_source"]="user_verified" if int(folder.name) in (site_labels or {}) else "unverified"
        if len(doses)==1:
            try:
                dose=read_dicom(doses[0],"RTDOSE")
                row.update(dose_link_status(plan,dose))
            except Exception:
                row.update({"dose_rtplan_uid_link":"parse_error","dose_patient_id_link":"parse_error"})
        else:
            row.update({"dose_rtplan_uid_link":"missing_or_ambiguous",
                        "dose_patient_id_link":"missing_or_ambiguous"})
        cases.append(row)
    counts=Counter(r["patient_group"] for r in cases)
    for row in cases:
        row["repeated_patient_group"]=counts[row["patient_group"]]>1
    return cases


SITE_VALUES={"unknown","brain","head_neck","thorax","breast","abdomen","pelvis",
             "prostate","gynecologic","spine","extremity","other"}
def read_site_mapping(path: Path) -> dict[int,str]:
    """Only predetermined site categories; never copy raw patient annotations."""
    result={}
    with path.open(newline="",encoding="utf-8-sig") as f:
        for line in csv.DictReader(f):
            index=int(line["case_index"])
            label=line["site_label"].strip().lower()
            if label not in SITE_VALUES:
                raise ValueError("Site must be an allowed non-identifying category")
            if index in result: raise ValueError("Duplicate site mapping case index")
            result[index]=label
    return result


def main() -> None:
    p=argparse.ArgumentParser(description="Private RTPLAN cohort strata, no PHI output")
    p.add_argument("case_root",type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--sites",type=Path,default=None,help="Verified case_index,site_label CSV")
    args=p.parse_args()
    if args.output.resolve().is_relative_to(args.case_root.resolve()):
        p.error("Output must be outside raw case tree")
    sites=read_site_mapping(args.sites) if args.sites else None
    rows=classify_private_folders(args.case_root,sites)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Stratified {len(rows)} cases; no patient IDs, names or file names exported.")

if __name__=="__main__":
    main()
