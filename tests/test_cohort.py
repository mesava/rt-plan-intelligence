"""Synthetic stratification tests. No patient identifiers, dose or QA files."""
from pathlib import Path
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from rtpi.cohort import (
    beam_stratum, summarize_rtplan, dose_link_status,
    classify_private_folders, read_site_mapping,
)

def plan(patient_id="TEST",angle2=120,beam_type="DYNAMIC",energy=6.,mode="STANDARD"):
    d=Dataset()
    d.Modality="RTPLAN"
    d.SOPInstanceUID="1.2.826.0.1.3680043.2.1125.1"
    d.PatientID=patient_id
    b=Dataset()
    b.BeamNumber=1
    b.BeamType=beam_type
    b.TreatmentMachineName="INTERNAL_ROOM"
    m=Dataset();m.FluenceMode=mode
    b.PrimaryFluenceModeSequence=Sequence([m])
    cp0=Dataset();cp0.GantryAngle=100;cp0.NominalBeamEnergy=energy
    cp1=Dataset();cp1.GantryAngle=angle2
    b.ControlPointSequence=Sequence([cp0,cp1])
    d.BeamSequence=Sequence([b])
    return d

def test_distinguish_rotating_and_fixed_dynamic_gantry():
    arc=plan()
    result=summarize_rtplan(arc)
    assert result["delivery_group"]=="ROTATING_ARC"
    assert result["nominal_energy_values_mv"]=="6.0"
    assert result["fluence_modes"]=="STANDARD"
    assert result["cp_count"]==2
    assert "TEST" not in str(result)
    sliding=plan(angle2=100,beam_type="DYNAMIC")
    assert beam_stratum(sliding.BeamSequence[0])["technique"]=="DYNAMIC_MLC_FIXED_GANTRY"
    segmented=plan(angle2=100,beam_type="STATIC")
    assert beam_stratum(segmented.BeamSequence[0])["technique"]=="STATIC_MLC_FIXED_GANTRY"

def test_missing_angle_is_not_auto_static():
    d=plan()
    for cp in d.BeamSequence[0].ControlPointSequence:
        del cp.GantryAngle
    assert beam_stratum(d.BeamSequence[0])["technique"]=="UNDETERMINED"

def test_normal_energy_10_and_fluence_nonstandard():
    d=plan(energy=10.,mode="NON_STANDARD")
    b=beam_stratum(d.BeamSequence[0])
    assert b["nominal_energy_mv"]==10.
    assert b["fluence_mode"]=="NON_STANDARD"

def test_dose_matching_uid_and_patient():
    p=plan(patient_id="NEVER_EXPORT")
    d=Dataset()
    d.Modality="RTDOSE"
    d.PatientID="NEVER_EXPORT"
    ref=Dataset()
    ref.ReferencedSOPInstanceUID=p.SOPInstanceUID
    d.ReferencedRTPlanSequence=Sequence([ref])
    result=dose_link_status(p,d)
    assert result=={"dose_rtplan_uid_link":"match","dose_patient_id_link":"match"}
    assert "NEVER_EXPORT" not in str(result)
    ref.ReferencedSOPInstanceUID="1.2.3.4"
    assert dose_link_status(p,d)["dose_rtplan_uid_link"]=="mismatch"

def test_private_patient_groups_do_not_leak_ids(tmp_path,monkeypatch):
    p=plan(patient_id="IDENTIFIER_THAT_MUST_NEVER_LEAVE_MEMORY")
    dose=Dataset();dose.Modality="RTDOSE"
    dose.PatientID=p.PatientID
    r=Dataset();r.ReferencedSOPInstanceUID=p.SOPInstanceUID
    dose.ReferencedRTPlanSequence=Sequence([r])
    for index in (1,2):
        c=tmp_path/str(index);c.mkdir()
        (c/"synthetic.dcm").write_bytes(b"x")
        (c/"synthetic_Dose.dcm").write_bytes(b"x")
    monkeypatch.setattr("rtpi.cohort.read_dicom",
                        lambda path,modality:p if modality=="RTPLAN" else dose)
    rows=classify_private_folders(tmp_path)
    assert rows[0]["patient_group"]==rows[1]["patient_group"]
    assert all(row["repeated_patient_group"] for row in rows)
    assert all(row["site_label"]=="unknown" for row in rows)
    assert all(row["dose_rtplan_uid_link"]=="match" for row in rows)
    assert "IDENTIFIER_THAT_MUST_NEVER_LEAVE_MEMORY" not in str(rows)
    assert "INTERNAL_ROOM" not in str(rows)

def test_site_mapping_is_validated_nonidentifying(tmp_path):
    p=tmp_path/"site.csv"
    p.write_text("case_index,site_label\n1,head_neck\n2,unknown\n",encoding="utf8")
    assert read_site_mapping(p)=={1:"head_neck",2:"unknown"}
    p.write_text("case_index,site_label\n1,Patient ABC\n",encoding="utf8")
    try: read_site_mapping(p)
    except ValueError: pass
    else: raise AssertionError("Potential patient label must be rejected")
