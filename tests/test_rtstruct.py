"""Only synthetic RTSTRUCTs and plan metadata. Never commit patient files."""
from io import BytesIO
from zipfile import ZipFile
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.filewriter import dcmwrite
from rtpi.rtstruct import (
    read_rtstruct_archive, infer_site_candidate, link_rtplan_rtstruct, link_rtdose_rtstruct_frame
)


def _rtstruct(roi_names,patient="SYNTHETIC_INTERNAL"):
    d=Dataset()
    d.SOPClassUID="1.2.840.10008.5.1.4.1.1.481.3"
    d.SOPInstanceUID="1.2.826.0.1.3680043.2.1125.44"
    d.Modality="RTSTRUCT"
    d.PatientID=patient
    d.StudyInstanceUID="1.2.826.0.1.3680043.2.1125.1"
    d.FrameOfReferenceUID="1.2.826.0.1.3680043.2.1125.2"
    rois=[]
    for i,name in enumerate(roi_names,1):
        r=Dataset();r.ROINumber=i;r.ROIName=name;rois.append(r)
    d.StructureSetROISequence=Sequence(rois)
    return d


def _plan(label="Esoph"):
    d=Dataset()
    d.PatientID="SYNTHETIC_INTERNAL"
    d.StudyInstanceUID="1.2.826.0.1.3680043.2.1125.1"
    d.FrameOfReferenceUID="1.2.826.0.1.3680043.2.1125.999"
    d.RTPlanLabel=label
    s=Dataset();s.ReferencedSOPInstanceUID="1.2.826.0.1.3680043.2.1125.333"
    d.ReferencedStructureSetSequence=Sequence([s])
    return d


def test_zip_matches_patient_and_study_but_rejects_unsafe_geometry(tmp_path):
    rs=_rtstruct(["PTV_Esoph", "PTV_Larynx", "SpinalCord"])
    buffer=BytesIO()
    dcmwrite(buffer,rs,implicit_vr=True,little_endian=True,enforce_file_format=False)
    archive=tmp_path/"synthetic.zip"
    with ZipFile(archive,"w") as z:z.writestr("synthetic.dcm",buffer.getvalue())
    index=read_rtstruct_archive(archive)
    assert len(index)==1
    a=link_rtplan_rtstruct(_plan("PTVEsoph"),index)
    b=link_rtplan_rtstruct(_plan("PTVLarynx"),index)
    assert a["rtstruct_match"]=="patient_study_match"
    assert a["rtstruct_sop_reference"]=="mismatch"
    assert a["rtstruct_frame_of_reference"]=="mismatch"
    assert a["geometry_safe_to_combine"] is False
    assert a["site_candidate"]=="thorax_esophagus"
    assert b["site_candidate"]=="head_neck"
    assert a["site_verified"] is False and b["site_verified"] is False
    assert "SYNTHETIC_INTERNAL" not in str(a)
    assert "SYNTHETIC_INTERNAL" not in str(b)


def test_ambiguous_rtstruct_patient_study_fails_closed(tmp_path):
    rs=_rtstruct(["PTV_Esoph"])
    buffer=BytesIO()
    dcmwrite(buffer,rs,implicit_vr=True,little_endian=True,enforce_file_format=False)
    zpath=tmp_path/"two.zip"
    with ZipFile(zpath,"w") as z:
        z.writestr("a.dcm",buffer.getvalue())
        z.writestr("b.dcm",buffer.getvalue())
    index=read_rtstruct_archive(zpath)
    result=link_rtplan_rtstruct(_plan(),index)
    assert result["rtstruct_match"]=="ambiguous"
    assert result["site_candidate"]=="unknown"
    assert result["geometry_safe_to_combine"] is False


def test_explicit_target_plan_precedence_for_shared_struct():
    mixed=("PTV_Esoph","PTV_Larynx","Esophagus","Larynx")
    assert infer_site_candidate("PTVEsoph50Gy",mixed)[0]=="thorax_esophagus"
    assert infer_site_candidate("PTVLarynx60Gy",mixed)[0]=="head_neck"
    assert infer_site_candidate("NoTarget",("CTV","Lung_L","Heart"))[0]=="unknown"


def test_zip_blocks_non_dicom_and_zero_files(tmp_path):
    zpath=tmp_path/"bad.zip"
    with ZipFile(zpath,"w") as z:z.writestr("notes.txt","not a DICOM")
    try:read_rtstruct_archive(zpath)
    except ValueError as exc: assert "DICOM" in str(exc)
    else:raise AssertionError("Accepted non-DICOM")


def test_dose_structure_frame_check_is_independent(tmp_path):
    rs=_rtstruct(["PTV_Esoph"])
    bio=BytesIO()
    dcmwrite(bio,rs,implicit_vr=True,little_endian=True,enforce_file_format=False)
    zp=tmp_path/"rs.zip"
    with ZipFile(zp,"w") as z:z.writestr("synthetic.dcm",bio.getvalue())
    structures=read_rtstruct_archive(zp)
    dose=Dataset()
    dose.PatientID=rs.PatientID
    dose.StudyInstanceUID=rs.StudyInstanceUID
    dose.FrameOfReferenceUID="1.2.826.0.1.3680043.2.1125.999"
    assert link_rtdose_rtstruct_frame(dose,structures)["rtdose_rtstruct_frame_of_reference"]=="mismatch"
    dose.FrameOfReferenceUID=rs.FrameOfReferenceUID
    assert link_rtdose_rtstruct_frame(dose,structures)["rtdose_rtstruct_frame_of_reference"]=="match"
    dose.PatientID="OTHER"
    assert link_rtdose_rtstruct_frame(dose,structures)["rtdose_rtstruct_frame_of_reference"]=="missing"
