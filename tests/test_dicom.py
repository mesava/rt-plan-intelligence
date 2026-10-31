"""Synthetic DICOM integration. No patient data."""
from io import BytesIO
import pytest
pydicom=pytest.importorskip("pydicom")
from pydicom.dataset import Dataset,FileDataset,FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import RTPlanStorage,ExplicitVRLittleEndian,generate_uid
from rtpi.dicom import read_rtplan

def synthetic_plan(modality="RTPLAN",orientation="MLCX",bad_boundaries=False,arc=False):
    meta=FileMetaDataset()
    meta.MediaStorageSOPClassUID=RTPlanStorage
    meta.MediaStorageSOPInstanceUID=generate_uid()
    meta.TransferSyntaxUID=ExplicitVRLittleEndian
    ds=FileDataset("synthetic",{},file_meta=meta,preamble=b"\0"*128)
    ds.Modality=modality;ds.SOPClassUID=RTPlanStorage;ds.SOPInstanceUID=meta.MediaStorageSOPInstanceUID
    ds.PatientName="Synthetic^Only"
    beam=Dataset();beam.BeamNumber=1;beam.BeamName="Synthetic"
    device=Dataset();device.RTBeamLimitingDeviceType=orientation
    device.NumberOfLeafJawPairs=2
    device.LeafPositionBoundaries=[-10,0,10] if not bad_boundaries else [0,0,10]
    beam.BeamLimitingDeviceSequence=Sequence([device])
    beam.FinalCumulativeMetersetWeight=1
    cp0=Dataset();cp0.ControlPointIndex=0;cp0.CumulativeMetersetWeight=0;cp0.GantryAngle=90
    pos=Dataset();pos.RTBeamLimitingDeviceType=orientation
    pos.LeafJawPositions=[-3,-5,3,5]
    cp0.BeamLimitingDevicePositionSequence=Sequence([pos])
    cp1=Dataset();cp1.ControlPointIndex=1;cp1.CumulativeMetersetWeight=1
    if arc: cp1.GantryAngle=100
    beam.ControlPointSequence=Sequence([cp0,cp1])
    ds.BeamSequence=Sequence([beam])
    ref=Dataset();ref.ReferencedBeamNumber=1;ref.BeamMeterset=100
    fg=Dataset();fg.ReferencedBeamSequence=Sequence([ref])
    ds.FractionGroupSequence=Sequence([fg])
    stream=BytesIO();ds.save_as(stream,write_like_original=False);stream.seek(0)
    return stream

def test_read_and_deidentify():
    result=read_rtplan(synthetic_plan(),include_control_points=True)[0]
    assert result["control_points"]==2
    assert result["monitor_units"]==100
    assert result["aperture_area_mean_mm2"]==160
    assert result["control_point_data"][1]["gantry_deg"]==90
    assert "Synthetic^Only" not in str(result)

def test_reject_other_modality():
    with pytest.raises(ValueError,match="RTPLAN"):
        read_rtplan(synthetic_plan(modality="RTDOSE"))

def test_reject_mlc_orientation():
    with pytest.raises(ValueError,match="MLCY"):
        read_rtplan(synthetic_plan(orientation="MLCY"))

def test_reject_bad_boundaries():
    with pytest.raises(ValueError,match="boundaries"):
        read_rtplan(synthetic_plan(bad_boundaries=True))

def test_static_gantry_does_not_get_vmat_mcs():
    result=read_rtplan(synthetic_plan())[0]
    assert result["mcs_vmat_jaw_clipped"] is None
    assert result["lsv_cp_mu_weighted"] is not None

def test_dynamic_arc_gets_vmat_mcs():
    result=read_rtplan(synthetic_plan(arc=True))[0]
    assert result["mcs_vmat_jaw_clipped"] is not None
    assert 0 <= result["mcs_vmat_jaw_clipped"] <= 1
