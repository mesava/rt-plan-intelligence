"""Synthetic-only tests, with no clinical XCC, PDF or patient information."""
import base64
from io import BytesIO, StringIO
import struct
import pytest

from rtpi.xcc import read_xcc_metadata
from rtpi.qa import parse_qa_text, cohort_summary_csv


def xcc_fixture(times=(0.2, 0.4, 0.6)):
    frames = []
    for i, t in enumerate(times):
        data = base64.b64encode(struct.pack("<3f", 0.1*i, 0.2*i, 0.3*i)).decode()
        frames.append(
            f"<Measurement><Time>{t}</Time><Angle>{140-i*5}</Angle>"
            f"<Data>{data}</Data></Measurement>"
        )
    return (
        "<PTW><Content><DetectorArray><DeviceStoreMode>DoseAccumulated"
        "</DeviceStoreMode><MatrixResolutionLR>5</MatrixResolutionLR>"
        "<MatrixResolutionGT>5</MatrixResolutionGT></DetectorArray>"
        "<MeasurementPreset><MeasTime>0.6</MeasTime></MeasurementPreset>"
        "<MeasData>" + "".join(frames) + "</MeasData></Content></PTW>"
    ).encode()


def test_xcc_summary_has_no_patient_data():
    info = read_xcc_metadata(BytesIO(xcc_fixture()))
    assert info.frame_count == 3
    assert info.detectors_per_frame == 3
    assert info.median_interval_s == pytest.approx(.2)
    assert info.last_angle_deg == 130
    assert info.device_store_mode == "DoseAccumulated"
    assert "patient" not in str(info.to_dict()).lower()


def test_xcc_bad_frame_count():
    text = xcc_fixture().replace(b"<Time>0.6</Time>", b"<Time>0.3</Time>")
    with pytest.raises(ValueError, match="increase"):
        read_xcc_metadata(BytesIO(text))


def test_xcc_rejects_xml_entities():
    s = b'<!DOCTYPE PTW [<!ENTITY exploit "expanded">]><PTW>&exploit;</PTW>'
    with pytest.raises(Exception):
        read_xcc_metadata(BytesIO(s))


def test_qa_pdf_text_template_gamma3d():
    pdf_text = """
    Patient Name NotForExport
    Gamma 3D - Parameters
    2,0 mm Distance- To- Agreement
    3,0 % Dose difference with ref. to max. dose of calculated volume
    Suppress dose below 10,0 % of max. dose
    Passed 5 386 ( 94,5 %)
    Failed 315 ( 5,5 %)
    Result 94,5 % (Yellow)
    """
    qa = parse_qa_text(pdf_text)
    assert qa.gamma_pass_percent == 94.5
    assert qa.dose_difference_percent == 3
    assert qa.dta_mm == 2
    assert qa.low_dose_threshold_percent == 10
    assert qa.passed_count == 5386
    assert qa.failed_count == 315
    assert "NotForExport" not in str(qa.to_dict())


def test_volume_qa_template():
    qa = parse_qa_text(
        "Volume Analysis - Parameters\n"
        "2,0 mm Distance- To- Agreement\n"
        "3,0 % Dose difference\n"
        "Result 91,1 % (Yellow)"
    )
    assert qa.gamma_pass_percent == 91.1
    assert qa.analysis_kind == "Volume Analysis"


def test_scanned_report_no_invented_gamma():
    qa = parse_qa_text("")
    assert qa.gamma_pass_percent is None
    assert qa.dta_mm is None


def test_csv_privacy_minimized_and_criteria_coverage():
    text = (
        "case,patient_id,patient_name,gamma_result_percent,"
        "dta_mm,dose_diff_percent,threshold_percent\n"
        "1,00001,NotForExport,94.5,2,3,10\n"
        "2,00002,PrivateName,98.2,,,\n"
    )
    cohort = cohort_summary_csv(StringIO(text))
    assert cohort["plans"] == 2
    assert cohort["gamma_below_95_count"] == 1
    assert cohort["explicit_dta_count"] == 1
    assert "NotForExport" not in str(cohort)
    assert "PrivateName" not in str(cohort)


def test_csv_duplicate_case_rejected():
    src = StringIO("case,gamma_result_percent\n1,97.1\n1,95.3\n")
    with pytest.raises(ValueError, match="Duplicate"):
        cohort_summary_csv(src)
