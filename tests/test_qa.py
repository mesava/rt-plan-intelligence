"""Only synthetic records; no clinical identifiers or reports committed."""
import base64
from rtpi.qa import parse_verisoft_pdf_text,parse_xcc_metadata

REPORT_A="""Gamma 3D - Parameters
2,0 mm Distance- To- Agreement
3,0 % Dose difference with ref. to max. dose of calculated volume
Suppress dose below 10,0 % of max. dose of calculated volume
Evaluated Dose Points 7 006 ( 64,8 %)
Passed 6 783 ( 96,8 %)
Failed 223 ( 3,2 %)
Passing criteria Gamma ≤ 1,0
"""
REPORT_B="""Volume Analysis - Parameters
2,0 mm Distance- To- Agreement
3,0 % Dose difference with ref. to max. dose of calculated volume
Suppress dose below 10,0 % of the normalization value (= 1.854 Gy) of Data Set A
Evaluated Voxels 591 421 ( 53,1 %)
Passed 538 682 ( 91,1 %)
Failed 52 739 ( 8,9 %)
Passing criteria Gamma ≤ 1,0
"""
def test_parses_dose_points():
    r=parse_verisoft_pdf_text(REPORT_A)
    assert r.gamma_pass_percent==96.8
    assert r.evaluated_points==7006
    assert r.passed_points==6783 and r.failed_points==223
    assert (r.dose_difference_percent,r.dta_mm,r.low_dose_threshold_percent)==(3.,2.,10.)
    assert r.low_dose_reference=="calculated_volume_max"
    assert not r.warnings

def test_parses_volume_reference_separately():
    r=parse_verisoft_pdf_text(REPORT_B)
    assert r.gamma_pass_percent==91.1
    assert r.evaluated_points==591421
    assert r.low_dose_reference=="normalization_value"
    assert r.gamma_passing_limit==1.
    assert not r.warnings

def test_does_not_guess_empty_or_conflicting_results():
    assert "pdf_text_empty_needs_manual_review_or_ocr" in parse_verisoft_pdf_text("").warnings
    r=parse_verisoft_pdf_text(REPORT_A+"Passed 999 ( 99,9 %)\n")
    assert r.gamma_pass_percent is None
    assert "multiple_passed_results" in r.warnings

def test_detects_bad_count_totals():
    r=parse_verisoft_pdf_text(REPORT_A.replace("Failed 223","Failed 224"))
    assert "point_counts_inconsistent" in r.warnings

def test_synthetic_xcc_metadata(tmp_path):
    numbers=base64.b64encode(bytes([1,0,2,0,3,0])).decode()
    data=base64.b64encode(bytes(12)).decode()
    xml=f"""<?xml version="1.0"?><PTW><Version>1.1.0.0</Version><Content>
    <Adminsitrative><Module>VeriSoft</Module></Adminsitrative>
    <Detector>PTW_DETECTOR_OCTAVIUS_1500</Detector>
    <DeviceStoreMode>DoseAccumulated</DeviceStoreMode><IntervalTime>0.2</IntervalTime>
    <DetectorNumbers>{numbers}</DetectorNumbers><MeasData>
    <Measurement><Time>0.2</Time><Data>{data}</Data></Measurement>
    <Measurement><Time>0.4</Time><Data>{data}</Data></Measurement>
    </MeasData></Content></PTW>"""
    f=tmp_path/"synthetic.xcc";f.write_text(xml)
    x=parse_xcc_metadata(f)
    assert x.measurement_count==2
    assert x.detector_count==3 and x.values_per_measurement==3
    assert x.last_sample_time_s==0.4 and x.sample_interval_s==0.2
    assert x.storage_mode=="DoseAccumulated"
    assert not x.warnings

def test_rejects_xml_doctype(tmp_path):
    f=tmp_path/"bad.xcc";f.write_text("<!DOCTYPE A><PTW/>")
    try:parse_xcc_metadata(f)
    except ValueError as e: assert "DOCTYPE" in str(e)
    else:assert False

def test_local_case_summary_never_exports_names(tmp_path,monkeypatch):
    from rtpi.qa_cli import summarize_case
    d=tmp_path/"7"
    d.mkdir()
    (d/"Patient Family.pdf").write_bytes(b"synthetic")
    (d/"Measurement Patient.xcc").write_bytes(b"synthetic")
    (d/"123456_QA1.dcm").write_bytes(b"synthetic")
    (d/"123456_QA1_Dose.dcm").write_bytes(b"synthetic")
    monkeypatch.setattr("rtpi.qa_cli._pdf_text",lambda _:REPORT_A)
    class X:
        detectors_per_frame=1405
        frame_count=472
        median_interval_s=.2
        last_time_s=94.4
        device_store_mode="DoseAccumulated"
    monkeypatch.setattr("rtpi.qa_cli.read_xcc_metadata",lambda _:X())
    r=summarize_case(d,with_plan=False)
    assert r["case_index"]==7 and r["bundle_status"]=="complete"
    assert r["gamma_pass_percent"]==96.8
    assert r["xcc_detector_count"]==1405
    assert r["qa_analysis_kind"]=="Gamma 3D"
    assert "Patient" not in str(r) and "123456" not in str(r)

def test_ambiguous_bundle_requires_review(tmp_path):
    from rtpi.qa_cli import summarize_case
    d=tmp_path/"10"
    d.mkdir()
    (d/"one.pdf").write_bytes(b"")
    (d/"two.pdf").write_bytes(b"")
    r=summarize_case(d,with_plan=False)
    assert r["bundle_status"]=="requires_review"
    assert "pdf_count_2" in r["quality_flags"]
