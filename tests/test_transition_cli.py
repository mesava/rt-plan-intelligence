"""Synthetic-only anonymized CSV export test."""
import csv
from pathlib import Path

from rtpi.transition_cli import export_transition_tables


def test_exports_numeric_records_without_original_filename(tmp_path,monkeypatch):
    fake={
        "beam_number":1,
        "monitor_units":100.0,
        "mcs_vmat_jaw_clipped":.7,
        "transition_summary":{
            "transition_count":1,
            "max_single_leaf_step_mm":2.0,
            "kinematic_min_total_seconds":None,
        },
        "transition_data":[{
            "cp_from":0,"cp_to":1,"delta_mu":100.0,
            "max_leaf_mm_per_mu":.02,
            "flags":"",
        }],
        "sector_data":[{
            "sector_start_deg":0,"sector_end_deg":30,"sum_delta_mu":100.0,
        }],
    }
    monkeypatch.setattr("rtpi.transition_cli.read_rtplan",lambda *a,**k:[fake])
    out=tmp_path/"out"
    counts=export_transition_tables(
        [Path("/secure/DO_NOT_EXPORT_IDENTIFIERS_PATIENTNAME.dcm")],out
    )
    assert counts=={"beams":1,"transitions":1,"sectors":1}
    records=list(csv.DictReader((out/"rtpi_cp_transitions.csv").open(encoding="utf8")))
    assert records[0]["plan_index"]=="1"
    assert records[0]["delta_mu"]=="100.0"
    for f in out.iterdir():
        s=f.read_text(encoding="utf8")
        assert "DO_NOT_EXPORT" not in s
        assert "PATIENTNAME" not in s


def test_static_beam_empty_sector_file_is_valid_csv(tmp_path,monkeypatch):
    beam={
        "beam_number":1,"monitor_units":10,"mcs_vmat_jaw_clipped":None,
        "transition_summary":{},"transition_data":[],"sector_data":[]
    }
    monkeypatch.setattr("rtpi.transition_cli.read_rtplan",lambda *a,**k:[beam])
    export_transition_tables([tmp_path/"fake.dcm"],tmp_path/"output")
    header=(tmp_path/"output"/"rtpi_vmat_sectors.csv").read_text()
    assert "plan_index,beam_number" in header
