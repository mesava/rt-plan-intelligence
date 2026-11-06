"""Local-only anonymized PTW VeriSoft QA batch import. No upload."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
from .qa import parse_verisoft_pdf_text, parse_xcc_metadata

FIELDS = ("case_index","bundle_status","plan_present","dose_present","xcc_present",
          "pdf_present","beam_count","control_points","monitor_units",
          "gamma_pass_percent","gamma_fail_percent","evaluated_points","passed_points",
          "failed_points","dose_difference_percent","dta_mm","low_dose_threshold_percent",
          "low_dose_reference","gamma_passing_limit","xcc_detector_count",
          "xcc_measurement_count","xcc_interval_s","xcc_last_time_s",
          "xcc_storage_mode","quality_flags")

def _pdf_text(path: Path) -> str:
    try:
        import fitz
    except ImportError as exc:
        raise RuntimeError("Local QA option missing: pip install -e .[qa]") from exc
    with fitz.open(path) as doc:
        return "\n".join(page.get_text(sort=True) for page in doc)

def summarize_case(folder: Path, *, with_plan: bool=True) -> dict:
    """Only pseudonymous numeric case indices leave the parser."""
    if not folder.is_dir() or not folder.name.isdecimal():
        raise ValueError("Case directory must be numeric")
    files=[f for f in folder.iterdir() if f.is_file()]
    kinds={
        "plan":[f for f in files if f.suffix.lower()==".dcm" and not f.name.lower().endswith("_dose.dcm")],
        "dose":[f for f in files if f.name.lower().endswith("_dose.dcm")],
        "xcc":[f for f in files if f.suffix.lower()==".xcc"],
        "pdf":[f for f in files if f.suffix.lower()==".pdf"],
    }
    flags=[f"{k}_count_{len(v)}" for k,v in kinds.items() if len(v)!=1]
    row={k:None for k in FIELDS}
    row["case_index"]=int(folder.name)
    row["bundle_status"]="complete" if all(len(v)==1 for v in kinds.values()) else "requires_review"
    for kind in kinds: row[kind+"_present"]=int(len(kinds[kind])==1)
    if len(kinds["pdf"])==1:
        try:
            result=parse_verisoft_pdf_text(_pdf_text(kinds["pdf"][0]))
            for k in ("gamma_pass_percent","gamma_fail_percent","evaluated_points",
                      "passed_points","failed_points","dose_difference_percent","dta_mm",
                      "low_dose_threshold_percent","low_dose_reference","gamma_passing_limit"):
                row[k]=getattr(result,k)
            flags.extend(result.warnings)
        except Exception as e:
            flags.append("pdf_error_"+type(e).__name__)
    if len(kinds["xcc"])==1:
        try:
            x=parse_xcc_metadata(kinds["xcc"][0])
            row["xcc_detector_count"]=x.detector_count
            row["xcc_measurement_count"]=x.measurement_count
            row["xcc_interval_s"]=x.sample_interval_s
            row["xcc_last_time_s"]=x.last_sample_time_s
            row["xcc_storage_mode"]="DoseAccumulated" if x.storage_mode=="DoseAccumulated" else "other"
            flags.extend(x.warnings)
        except Exception as e:
            flags.append("xcc_error_"+type(e).__name__)
    if with_plan and len(kinds["plan"])==1:
        try:
            from .dicom import read_rtplan
            beams=read_rtplan(kinds["plan"][0])
            row["beam_count"]=len(beams)
            row["control_points"]=sum(b["control_points"] for b in beams)
            if all(b["monitor_units"] is not None for b in beams):
                row["monitor_units"]=sum(b["monitor_units"] for b in beams)
        except Exception as e:
            flags.append("rtplan_error_"+type(e).__name__)
    row["quality_flags"]=";".join(sorted(set(flags)))
    return row

def main() -> None:
    parser=argparse.ArgumentParser(description="Private local research import; not clinical acceptance")
    parser.add_argument("case_root",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--no-plan-metrics",action="store_true")
    a=parser.parse_args()
    if not a.case_root.is_dir():parser.error("Case root must be a directory")
    cases=sorted((p for p in a.case_root.iterdir() if p.is_dir() and p.name.isdecimal()),key=lambda p:int(p.name))
    if not cases:parser.error("No numbered case folders found")
    if a.output.resolve().is_relative_to(a.case_root.resolve()):
        parser.error("Write the confidential CSV outside the raw data tree")
    rows=[summarize_case(p,with_plan=not a.no_plan_metrics) for p in cases]
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"Cases: {len(rows)}; PDF gamma: {sum(r['gamma_pass_percent'] is not None for r in rows)}; "
          f"quality flags: {sum(bool(r['quality_flags']) for r in rows)}")
    print("No original names, identifiers or source paths in the export.")
if __name__=="__main__":
    main()
