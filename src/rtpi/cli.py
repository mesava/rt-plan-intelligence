"""Anonymous batch CSV export from RTPLAN inputs."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
from .dicom import read_rtplan

def main() -> None:
    p=argparse.ArgumentParser(description="Research RTPLAN metrics, not clinical QA approval")
    p.add_argument("plans",nargs="+",type=Path)
    p.add_argument("--output",type=Path,default=Path("rtpi_metrics.csv"))
    args=p.parse_args()
    rows=[]
    for idx,path in enumerate(args.plans,1):
        try:
            for beam in read_rtplan(path): rows.append({"plan_index":idx,**beam})
        except (ValueError,OSError,KeyError,AttributeError) as exc:
            p.error(f"Plan #{idx} could not be analyzed: {exc}")
    if not rows: p.error("No supported beams found")
    with args.output.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Exported {len(rows)} anonymous beam summaries to {args.output}")

if __name__=="__main__":
    main()
