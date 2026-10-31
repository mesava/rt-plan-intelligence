"""Conservative DICOM RTPLAN MLCX extraction, with no patient identifiers."""
from __future__ import annotations
from io import BytesIO
from pathlib import Path
from typing import BinaryIO
import numpy as np
from .metrics import Aperture, calculate_beam, control_point_records

def _positions(sequence,previous:dict) -> dict:
    result=previous.copy()
    for item in sequence:
        result[str(item.RTBeamLimitingDeviceType)]=np.asarray(item.LeafJawPositions,dtype=float)
    return result

def _jaw(positions:dict, keys:tuple[str,...], default:tuple[float,float]) -> tuple[float,float]:
    for key in keys:
        if key in positions:
            jaw=positions[key]
            if len(jaw)!=2 or not np.all(np.isfinite(jaw)):
                raise ValueError(f"Invalid {key} jaw positions")
            return float(jaw[0]),float(jaw[1])
    return default

def read_rtplan(source:str|Path|BinaryIO|BytesIO, *, include_control_points:bool=False) -> list[dict]:
    """Read supported MLCX beams; skip non-MLCX; reject malformed structures."""
    import pydicom
    ds=pydicom.dcmread(source,stop_before_pixels=True)
    if str(getattr(ds,"Modality",""))!="RTPLAN":
        raise ValueError("Expected a DICOM RTPLAN file")
    if not getattr(ds,"BeamSequence",None):
        raise ValueError("RTPLAN has no BeamSequence")
    mu_by_beam={}
    for fg in getattr(ds,"FractionGroupSequence",[]):
        for ref in getattr(fg,"ReferencedBeamSequence",[]):
            number=int(ref.ReferencedBeamNumber)
            if not hasattr(ref,"BeamMeterset"): continue
            mu=float(ref.BeamMeterset)
            if number in mu_by_beam and not np.isclose(mu_by_beam[number],mu):
                raise ValueError("Conflicting BeamMeterset between fraction groups")
            mu_by_beam[number]=mu
    summaries=[]
    for beam in ds.BeamSequence:
        devices={str(d.RTBeamLimitingDeviceType):d for d in getattr(beam,"BeamLimitingDeviceSequence",[])}
        if "MLCX" not in devices: continue
        boundaries=np.asarray(devices["MLCX"].LeafPositionBoundaries,dtype=float)
        n=len(boundaries)-1
        if n<1 or not np.all(np.isfinite(boundaries)) or not np.all(np.diff(boundaries)>0):
            raise ValueError("Invalid MLCX leaf boundaries")
        cps=getattr(beam,"ControlPointSequence",None)
        if not cps: raise ValueError("MLCX beam has no control points")
        prior={}; aps=[]; cmw=[]; angles=[]; angle=None
        for cp in cps:
            prior=_positions(getattr(cp,"BeamLimitingDevicePositionSequence",[]),prior)
            if "MLCX" not in prior:
                raise ValueError("Missing MLCX positions at initial control point")
            leaf=prior["MLCX"]
            if len(leaf)!=2*n:
                raise ValueError("MLCX leaf position count mismatch")
            xjaw=_jaw(prior,("ASYMX","X"),(-np.inf,np.inf))
            yjaw=_jaw(prior,("ASYMY","Y"),(float(boundaries[0]),float(boundaries[-1])))
            aps.append(Aperture(leaf[:n].copy(),leaf[n:].copy(),boundaries,xjaw,yjaw))
            if not hasattr(cp,"CumulativeMetersetWeight"):
                raise ValueError("Control point missing CumulativeMetersetWeight")
            cmw.append(float(cp.CumulativeMetersetWeight))
            if hasattr(cp,"GantryAngle"): angle=float(cp.GantryAngle)
            angles.append(angle)
        if hasattr(beam,"FinalCumulativeMetersetWeight") and not np.isclose(float(beam.FinalCumulativeMetersetWeight),cmw[-1],atol=1e-5):
            raise ValueError("Final cumulative meterset weight does not match last CP")
        num=int(beam.BeamNumber)
        out=calculate_beam(num,f"Beam {num}",aps,cmw,mu_by_beam.get(num)).to_dict()
        # The VMAT interval MCSv is not a score for static-gantry IMRT.
        # Infer arc rotation only from actual changes in gantry angle at CPs.
        rotates = any(
            a is not None and b is not None and
            abs(((b-a+180.) % 360.) - 180.) > 1e-5
            for a,b in zip(angles[:-1],angles[1:])
        )
        if not rotates:
            out["mcs_vmat_jaw_clipped"] = None
            out["notes"] = (
                "No changing gantry angle: VMAT MCSv not applicable. "
                "LSV/beamAAV are geometric CP descriptors only."
            )

        if include_control_points: out["control_point_data"]=control_point_records(aps,cmw,angles)
        summaries.append(out)
    if not summaries:
        raise ValueError("No supported MLCX beams found (MLCY/non-MLC unsupported in v0.1)")
    return summaries
