import numpy as np
import pytest
from rtpi.metrics import Aperture,calculate_beam,control_point_records

def ap(left,right,xjaw=(-20.,20.),yjaw=(-10.,10.)):
    return Aperture(np.asarray(left,dtype=float),np.asarray(right,dtype=float),
                    np.array([-10.,0.,10.]),xjaw,yjaw)

def test_area_and_jaw_clipping():
    assert ap([-3,-5],[3,5]).area_mm2()==160.
    assert ap([-30,-30],[30,30]).area_mm2()==800.
    assert ap([-3,-5],[3,5],yjaw=(-10,0)).area_mm2()==60.
    assert ap([-30,-30],[30,30],xjaw=(-2,2)).area_mm2()==80.

def test_travel_and_sas():
    a,b=ap([-2,-2],[2,2]),ap([-3,-3],[3,3])
    m=calculate_beam(1,"synthetic",[a,b],[0,1],100)
    assert m.leaf_travel_total_mm==4.
    assert m.leaf_travel_max_per_transition_mm==4.
    assert m.sas5_mu_weighted==0.
    assert m.sas10_mu_weighted==1.
    assert m.aperture_area_mean_mm2==b.area_mm2()

def test_closed_and_zero_weight():
    a=ap([0,0],[0,0])
    m=calculate_beam(1,"closed",[a],[1],100)
    assert m.sas5_mu_weighted is None
    assert m.aperture_area_cv is None
    assert m.leaf_gap_mean_mm is None
    b=ap([-2,-2],[2,2])
    m=calculate_beam(1,"zero",[b,b],[0,0],0)
    assert m.aperture_area_mean_mm2 is None
    assert m.sas10_mu_weighted is None

def test_invalid_meterset_and_boundaries():
    a=ap([-1,-1],[1,1])
    with pytest.raises(ValueError,match="decreasing"):
        calculate_beam(1,"bad",[a,a],[1,0],100)
    with pytest.raises(ValueError,match="align"):
        calculate_beam(1,"bad",[a],[],100)
    with pytest.raises(ValueError,match="strictly increasing"):
        Aperture(np.array([-1.,-1.]),np.array([1.,1.]),np.array([0.,0.,10.]),(-20,20),(-10,10))

def test_anonymous_cp_data():
    a=ap([-1,-1],[1,1])
    cps=control_point_records([a,a],[0,1],[0,45])
    assert cps[1]["gantry_deg"]==45
    assert cps[1]["travel_from_previous_cp_mm"]==0
    assert "PatientName" not in str(cps)
