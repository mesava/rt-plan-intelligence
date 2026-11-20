"""Pure synthetic MLCX/CP physics; not validation against measured QA."""
import numpy as np
import pytest
from rtpi.metrics import Aperture
from rtpi.transitions import MachineLimits, analyze_transitions


def ap(left, right, xjaw=(-20., 20.), yjaw=(-15., 15.)):
    return Aperture(np.array(left, dtype=float), np.array(right, dtype=float),
                    np.array([-15.,-5.,5.,15.]), xjaw, yjaw)


def test_wrap_around_and_mu_normalization_and_sectors():
    a=ap([-5]*3,[5]*3)
    b=ap([-3]*3,[7]*3)
    c=ap([-2]*3,[8]*3)
    r=analyze_transitions([a,b,c],[0,.5,1.],[359.,1.,3.],120.)
    assert r.summary["rotating_transition_count"]==2
    assert r.summary["total_delta_mu"]==pytest.approx(120.)
    assert r.summary["time_axis_available"] is False
    assert r.summary["measured_deliverability_or_gamma"] is False
    assert len(r.sectors)==1
    sector=r.sectors[0]
    assert sector["sector_start_deg"]==0
    assert sector["sum_delta_mu"]==pytest.approx(120.)
    assert sector["sum_gantry_sweep_deg"]==pytest.approx(4.)
    assert r.transitions[0]["gantry_sweep_shortest_deg"]==pytest.approx(2.)
    assert r.transitions[0]["max_single_leaf_step_mm"]==pytest.approx(2.)
    assert r.transitions[0]["max_leaf_mm_per_mu"]==pytest.approx(2/60)


def test_zero_mu_reposition_is_flagged_but_not_failed_delivery():
    a=ap([-5]*3,[5]*3)
    b=ap([-1]*3,[5]*3)
    c=ap([-1]*3,[5]*3)
    r=analyze_transitions([a,b,c],[0,0,1],[0,2,4],100.)
    assert r.transitions[0]["delta_mu"]==0
    assert "leaf_motion_without_incremental_mu" in r.transitions[0]["flags"]
    assert "gantry_motion_without_incremental_mu" in r.transitions[0]["flags"]
    assert r.transitions[0]["max_leaf_mm_per_mu"] is None
    assert r.summary["zero_mu_leaf_motion_count"]==1
    assert not any(k in r.summary for k in ("passed","gamma","delivered"))


def test_fixed_gantry_is_not_a_vmat_sector():
    a=ap([-5]*3,[5]*3)
    b=ap([-4]*3,[5]*3)
    r=analyze_transitions([a,b],[0,1],[90,90],100)
    assert r.summary["rotating_transition_count"]==0
    assert r.sectors==[]
    assert r.transitions[0]["max_leaf_mm_per_degree"] is None
    assert r.transitions[0]["max_leaf_mm_per_mu"]==pytest.approx(.01)


def test_true_mu_from_final_cmw_not_assumed_to_equal_one():
    a=ap([-5]*3,[5]*3)
    b=ap([-4]*3,[6]*3)
    r=analyze_transitions([a,b],[0,100],[0,2],200,final_cmw=100)
    assert r.transitions[0]["delta_mu"]==pytest.approx(200)
    with pytest.raises(ValueError,match="Final CMW"):
        analyze_transitions([a,b],[0,100],[0,2],200,final_cmw=1)


def test_missing_mu_is_null_not_division_by_zero():
    a=ap([-5]*3,[5]*3)
    r=analyze_transitions([a,a],[0,1],[None,None],None)
    assert r.transitions[0]["delta_mu"] is None
    assert r.transitions[0]["max_leaf_mm_per_mu"] is None
    assert "gantry_angle_unavailable" in r.transitions[0]["flags"]
    assert r.summary["total_delta_mu"] is None


def test_lower_bound_only_from_supplied_limits():
    a=ap([-5]*3,[5]*3)
    b=ap([-3]*3,[7]*3)
    bare=analyze_transitions([a,b],[0,1],[0,2],60)
    assert bare.summary["kinematic_min_total_seconds"] is None
    limited=analyze_transitions(
        [a,b],[0,1],[0,2],60,
        limits=MachineLimits(max_leaf_speed_mm_s=1,
                             max_gantry_speed_deg_s=2,
                             max_dose_rate_mu_min=600))
    # Lower bound=max(2 s leaf,1 s gantry,6 s dose-rate)
    assert limited.transitions[0]["kinematic_min_seconds_from_limits"]==pytest.approx(6.)
    assert limited.summary["kinematic_min_total_seconds"]==pytest.approx(6.)
    with pytest.raises(ValueError,match="positive"):
        MachineLimits(max_leaf_speed_mm_s=0)


def test_active_leaves_do_not_count_motion_behind_closed_aperture():
    a=ap([-5,-5,-5],[5,5,5],yjaw=(-5,5))
    b=ap([10,-4,10],[10,6,10],yjaw=(-5,5))
    r=analyze_transitions([a,b],[0,1],[0,3],100)
    assert r.transitions[0]["max_single_leaf_step_mm"]==pytest.approx(15.)
    assert r.transitions[0]["max_active_leaf_step_mm"]==pytest.approx(1.)


def test_invalid_alignment_and_sector_width():
    a=ap([-5]*3,[5]*3)
    with pytest.raises(ValueError,match="align"):
        analyze_transitions([a,a],[0,1],[0],60)
    with pytest.raises(ValueError,match="divisor"):
        analyze_transitions([a,a],[0,1],[0,10],60,sector_width_deg=17)


def test_sector_midpoint_not_spatial_dose_contribution():
    a=ap([-5]*3,[5]*3)
    r=analyze_transitions([a,a],[0,1],[28,32],40,sector_width_deg=30)
    assert len(r.sectors)==1
    assert r.sectors[0]["sector_start_deg"]==30
    assert r.summary["angular_assignment"]=="interval_midpoint_shortest_path"
