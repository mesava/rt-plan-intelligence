"""Reference-like, completely synthetic VMAT complexity tests.

Checks exact numerical expectations for the chosen jaw-aware MCS variant.
They do NOT establish agreement with an independent clinical implementation.
"""
import numpy as np
import pytest
from rtpi.metrics import Aperture, calculate_beam, control_point_records

def aperture(left, right, bounds=(-15., -5., 5., 15.),
             xjaw=(-20., 20.), yjaw=(-15., 15.)):
    return Aperture(np.array(left, dtype=float), np.array(right, dtype=float),
                    np.array(bounds, dtype=float), xjaw, yjaw)

def test_identical_regular_apertures_yield_unity():
    a = aperture([-5]*3, [5]*3)
    v = calculate_beam(1, "regular", [a,a,a], [0., .3, 1.], 200)
    assert v.mcs_vmat_jaw_clipped == pytest.approx(1.)
    assert v.lsv_cp_mu_weighted == pytest.approx(1.)
    assert v.aav_beam_cp_mu_weighted == pytest.approx(1.)

def test_aav_reduced_at_smaller_control_point():
    wide = aperture([-10]*3,[10]*3)
    narrow = aperture([-5]*3,[5]*3)
    r = calculate_beam(1,"varying",[wide,narrow],[0.,1.],120)
    assert r.lsv_cp_mu_weighted == pytest.approx(1.)
    assert r.aav_beam_cp_mu_weighted == pytest.approx(.75)
    assert r.mcs_vmat_jaw_clipped == pytest.approx(.75)
    records=control_point_records([wide,narrow],[0,1],[0,10])
    assert records[0]["aav_beam_cp"] == pytest.approx(1.)
    assert records[1]["aav_beam_cp"] == pytest.approx(.5)

def test_leaf_irregularity_reduces_lsv_but_not_area_envelope():
    jagged = aperture([-5,-5,-10],[5,5,5])
    r = calculate_beam(1,"jagged",[jagged,jagged],[0.,1.],100)
    assert r.lsv_cp_mu_weighted == pytest.approx(.5)
    assert r.aav_beam_cp_mu_weighted == pytest.approx(1.)
    assert r.mcs_vmat_jaw_clipped == pytest.approx(.5)

def test_unequal_leaf_widths_and_jaw_clipping():
    a = aperture([-10,-10,-10],[10,10,10],
                 bounds=(-12,-8,0,12),xjaw=(-4,4),yjaw=(-10,10))
    b = aperture([-2,-2,-2],[2,2,2],
                 bounds=(-12,-8,0,12),xjaw=(-4,4),yjaw=(-10,10))
    r=calculate_beam(1,"jaws",[a,b],[0,1],80)
    assert a.area_mm2() == pytest.approx(8*20)
    assert b.area_mm2() == pytest.approx(4*20)
    assert r.aav_beam_cp_mu_weighted == pytest.approx(.75)
    assert r.mcs_vmat_jaw_clipped == pytest.approx(.75)

def test_zero_meterset_and_fully_closed_no_fake_values():
    a=aperture([-5]*3,[5]*3)
    b=aperture([0]*3,[0]*3)
    zero=calculate_beam(1,"idle",[a,a],[0,0],0)
    assert zero.mcs_vmat_jaw_clipped is None
    assert zero.lsv_cp_mu_weighted is None
    assert zero.aav_beam_cp_mu_weighted is None
    closed=calculate_beam(1,"closed",[b,b],[0,1],100)
    assert closed.mcs_vmat_jaw_clipped is None

def test_zero_mu_transition_does_not_change_score():
    a=aperture([-10]*3,[10]*3)
    b=aperture([-5]*3,[5]*3)
    r=calculate_beam(1,"idletransition",[a,a,b],[0,0,1],100)
    assert r.mcs_vmat_jaw_clipped == pytest.approx(.75)

def test_single_control_point_is_not_a_vmat_mcs():
    a=aperture([-5]*3,[5]*3)
    r=calculate_beam(1,"static",[a],[1],100)
    assert r.mcs_vmat_jaw_clipped is None

