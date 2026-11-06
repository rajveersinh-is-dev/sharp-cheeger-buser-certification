import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from fractions import Fraction
import mpmath as mp
from sharp_constants import (
    CertifiedInterval,
    two_point_ratio,
    sup_two_point_ratio,
    sup_two_point_ratio_certified,
    one_edge_graph,
    sharpness_family
)

def test_certified_interval():
    ci = CertifiedInterval(1.0, 2.0)
    assert ci.lo == 1.0
    assert ci.hi == 2.0
    assert ci.mid == 1.5
    assert ci.width == 1.0
    assert repr(ci) == "[1.0, 2.0]"

def test_two_point_ratio():
    assert abs(two_point_ratio(1.0) - 1.0) < 1e-12
    assert abs(two_point_ratio(0.5) - (2 * 0.5 / 1.5)) < 1e-12

def test_sup_two_point_ratio():
    ci = sup_two_point_ratio()
    assert ci.lo <= 1.0 <= ci.hi
    assert ci.width < 1e-10

def test_sup_two_point_ratio_certified():
    ci = sup_two_point_ratio_certified()
    assert ci.lo <= 2.0 <= ci.hi
    assert ci.width < 1e-10

def test_one_edge_graph():
    W = one_edge_graph()
    assert W == [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]

def test_sharpness_family():
    fam = sharpness_family()
    assert len(fam) == 6
    assert fam[0] == [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]
    assert fam[5] == [[Fraction(5), Fraction(1)], [Fraction(1), Fraction(5)]]
