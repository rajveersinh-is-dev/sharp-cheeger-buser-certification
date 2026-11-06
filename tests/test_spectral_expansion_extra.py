import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
import pytest
import numpy as np
from fractions import Fraction
import sympy as sp
import mpmath as mp
from spectral_expansion import (
    cheeger_constant, neighbourhood_cheeger_constant,
    normalized_laplacian, lambda2_numeric, lambda2_exact,
    lambda2_certified, verify_bounds_certified, volume
)

def test_volume():
    W = [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]
    assert volume(W, [0]) == Fraction(1)

def test_cheeger_constant_disconnected_graph_or_zero():
    # zero degrees
    W_zero = [[Fraction(0)]]
    assert cheeger_constant(W_zero) is None

def test_neighbourhood_cheeger_constant():
    W = [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]
    assert neighbourhood_cheeger_constant(W) == Fraction(1, 5)

def test_neighbourhood_cheeger_constant_none():
    W_zero = [[Fraction(0)]]
    assert neighbourhood_cheeger_constant(W_zero) is None

def test_normalized_laplacian_none():
    W_zero = [[Fraction(0)]]
    assert normalized_laplacian(W_zero) is None

def test_lambda2_numeric():
    W = [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]
    assert abs(lambda2_numeric(W) - 2.0) < 1e-10

def test_lambda2_numeric_none():
    W_zero = [[Fraction(0)]]
    assert lambda2_numeric(W_zero) is None

def test_lambda2_numeric_disconnected():
    W = [[Fraction(0), Fraction(1), Fraction(0), Fraction(0)],
         [Fraction(1), Fraction(0), Fraction(0), Fraction(0)],
         [Fraction(0), Fraction(0), Fraction(0), Fraction(1)],
         [Fraction(0), Fraction(0), Fraction(1), Fraction(0)]]
    assert lambda2_numeric(W) is None

def test_lambda2_certified_none():
    W_zero = [[Fraction(0)]]
    assert lambda2_certified(W_zero) is None

def test_lambda2_certified_disconnected():
    W = [[Fraction(0), Fraction(1), Fraction(0), Fraction(0)],
         [Fraction(1), Fraction(0), Fraction(0), Fraction(0)],
         [Fraction(0), Fraction(0), Fraction(0), Fraction(1)],
         [Fraction(0), Fraction(0), Fraction(1), Fraction(0)]]
    assert lambda2_certified(W) is None

def test_verify_bounds_certified_degenerate():
    W_zero = [[Fraction(0)]]
    assert verify_bounds_certified(W_zero) is None

def test_cheeger_constant_disconnected_b_zero():
    W = [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)]]
    assert cheeger_constant(W) is None

def test_neighbourhood_cheeger_constant_disconnected_b_zero():
    W = [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)]]
    assert neighbourhood_cheeger_constant(W) is None

def test_lambda2_exact_sort_typeerror(monkeypatch):
    W = [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]

    # Mock sp.N directly to raise TypeError when called
    def mock_N(*args, **kwargs):
        raise TypeError("mock error")

    monkeypatch.setattr(sp, 'N', mock_N)
    assert lambda2_exact(W) is None

def test_lambda2_exact_len_1_zero():
    W = [[Fraction(1)]]
    assert lambda2_exact(W) is None

def test_lambda2_certified_too_few_vals(monkeypatch):
    W = [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]
    def mock_eigsy(*args, **kwargs):
        return [mp.mpf(1)] # Only 1 eigenvalue for a 2x2 matrix
    monkeypatch.setattr(mp, 'eigsy', mock_eigsy)
    assert lambda2_certified(W) is None

def test_neighbourhood_cheeger_constant_vS_zero_or_large():
    # S = {0} has vS = 11, which is > 12/2.
    W = [[Fraction(10), Fraction(1)], [Fraction(1), Fraction(0)]]
    assert neighbourhood_cheeger_constant(W) == Fraction(1, 25)

def test_neighbourhood_cheeger_constant_vS_zero():
    # We want a valid graph where some subset S has vS = 0.
    W = [[Fraction(0), Fraction(0), Fraction(0)],
         [Fraction(0), Fraction(0), Fraction(1)],
         [Fraction(0), Fraction(1), Fraction(0)]]
    assert neighbourhood_cheeger_constant(W) == Fraction(1, 5)
