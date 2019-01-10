"""
    File name: utils
    Author: Artelys
    Creation date: 08/01/2019
    Python Version: 3.6
"""
from nose.tools import assert_true


def assert_almost_intersect(shape1, shape2, msg='', tol=1e-4):
    assert_true(shape1.buffer(tol).intersects(shape2), msg)
