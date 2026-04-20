"""
Unit and regression test for the CRISPRSCope package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import CRISPRSCope


def test_crisprscope_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "CRISPRSCope" in sys.modules
