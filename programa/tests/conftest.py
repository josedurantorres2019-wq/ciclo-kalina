"""Agrega programa/ al sys.path para que los tests importen kalina y nh3h2o."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))