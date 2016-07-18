#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Datasets and parameter sets for testing physical property estimation.

Authors
-------
* John D. Chodera <john.chodera@choderalab.org>

TODO
----
* Implement methods

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import sys
import time
import copy
import numpy as np

from simtk import openmm, unit

#=============================================================================================
# TEST DATASETS
#=============================================================================================

from openforcefield import PhysicalPropertyDataset

TestDataset = PhysicalPropertyDataset()

#=============================================================================================
# TEST PARAMETER SETS
#=============================================================================================

from openforcefield import ParameterSet

TestParameterSet = ParameterSet()
