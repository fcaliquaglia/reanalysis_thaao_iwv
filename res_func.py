#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------------
#
"""
Brief description
"""

# =============================================================
# CREATED: 
# AFFILIATION: INGV
# AUTHORS: Filippo Cali' Quaglia
# =============================================================
#
# -------------------------------------------------------------------------------
__author__ = "Filippo Cali' Quaglia"
__credits__ = ["??????"]
__license__ = "GPL"
__version__ = "0.1"
__email__ = "filippo.caliquaglia@ingv.it"
__status__ = "Research"
__lastupdate__ = ""

import pandas as pd

from inputs import *


def data_resampling(vr, var_c, var_e, var_l, var_t, var_t1, var_t2):
    try:
        var_c_res = var_c.resample(tres).mean()
    except (TypeError, NameError):
        var_c_res = pd.DataFrame()
    try:
        var_e_res = var_e.resample(tres).mean()
    except (TypeError, NameError):
        var_e_res = pd.DataFrame()
    try:
        var_l_res = var_l.resample(tres).mean()
    except (TypeError, NameError):
        var_l_res = pd.DataFrame()
    try:
        var_t_res = var_t.resample('1h').mean()
    except (TypeError, NameError):
        var_t_res = pd.DataFrame()
    try:
        var_t1_res = var_t1.resample(tres).mean()
    except (TypeError, NameError):
        var_t1_res = pd.DataFrame()
    try:
        var_t2.index = var_t2.index.round('h')  # oppure var_t2.index.round('h')
        var_t2_res = var_t2.resample(tres_rs).mean()
    except (TypeError, NameError):
        var_t2_res = pd.DataFrame()

    return var_c_res, var_e_res, var_l_res, var_t_res, var_t1_res, var_t2_res
