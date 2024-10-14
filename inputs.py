#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------------
#
"""
Brief description
"""

# =============================================================
# CREATED: 
# AFFILIATION: UNIVE, INGV
# AUTHORS: Filippo Cali' Quaglia
# =============================================================
#
# -------------------------------------------------------------------------------
__author__ = "Filippo Cali' Quaglia"
__credits__ = ["??????"]
__license__ = "GPL"
__version__ = "0.1"
__email__ = "filippo.caliquaglia@gmail.com"
__status__ = "Research"
__lastupdate__ = ""

import datetime as dt
import os

import matplotlib.dates as mdates
import numpy as np
import pandas as pd

## FOLDERS
basefol_c = os.path.join('H:\\Shared drives', 'Dati_elab_docs', 'reanalysis', 'carra', 'thaao')
basefol_e = os.path.join('H:\\Shared drives', 'Dati_elab_docs', 'reanalysis', 'era5', 'thaao')
basefol_l = os.path.join('H:\\Shared drives', 'Dati_elab_docs', 'reanalysis', 'era5-land', 'thaao')
basefol_t = os.path.join('H:\\Shared drives', 'Dati')
basefol_out = os.path.join('H:\\Shared drives', 'Dati_elab_docs', 'thaao_comparisons')

##
tres = '24h'
var_list = ['temp', 'sw_down', 'sw_up', 'lw_up', 'lw_down', 'alb', 'lwp', 'iwv', 'precip', 'windd', 'winds', 'alb',
            'iwv', 'rh', 'surf_pres']
# 'tcc', 'cbh'

years = np.arange(2016, 2025, 1)

seass = {'all': {'name': 'all', 'months': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]},
         'DJF': {'name': 'DJF', 'months': [12, 1, 2], 'col': 'blue'},
         'MAM': {'name': 'MAM', 'months': [3, 4, 5], 'col': 'green'},
         'JJA': {'name': 'JJA', 'months': [6, 7, 8], 'col': 'orange'},
         'SON': {'name': 'SON', 'months': [9, 10, 11], 'col': 'brown'},
         'MA' : {'name': 'MA', 'months': [3, 4], 'col': 'yellow'},
         'MJ' : {'name': 'MJ', 'months': [5, 6], 'col': 'cyan'}, 'JA': {'name': 'JA', 'months': [7, 8], 'col': 'grey'},
         'SO' : {'name': 'SO', 'months': [9, 10], 'col': 'purple'}}

SMALL_SIZE = 12
c_col = 'red'
e_col = 'blue'
l_col = 'darkgreen'
t_col = 'black'
t1_col = 'green'
t2_col = 'purple'

c_col_ori = 'orange'
e_col_ori = 'cyan'
l_col_ori = 'lightgreen'
t_col_ori = 'grey'
t1_col_ori = 'lightgreen'
t2_col_ori = 'violet'

myFmt = mdates.DateFormatter('%d-%b')

extr = {'temp'     : {'min': -40, 'max': 20, 'res_min': -10, 'res_max': 10},
        'lwp'      : {'min': 0, 'max': 50, 'res_min': -20, 'res_max': 20},
        'rh'       : {'min': 0, 'max': 100, 'res_min': -10, 'res_max': 10},
        'tcc'      : {'min': 0, 'max': 100, 'res_min': -10, 'res_max': 10},
        'windd'    : {'min': 0, 'max': 360, 'res_min': -90, 'res_max': 90},
        'winds'    : {'min': 0, 'max': 30, 'res_min': -10, 'res_max': 10},
        'precip'   : {'min': 0, 'max': 25, 'res_min': -10, 'res_max': 10},
        'surf_pres': {'min': 925, 'max': 1013, 'res_min': -20, 'res_max': 20},
        'alb'      : {'min': 0, 'max': 1, 'res_min': -0.5, 'res_max': 0.5},
        'iwv'      : {'min': 0, 'max': 20, 'res_min': -5, 'res_max': 5},
        'lw_up'    : {'min': 0, 'max': 500, 'res_min': -20, 'res_max': 20},
        'lw_down'  : {'min': 0, 'max': 500, 'res_min': -20, 'res_max': 20},
        'sw_up'    : {'min': 0, 'max': 500, 'res_min': -20, 'res_max': 20},
        'sw_down'  : {'min': 0, 'max': 700, 'res_min': -20, 'res_max': 20}}

aws_ecapac_daterange = pd.date_range(start=dt.datetime(2023, 4, 1), end=dt.datetime(2024, 12, 31), freq='1D')
