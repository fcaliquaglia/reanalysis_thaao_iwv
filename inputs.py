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

import os

import matplotlib.dates as mdates
import numpy as np

uom = ' [kg/m2]'

## FOLDERS
basefol_r = os.path.join('H:\\Shared drives', 'Reanalysis')
basefol_t = os.path.join('H:\\Shared drives', 'Dati_THAAO')
basefol_t_elab = os.path.join('H:\\Shared drives', 'Dati_elab_docs')

basefol_c = os.path.join(basefol_r, 'carra', 'thaao', 'v1')
basefol_e = os.path.join(basefol_r, 'era5', 'thaao', 'v1')
basefol_l = os.path.join(basefol_r, 'era5-land', 'thaao', 'v1')

basefol_out = os.path.join(basefol_t_elab, 'thaao_reanalysis')

##
tres = '3h'
tres_rs = '1h'  # only for radiosoundings
var_name = 'iwv'
var_name_u = var_name.upper()

years = np.arange(2016, 2024, 1)

bin_nr=20

SMALL_SIZE = 12

seass = {'all': {'name'      : 'all', 'months': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], 'col': 'pink',
                 'col_CARRA' : 'red', 'col_ERA5': 'blue', 'col_ERA5-L': 'purple', 'col_THAAO': 'grey',
                 'col_HATPRO': 'grey', 'col_VESPA': 'grey', 'col_AWS_ECAPAC': 'purple'},
         'DJF': {'name': 'DJF', 'months': [12, 1, 2], 'col': 'blue'},
         'MAM': {'name': 'MAM', 'months': [3, 4, 5], 'col': 'green'},
         'JJA': {'name': 'JJA', 'months': [6, 7, 8], 'col': 'orange'},
         'SON': {'name': 'SON', 'months': [9, 10, 11], 'col': 'brown'},
         # 'MA' : {'name': 'MA', 'months': [3, 4], 'col': 'yellow'},
         # 'MJ' : {'name': 'MJ', 'months': [5, 6], 'col': 'cyan'},
         # 'JA': {'name': 'JA', 'months': [7, 8], 'col': 'grey'},
         # 'SO' : {'name': 'SO', 'months': [9, 10], 'col': 'purple'}
         }

comps = ['vr_c', 'vr_e', 'vr_t1', 'vr_t2']
var_names = ['vr_c', 'vr_e', 'vr_l', 'vr_t', 'vr_t1', 'vr_t2']

var_dict = {'vr_c': {'name'     : 'vr_c', 'col': 'red', 'col_ori': 'orange', 'label': 'CARRA',
                     'label_uom': f'{var_name} CARRA {uom}'},
    'vr_e': {'name': 'vr_e', 'col': 'blue', 'col_ori': 'cyan', 'label': 'ERA5', 'label_uom': f'{var_name} ERA5 {uom}'},
    'vr_l': {'name'     : 'vr_l', 'col': 'darkgreen', 'col_ori': 'lightgreen', 'label': 'ERA5-L',
             'label_uom': f'{var_name} ERA5-L {uom}'},
    'vr_t': {'name'     : 'vr_t', 'col': 'black', 'col_ori': 'grey', 'label': 'VESPA',
             'label_uom': f'{var_name} VESPA {uom}'},
    'vr_t1': {'name'     : 'vr_t1', 'col': 'green', 'col_ori': 'lightgreen', 'label': 'HATPRO',
              'label_uom': f'{var_name} HATPRO {uom}'},
    'vr_t2': {'name'     : 'vr_t2', 'col': 'purple', 'col_ori': 'violet', 'label': 'RS',
              'label_uom': f'{var_name} RS {uom}'}}

myFmt = mdates.DateFormatter('%b')

extr = {'iwv': {'min': 0, 'max': 20, 'res_min': -5, 'res_max': 5}}
