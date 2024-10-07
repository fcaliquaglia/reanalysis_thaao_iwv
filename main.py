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

import julian
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
from metpy.calc import relative_humidity_from_dewpoint, wind_components, wind_direction, wind_speed
from metpy.units import units

seass = {'all': {'name': 'all', 'months': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]},
         'DJF': {'name': 'DJF', 'months': [12, 1, 2], 'col': 'blue'},
         'MAM': {'name': 'MAM', 'months': [3, 4, 5], 'col': 'green'},
         'JJA': {'name': 'JJA', 'months': [6, 7, 8], 'col': 'orange'},
         'SON': {'name': 'SON', 'months': [9, 10, 11], 'col': 'brown'},
         'MA' : {'name': 'MA', 'months': [3, 4], 'col': 'darkyellow'},
         'MJ' : {'name': 'MJ', 'months': [5, 6], 'col': 'cyan'}, 'JA': {'name': 'JA', 'months': [7, 8], 'col': 'grey'},
         'SO' : {'name': 'SO', 'months': [9, 10], 'col': 'purple'}}

SMALL_SIZE = 12
c_col = 'red'
e_col = 'blue'
l_col = 'purple'
t_col = 'black'
t1_col = 'green'
t2_col = 'purple'

c_col_ori = 'orange'
e_col_ori = 'cyan'
l_col_ori = 'violet'
t_col_ori = 'grey'
t1_col_ori = 'lightgreen'
t2_col_ori = 'violet'

tres = '24h '
var_list = ['temp', 'lwp', 'precip', 'windd', 'winds', 'alb', 'iwv', 'rh', 'surf_pres']  # 'tcc', 'cbh'
# 'cbh'
years = np.arange(2016, 2025, 1)

basefol_c = os.path.join('H:\\Shared drives', 'Dati_elab_docs', 'reanalysis', 'carra', 'thaao')
basefol_e = os.path.join('H:\\Shared drives', 'Dati_elab_docs', 'reanalysis', 'era5', 'thaao')
basefol_l = os.path.join('H:\\Shared drives', 'Dati_elab_docs', 'reanalysis', 'era5-land', 'thaao')
basefol_t = os.path.join('H:\\Shared drives', 'Dati')

basefol_out = os.path.join('H:\\Shared drives', 'Dati_elab_docs', 'thaao_comparisons')

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
        'iwv'      : {'min': 0, 'max': 20, 'res_min': -5, 'res_max': 5}}


def read_temp():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()

    # CARRA
    fn = 'thaao_carra_2m_temperature_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c[2] = c.values - 273.15
    c.columns = [var]

    # ERA5
    fn = 'thaao_era5_2m_temperature_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)

            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e[2] = e[2].values - 273.15
    e.columns = [var]

    # ERA5-L
    fn = 'thaao_era5-land_2m_temperature_'
    for yy, year in enumerate(years):
        try:
            l_tmp = pd.read_table(
                    os.path.join(basefol_l, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            l_tmp[l_tmp == -32767.0] = np.nan
            l = pd.concat([l, l_tmp], axis=0)

            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    l.index = pd.to_datetime(l[0] + ' ' + l[1], format='%Y-%m-%d %H:%M:%S')
    l.drop(columns=[0, 1], inplace=True)
    l[2] = l[2].values - 273.15
    l.columns = [var]

    # THAAO
    import xarray as xr
    fn = 'Meteo_weekly_all'
    try:
        t = xr.open_dataset(os.path.join(basefol_t, 'thule_phaao_meteo', fn + '.nc'), engine='netcdf4').to_dataframe()
        print('OK: ' + fn + '.nc')
    except FileNotFoundError:
        print('NOT FOUND: ' + fn + '.nc')
    t.drop(columns=['BP_hPa', 'RH_%'], inplace=True)
    t['Air_K'] = t.values - 273.15
    t.columns = [var]

    # AWS ECAPAC
    fn = 'AWS_THAAO_'
    for i in pd.date_range(start=dt.datetime(2023, 4, 1), end=dt.datetime(2024, 6, 30), freq='1D'):
        try:
            file = os.path.join(
                    basefol_t, 'thule_phaao_ecapac_aws_snow', 'AWS_ECAPAC', 'Dati_giornalieri_ftp',
                    fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
            t2_tmp = pd.read_csv(
                    file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
                    index_col='TIMESTAMP').iloc[1:, :]
            t2 = pd.concat([t2, t2_tmp], axis=0)
            print('OK: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
    t2.index = pd.DatetimeIndex(t2.index)
    t2.index.name = 'datetime'
    t2 = t2.iloc[:, :].filter(["AirTC"]).astype(float)
    t2.columns = [var]

    return [c, e, l, t, t1, t2]


def read_rh():
    c = pd.DataFrame()
    e_p = pd.DataFrame()
    e_t = pd.DataFrame()
    e_td = pd.DataFrame()
    e = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()

    # CARRA
    fn = 'thaao_carra_2m_relative_humidity_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [var]

    # ERA5
    fn1 = 'thaao_era5_2m_dewpoint_temperature_'
    fn2 = 'thaao_era5_surface_pressure_'
    fn3 = 'thaao_era5_2m_temperature_'
    for yy, year in enumerate(years):
        try:
            e_t_tmp = pd.read_table(
                    os.path.join(basefol_e, fn3 + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_t_tmp[e_t_tmp == -32767.0] = np.nan
            e_t = pd.concat([e_t, e_t_tmp], axis=0)

            print('OK: ' + fn3 + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn3 + str(year) + '.txt')
    e_t.index = pd.to_datetime(e_t[0] + ' ' + e_t[1], format='%Y-%m-%d %H:%M:%S')
    e_t.drop(columns=[0, 1], inplace=True)
    e_t[2].name = var
    e_t.columns = ['e_t']

    for yy, year in enumerate(years):
        try:
            e_td_tmp = pd.read_table(
                    os.path.join(basefol_e, fn1 + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_td[e_td_tmp == -32767.0] = np.nan
            e_td = pd.concat([e_td, e_td_tmp], axis=0)

            print('OK: ' + fn1 + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn1 + str(year) + '.txt')
    e_td.index = pd.to_datetime(e_td[0] + ' ' + e_td[1], format='%Y-%m-%d %H:%M:%S')
    e_td.drop(columns=[0, 1], inplace=True)
    e_td[2].name = var
    e_td.columns = ['e_td']

    e = pd.concat([e_td, e_t], axis=1)

    e['rh'] = relative_humidity_from_dewpoint(e['e_t'].values * units.K, e['e_td'].values * units.K).to('percent')
    e.drop(columns=['e_t', 'e_td'], inplace=True)
    e.columns = [var]

    # THAAO
    import xarray as xr
    fn = 'Meteo_weekly_all'
    try:
        t = xr.open_dataset(os.path.join(basefol_t, 'thule_phaao_meteo', fn + '.nc'), engine='netcdf4').to_dataframe()
        print('OK: ' + fn + '.nc')
    except FileNotFoundError:
        print('NOT FOUND: ' + fn + '.nc')
    t.drop(columns=['BP_hPa', 'Air_K'], inplace=True)
    t.columns = [var]
    #    t.drop(columns=['BP_hPa','Air_K', 'RH_%'], inplace=True)

    # AWS ECAPAC
    fn = 'AWS_THAAO_'
    for i in pd.date_range(start=dt.datetime(2023, 4, 1), end=dt.datetime(2024, 6, 30), freq='1D'):
        try:
            file = os.path.join(
                    basefol_t, 'thule_phaao_ecapac_aws_snow', 'AWS_ECAPAC', 'Dati_giornalieri_ftp',
                    fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
            t2_tmp = pd.read_csv(
                    file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
                    index_col='TIMESTAMP').iloc[1:, :]
            t2 = pd.concat([t2, t2_tmp], axis=0)
            print('OK: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
    t2.index = pd.DatetimeIndex(t2.index)
    t2.index.name = 'datetime'
    t2 = t2.iloc[:, :].filter(["RH"]).astype(float)
    t2.columns = [var]
    # "BP_mbar", "AirTC", "RH", "WS_aws", "WD_aws"
    return [c, e, t, t1, t2]


def read_msl_pres():
    c = pd.DataFrame()
    e = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()

    # CARRA
    fn = 'thaao_carra_mean_sea_level_pressure_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [var]

    # # AWS ECAPAC
    # fn = 'AWS_THAAO_'
    # for i in pd.date_range(start=dt.datetime(2023, 4, 1), end=dt.datetime(2024, 6, 30), freq='1D'):
    #     try:
    #         file = os.path.join(
    #                 basefol_t, 'thule_phaao_ecapac_aws_snow', 'AWS_ECAPAC', 'Dati_giornalieri_ftp', fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
    #         t2_tmp = pd.read_csv(
    #                 file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
    #                 index_col='TIMESTAMP').iloc[1:, :]
    #         t2 = pd.concat([t2, t2_tmp], axis=0)
    #         print('OK: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
    #     except FileNotFoundError:
    #         print('NOT FOUND: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
    # t2.index = pd.DatetimeIndex(t2.index)
    # t2.index.name = 'datetime'
    # t2 = t2.iloc[:, :].filter(["AirTC"]).astype(float)
    # t2.columns = [var]

    return [c, e, t, t1, t2]


def read_surf_pres():
    c = pd.DataFrame()
    e = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()

    # CARRA
    fn = 'thaao_carra_surface_pressure_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c[2] = c.values / 100.
    c.columns = [var]

    # ERA5
    fn = 'thaao_era5_surface_pressure_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e[2] = e.values / 100.
    e.columns = [var]

    # THAAO
    import xarray as xr
    fn = 'Meteo_weekly_all'
    try:
        t = xr.open_dataset(os.path.join(basefol_t, 'thule_phaao_meteo', fn + '.nc'), engine='netcdf4').to_dataframe()
        print('OK: ' + fn + '.nc')
    except FileNotFoundError:
        print('NOT FOUND: ' + fn + '.nc')
    t.drop(columns=['Air_K', 'RH_%'], inplace=True)
    t.columns = [var]

    # AWS ECAPAC
    fn = 'AWS_THAAO_'
    for i in pd.date_range(start=dt.datetime(2023, 4, 1), end=dt.datetime(2024, 6, 30), freq='1D'):
        try:
            file = os.path.join(
                    basefol_t, 'thule_phaao_ecapac_aws_snow', 'AWS_ECAPAC', 'Dati_giornalieri_ftp',
                    fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
            t2_tmp = pd.read_csv(
                    file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
                    index_col='TIMESTAMP').iloc[1:, :]
            t2 = pd.concat([t2, t2_tmp], axis=0)
            print('OK: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
    t2.index = pd.DatetimeIndex(t2.index)
    t2.index.name = 'datetime'
    t2 = t2.iloc[:, :].filter(["BP_mbar"]).astype(float)
    t2.columns = [var]
    # "BP_mbar", "AirTC", "RH", "WS_aws", "WD_aws"
    return [c, e, t, t1, t2]


def read_alb():
    c = pd.DataFrame()
    e = pd.DataFrame()
    t = pd.DataFrame()

    # CARRA
    fn = 'thaao_carra_albedo_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c[2] = c.values / 100.
    c.columns = [var]

    # ERA5
    fn = 'thaao_era5_forecast_albedo_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e.columns = [var]

    # # ERA5
    # fn = 'thaao_era5_snow_albedo_'
    # for yy, year in enumerate(years):
    #     try:
    #         e_tmp = pd.read_table(
    #                 os.path.join(basefol_e, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
    #                 engine='python')
    #         e_tmp[e_tmp == -32767.0] = np.nan
    #         e = pd.concat([e, e_tmp], axis=0)
    #         print('OK: ' + fn + str(year) + '.txt')
    #     except FileNotFoundError:
    #         print('NOT FOUND: ' + fn + str(year) + '.txt')
    # e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    # e.drop(columns=[0, 1], inplace=True)
    # e.columns = [var]

    # THAAO
    fn = 'ALBEDO_SW_'
    for yy, year in enumerate(years):
        try:
            t_tmp = pd.read_table(
                    os.path.join(basefol_t, 'thule_phaao_rad', fn + str(year) + '_5MIN.DAT'), engine='python',
                    skiprows=None, header=0, decimal='.', sep='\s+')
            tmp = np.empty(t_tmp['JDAY_UT'].shape, dtype=dt.datetime)
            for ii, el in enumerate(t_tmp['JDAY_UT']):
                new_jd_ass = el + julian.to_jd(dt.datetime(year - 1, 12, 31, 0, 0), fmt='jd')
                tmp[ii] = julian.from_jd(new_jd_ass, fmt='jd')
                tmp[ii] = tmp[ii].replace(microsecond=0)
            t_tmp.index = pd.DatetimeIndex(tmp)
            t_tmp.drop(['JDAY_UT', 'JDAY_LOC', 'SZA', 'SW_DOWN', 'SW_UP'], axis=1, inplace=True)
            t = pd.concat([t, t_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    t.columns = [var]
    return [c, e, t]


def read_iwv():
    c = pd.DataFrame()
    e = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()

    # CARRA
    fn = 'thaao_carra_total_column_integrated_water_vapour_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [var]

    # ERA5
    fn = 'thaao_era5_total_column_water_vapour_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e.columns = [var]

    # THAAO (vespa)
    fn = 'Vapor_20160712_20221130'
    try:
        t = pd.read_table(
                os.path.join(basefol_t, 'thule_phaao_vespa', fn + '.txt'), skipfooter=1, sep='\s+', header=None,
                skiprows=1, engine='python')
        print('OK: ' + fn + '.txt')
    except FileNotFoundError:
        print('NOT FOUND: ' + fn + str(year) + '.txt')
    t.index = pd.to_datetime(t[0] + ' ' + t[1], format='%Y-%m-%d %H:%M:%S')
    t.drop(columns=[0, 1], inplace=True)
    t.columns = [var]

    # THAAO (hatpro)
    fn = 'QC_IWV_15_min_'
    for yy, year in enumerate(years):
        try:
            t1_tmp = pd.read_table(
                    os.path.join(
                            basefol_t, 'thule_phaao_hatpro', 'definitivi_da_giando', fn + str(year),
                                                                                     fn + str(year) + '.DAT'),
                    sep='\s+', engine='python', header=None, skiprows=1)
            t1_tmp.columns = ['JD_rif', 'IWV', 'STD_IWV', 'RF', 'N']
            tmp = np.empty(t1_tmp['JD_rif'].shape, dtype=dt.datetime)
            for ii, el in enumerate(t1_tmp['JD_rif']):
                new_jd_ass = el + julian.to_jd(dt.datetime(year - 1, 12, 31, 0, 0), fmt='jd')
                tmp[ii] = julian.from_jd(new_jd_ass, fmt='jd')
                tmp[ii] = tmp[ii].replace(microsecond=0)
            t1_tmp.index = pd.DatetimeIndex(tmp)
            t1_tmp.drop(columns=['JD_rif', 'STD_IWV', 'RF', 'N'], axis=1, inplace=True)
            t1 = pd.concat([t1, t1_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.DAT')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.DAT')
    t1['IWV'] = t1['IWV'].values
    t1.columns = [var]

    return [c, e, t, t1]


def read_winds():
    c = pd.DataFrame()
    e = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()

    # CARRA
    fn = 'thaao_carra_10m_wind_speed_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [var]

    # ERA5
    fn_u = 'thaao_era5_10m_u_component_of_wind_'
    fn_v = 'thaao_era5_10m_v_component_of_wind_'
    e_u = pd.DataFrame()
    e_v = pd.DataFrame()
    for yy, year in enumerate(years):
        try:
            e_u_tmp = pd.read_table(
                    os.path.join(basefol_e, fn_u + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None,
                    skiprows=1, engine='python')
            e_u = pd.concat([e_u, e_u_tmp], axis=0)
            print('OK: ' + fn_u + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn_u + str(year) + '.txt')
    e_u.drop(columns=[0, 1], inplace=True)

    for yy, year in enumerate(years):
        try:
            e_v_tmp = pd.read_table(
                    os.path.join(basefol_e, fn_v + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None,
                    skiprows=1, engine='python')
            e_v = pd.concat([e_v, e_v_tmp], axis=0)
            print('OK: ' + fn_v + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn_v + str(year) + '.txt')
    e_v.index = pd.to_datetime(e_v[0] + ' ' + e_v[1], format='%Y-%m-%d %H:%M:%S')
    e_v.drop(columns=[0, 1], inplace=True)

    e_ws = wind_speed(e_u.values * units('m/s'), e_v.values * units('m/s'))

    e.index = e_v.index
    e[var] = e_ws.magnitude
    e.columns = [var]

    # AWS ECAPAC
    fn = 'AWS_THAAO_'
    for i in pd.date_range(start=dt.datetime(2023, 4, 1), end=dt.datetime(2024, 6, 30), freq='1D'):
        try:
            file = os.path.join(
                    basefol_t, 'thule_phaao_ecapac_aws_snow', 'AWS_ECAPAC', 'Dati_giornalieri_ftp',
                    fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
            t2_tmp = pd.read_csv(
                    file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
                    index_col='TIMESTAMP').iloc[1:, :]
            t2 = pd.concat([t2, t2_tmp], axis=0)
            print('OK: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
    t2.index = pd.DatetimeIndex(t2.index)
    t2.index.name = 'datetime'
    t2 = t2.iloc[:, :].filter(["WS_aws"]).astype(float)
    t2.columns = [var]
    # "BP_mbar", "AirTC", "RH", "WS_aws", "WD_aws"
    return [c, e, t, t1, t2]


def read_windd():
    c = pd.DataFrame()
    e = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()

    # CARRA
    fn = 'thaao_carra_10m_wind_direction_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [var]

    # ERA5
    fn_u = 'thaao_era5_10m_u_component_of_wind_'
    fn_v = 'thaao_era5_10m_v_component_of_wind_'
    e_u = pd.DataFrame()
    e_v = pd.DataFrame()
    for yy, year in enumerate(years):
        try:
            e_u_tmp = pd.read_table(
                    os.path.join(basefol_e, fn_u + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None,
                    skiprows=1, engine='python')
            e_u = pd.concat([e_u, e_u_tmp], axis=0)
            print('OK: ' + fn_u + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn_u + str(year) + '.txt')
    e_u.drop(columns=[0, 1], inplace=True)

    for yy, year in enumerate(years):
        try:
            e_v_tmp = pd.read_table(
                    os.path.join(basefol_e, fn_v + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None,
                    skiprows=1, engine='python')
            e_v = pd.concat([e_v, e_v_tmp], axis=0)
            print('OK: ' + fn_v + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn_v + str(year) + '.txt')
    e_v.index = pd.to_datetime(e_v[0] + ' ' + e_v[1], format='%Y-%m-%d %H:%M:%S')
    e_v.drop(columns=[0, 1], inplace=True)

    e_wd = wind_direction(e_u.values * units('m/s'), e_v.values * units('m/s'))
    e.index = e_v.index
    e[var] = e_wd.magnitude
    e.columns = [var]

    # AWS ECAPAC
    fn = 'AWS_THAAO_'
    for i in pd.date_range(start=dt.datetime(2023, 4, 1), end=dt.datetime(2024, 6, 30), freq='1D'):
        try:
            file = os.path.join(
                    basefol_t, 'thule_phaao_ecapac_aws_snow', 'AWS_ECAPAC', 'Dati_giornalieri_ftp',
                    fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
            t2_tmp = pd.read_csv(
                    file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
                    index_col='TIMESTAMP').iloc[1:, :]
            t2 = pd.concat([t2, t2_tmp], axis=0)
            print('OK: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
    t2.index = pd.DatetimeIndex(t2.index)
    t2.index.name = 'datetime'
    t2 = t2.iloc[:, :].filter(["WD_aws"]).astype(float)
    t2.columns = [var]
    # "BP_mbar", "AirTC", "RH", "WS_aws", "WD_aws"
    return [c, e, t, t1, t2]


def read_tcc():
    c = pd.DataFrame()
    e = pd.DataFrame()
    t = pd.DataFrame()

    # CARRA
    fn = 'thaao_carra_total_cloud_cover_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [var]

    # ERA5
    fn = 'thaao_era5_total_cloud_cover_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e[2] = e.values * 100.
    e.columns = [var]

    return [c, e, t]


def read_cbh():
    c = pd.DataFrame()
    e = pd.DataFrame()
    t = pd.DataFrame()

    # CARRA
    fn = 'thaao_carra_cloud_base_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [var]

    # ERA5
    fn = 'thaao_era5_cloud_base_height_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e.columns = [var]

    return [c, e, t]


def read_precip():
    c = pd.DataFrame()
    e = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()

    # CARRA
    fn = 'thaao_carra_total_precipitation_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [var]

    # ERA5
    fn = 'thaao_era5_total_precipitation_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e[2] = e.values * 1000.
    e.columns = [var]

    # AWS ECAPAC
    fn = 'AWS_THAAO_'
    for i in pd.date_range(start=dt.datetime(2023, 4, 1), end=dt.datetime(2024, 6, 30), freq='1D'):
        try:
            file = os.path.join(
                    basefol_t, 'thule_phaao_ecapac_aws_snow', 'AWS_ECAPAC', 'Dati_giornalieri_ftp',
                    fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
            t2_tmp = pd.read_csv(
                    file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
                    index_col='TIMESTAMP').iloc[1:, :]
            t2 = pd.concat([t2, t2_tmp], axis=0)
            print('OK: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
    t2.index = pd.DatetimeIndex(t2.index)
    t2.index.name = 'datetime'
    t2 = t2.iloc[:, :].filter(["PR"]).astype(float)
    t2.columns = [var]
    t2 = t2

    return [c, e, t, t1, t2]


def read_lwp():
    c = pd.DataFrame()
    e = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()

    # CARRA
    fn = 'thaao_carra_total_column_cloud_liquid_water_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c[2] = c.values * 10000
    c.columns = [var]

    # ERA5
    fn = 'thaao_era5_total_column_cloud_liquid_water_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, fn + str(year) + '.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print('OK: ' + fn + str(year) + '.txt')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e[2] = e.values * 1000
    e.columns = [var]

    # THAAO (hatpro)
    fn = 'LWP_15_min_'
    for yy, year in enumerate(years):
        try:
            t1_tmp = pd.read_table(
                    os.path.join(
                            basefol_t, 'thule_phaao_hatpro', 'definitivi_da_giando', fn + str(year) + '_SITO',
                                                                                     fn + str(year) + '_SITO.dat'),
                    sep='\s+', engine='python')
            tmp = np.empty(t1_tmp['JD_rif'].shape, dtype=dt.datetime)
            for ii, el in enumerate(t1_tmp['JD_rif']):
                new_jd_ass = el + julian.to_jd(dt.datetime(year - 1, 12, 31, 0, 0), fmt='jd')
                tmp[ii] = julian.from_jd(new_jd_ass, fmt='jd')
                tmp[ii] = tmp[ii].replace(microsecond=0)
            t1_tmp.index = pd.DatetimeIndex(tmp)
            t1_tmp['LWP_gm-2'] = t1_tmp['LWP_gm-2'].values
            t1_tmp.drop(columns=['JD_rif', 'RF', 'N', 'STD_LWP'], axis=1, inplace=True)
            t1 = pd.concat([t1, t1_tmp], axis=0)

            print('OK: ' + fn + str(year) + '.dat')
        except FileNotFoundError:
            print('NOT FOUND: ' + fn + str(year) + '.dat')
    t1.columns = [var]

    return [c, e, t, t1]


def read(var):
    """

    :param var:
    :return:
    """
    if var == 'temp':
        return read_temp()
    if var == 'rh':
        return read_rh()
    if var == 'surf_pres':
        return read_surf_pres()
    if var == 'msl_pres':
        return read_msl_pres()
    if var == 'iwv':
        return read_iwv()
    if var == 'winds':
        return read_winds()
    if var == 'windd':
        return read_windd()
    if var == 'alb':
        return read_alb()
    if var == 'precip':
        return read_precip()
    if var == 'cbh':
        return read_cbh()
    if var == 'tcc':
        return read_tcc()
    if var == 'lwp':
        return read_lwp()


def precip_res(pr, timeres):
    pr_res = pr.resample(timeres).sum()

    return pr_res


def windd_res(wd, ws, timeres):
    u_df = pd.DataFrame()
    v_df = pd.DataFrame()
    wd_res = pd.DataFrame()

    wind = pd.concat([ws, wd], axis=1)
    wind.columns = ['ws', 'wd']
    u, v = wind_components(wind['ws'].values * units('m/s'), wind['wd'].values * units.deg)
    u_df.index = wind.index
    u_df['u'] = u.magnitude
    v_df.index = wind.index
    v_df['v'] = v.magnitude

    u_df_res = u_df.resample(timeres).mean()
    v_df_res = v_df.resample(timeres).mean()

    wd_res.index = v_df_res.index
    wd_res['windd'] = wind_direction(u_df_res.values * units('m/s'), v_df_res.values * units('m/s')).magnitude

    return wd_res


def plot_ts(period_label):
    fig, ax = plt.subplots(len(years), 1, figsize=(12, 17), dpi=300)
    fig.suptitle(var.upper() + ' all ' + tres, fontweight='bold')
    for [yy, year] in enumerate(years):
        print('plotting ' + str(year))

        # original resolution
        try:
            ax[yy].plot(
                    var_c[var_c.index.year == year], color=c_col_ori, label='CARRA 3h', alpha=0.2, lw=0, marker='.',
                    ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    var_e[var_e.index.year == year], color=e_col_ori, label='ERA5 1h', alpha=0.2, lw=0, marker='.',
                    ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    var_l[var_l.index.year == year], color=l_col_ori, label='ERA5 1h', alpha=0.2, lw=0, marker='.',
                    ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    var_t[var_t.index.year == year], color=t_col_ori, label='THAAO ori', alpha=0.2, lw=0, marker='.',
                    ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    var_t1[var_t1.index.year == year], color=t1_col_ori, label='HATPRO ori', alpha=0.2, lw=0,
                    marker='.', ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    var_t2[var_t2.index.year == year], color=t2_col_ori, label='AWS ECAPAC 1 min', alpha=0.02, lw=0,
                    marker='.', ms=1)
        except AttributeError:
            pass

        # resampled resolution
        try:
            ax[yy].plot(
                    var_c_res[var_c_res.index.year == year], color=c_col, label='CARRA ' + tres, lw=0, marker='.', ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    var_e_res[var_e_res.index.year == year], color=e_col, label='ERA5 ' + tres, lw=0, marker='.', ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    var_l_res[var_l_res.index.year == year], color=l_col, label='ERA5-LAND ' + tres, lw=0, marker='.',
                    ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    var_t_res[var_t_res.index.year == year], color=t_col, label='THAAO ' + tres, lw=0, marker='.', ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    var_t1_res[var_t1_res.index.year == year], color=t1_col, label='HATPRO ' + tres, lw=0, marker='.',
                    ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    var_t2_res[var_t2_res.index.year == year], color=t2_col, label='AWS ECAPAC ' + tres, lw=0,
                    marker='.', ms=2)
        except AttributeError:
            pass

        if var == 'alb':
            range1 = pd.date_range(dt.datetime(year, 1, 1), dt.datetime(year, 2, 15), freq=tres)
            range2 = pd.date_range(dt.datetime(year, 11, 1), dt.datetime(year, 12, 31), freq=tres)
            ax[yy].vlines(range1.values, 0, 1, color='grey', alpha=0.3)
            ax[yy].vlines(range2.values, 0, 1, color='grey', alpha=0.3)
            ax[yy].set_ylim(extr[var]['min'], extr[var]['max'])
        else:
            pass
        ax[yy].set_ylim(extr[var]['min'], extr[var]['max'])
        ax[yy].text(0.45, 0.85, year, transform=ax[yy].transAxes)
        ax[yy].xaxis.set_major_formatter(myFmt)
        ax[yy].set_xlim(dt.datetime(year, 1, 1), dt.datetime(year, 12, 31))
    plt.xlabel('Time')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(basefol_out, tres + '_' + period_label + '_' + f'{var}.png'))
    plt.close('all')


def plot_residuals(period_label):
    fig, ax = plt.subplots(len(years), 1, figsize=(12, 17), dpi=300)
    fig.suptitle('residuals ' + var.upper() + ' all ' + tres, fontweight='bold')
    for [yy, year] in enumerate(years):
        print('plotting ' + str(year))

        if var == 'lwp':
            var_ref = var_t1_res
        elif var in ['precip', 'windd', 'winds']:
            var_ref = var_t2_res
        else:
            var_ref = var_t_res
        # resampled resolution
        daterange = pd.date_range(dt.datetime(year, 1, 1), dt.datetime(year, 12, 31))
        ax[yy].plot(daterange, np.repeat(0, len(daterange)), color='black', lw=2, ls='--')
        try:
            ax[yy].plot(
                    (var_c_res[var_c_res.index.year == year] - var_ref[var_ref.index.year == year]), color=c_col,
                    label='CARRA ' + tres, lw=1, marker='.', ms=0)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    (var_e_res[var_e_res.index.year == year] - var_ref[var_ref.index.year == year]), color=e_col,
                    label='ERA5 ' + tres, lw=1, marker='.', ms=0)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    (var_l_res[var_l_res.index.year == year] - var_ref[var_ref.index.year == year]), color=l_col,
                    label='ERA5-LAND ' + tres, lw=1, marker='.', ms=0)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    (var_t1_res[var_t1_res.index.year == year] - var_ref[var_ref.index.year == year]), color=t1_col,
                    label='HATPRO ' + tres, lw=1, marker='.', ms=0)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    (var_t2_res[var_t2_res.index.year == year] - var_ref[var_ref.index.year == year]), color=t2_col,
                    label='AWS ECAPAC ' + tres, lw=1, marker='.', ms=0)
        except AttributeError:
            pass

        if var == 'alb':
            range1 = pd.date_range(dt.datetime(year, 1, 1), dt.datetime(year, 2, 15), freq=tres)
            range2 = pd.date_range(dt.datetime(year, 11, 1), dt.datetime(year, 12, 31), freq=tres)
            ax[yy].vlines(range1.values, -0.5, 0.5, color='grey', alpha=0.3)
            ax[yy].vlines(range2.values, -0.5, 0.5, color='grey', alpha=0.3)
        else:
            pass

        ax[yy].set_ylim(extr[var]['res_min'], extr[var]['res_max'])
        ax[yy].text(0.45, 0.85, year, transform=ax[yy].transAxes)
        ax[yy].xaxis.set_major_formatter(myFmt)
        ax[yy].set_xlim(dt.datetime(year, 1, 1), dt.datetime(year, 12, 31))
    plt.xlabel('Time')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(basefol_out, tres + '_' + period_label + '_residuals_' + f'{var}.png'))
    plt.close('all')


def plot_scatter(period_label):
    fig, ax = plt.subplots(2, 2, figsize=(12, 12), dpi=300)
    axs = ax.ravel()
    if var == 'lwp':
        comps = ['c', 'e', 't', 't1']
        x = var_t1_res[var]
        xlabel = 'HATPRO'
    elif var in ['windd', 'winds', 'precip']:
        comps = ['c', 'e', 't', 't1']
        x = var_t2_res[var]
        xlabel = 'AWS_ECAPAC'
    else:
        comps = ['c', 'e', 't1', 't2']
        x = var_t_res[var]
        if var == 'iwv':
            xlabel = 'VESPA'
        else:
            xlabel = 'THAAO'

    for i, comp in enumerate(comps):
        axs[i].set_xlabel(xlabel)
        if comp == 'c':
            label = 'CARRA'
            axs[i].set_ylabel(label)
            try:
                y = var_c_res[var]
            except KeyError:
                print('error with ' + label)
                continue
        if comp == 'e':
            label = 'ERA5'
            axs[i].set_ylabel(label)
            try:
                y = var_e_res[var]
            except KeyError:
                print('error with ' + label)
                continue
        if comp == 't':
            label = 'THAAO'
            axs[i].set_ylabel(label)
            try:
                y = var_t_res[var]
            except KeyError:
                print('error with ' + label)
                continue
        if comp == 't1':
            label = 'HATPRO'
            axs[i].set_ylabel(label)
            try:
                y = var_t1_res[var]
            except KeyError:
                print('error with ' + label)
                continue
        if comp == 't2':
            label = 'AWS ECAPAC'
            axs[i].set_ylabel(label)
            try:
                y = var_t2_res[var]
            except KeyError:
                print('error with ' + label)
                continue
        try:
            print('plotting scatter THAAO-' + label)

            fig.suptitle(var.upper() + ' ' + seass[period_label]['name'] + ' ' + tres, fontweight='bold')
            axs[i].set_title(label)

            time_list = pd.date_range(start=dt.datetime(2016, 1, 1), end=dt.datetime(2024, 12, 31), freq=tres)
            if x.empty | y.empty:
                continue
            x_all = x.reindex(time_list)
            x = x_all.loc[(x_all.index.month.isin(seass[period_label]['months']))]
            y_all = y.reindex(time_list)
            y = y_all.loc[(y_all.index.month.isin(seass[period_label]['months']))]

            idx = np.isfinite(x) & np.isfinite(y)

            if seass[period_label]['name'] != 'all':
                axs[i].scatter(x[idx], y[idx], color=seass[period_label]['col'])
            else:
                axs[i].hist2d(x[idx], y[idx], bins=(100, 100), cmap=plt.cm.jet, cmin=1)

            b, a = np.polyfit(x[idx], y[idx], deg=1)
            xseq = np.linspace(extr[var]['min'], extr[var]['max'], num=1000)
            axs[i].plot(xseq, a + b * xseq, color='red', lw=2.5, ls='--')
            axs[i].plot(
                    [extr[var]['min'], extr[var]['max']], [extr[var]['min'], extr[var]['max']], color='black', lw=1.5,
                    ls='-')
            corcoef = ma.corrcoef(x[idx], y[idx])

            N = x[idx].shape[0]
            rmse = np.sqrt(np.sum((x[idx] - y[idx]) ** 2) / N)
            bias = np.nanmean(x[idx] - y[idx])
            axs[i].text(
                    0.60, 0.15, 'R=' + f"{corcoef[0, 1]:1.3}" + '\nrmse=' + f"{rmse:1.3}" + '\nN=' + str(
                            N) + '\nbias=' + f"{bias:1.3}", fontsize=14, transform=axs[i].transAxes)
            axs[i].set_xlim(extr[var]['min'], extr[var]['max'])
            axs[i].set_ylim(extr[var]['min'], extr[var]['max'])
        except:
            print('error with ' + label)
    plt.savefig(os.path.join(basefol_out, tres + '_scatter_' + seass[period_label]['name'] + '_' + f'{var}.png'))
    plt.close('all')


def data_resampling(vr, var_c, var_e, var_l, var_t, var_t1, var_t2):
    if vr == 'windd':
        [var_c_ws, var_e_ws, var_l_ws, var_t_ws, var_t1_ws, var_t2_ws] = read('winds')
        try:
            var_c_res = windd_res(var_c, var_c_ws, tres)
        except (TypeError, NameError, ValueError):
            var_c_res = pd.DataFrame()
        try:
            var_e_res = windd_res(var_e, var_e_ws, tres)
        except (TypeError, NameError, ValueError):
            var_e_res = pd.DataFrame()
        try:
            var_l_res = windd_res(var_l, var_l_ws, tres)
        except (TypeError, NameError, ValueError):
            var_l_res = pd.DataFrame()
        try:
            var_t_res = windd_res(var_t, var_t_ws, tres)
        except (TypeError, NameError, ValueError):
            var_t_res = pd.DataFrame()
        try:
            var_t1_res = windd_res(var_t1, var_t1_ws, tres)
        except (TypeError, NameError, ValueError):
            var_t1_res = pd.DataFrame()
        try:
            var_t2_res = windd_res(var_t2, var_t2_ws, tres)
        except (TypeError, NameError, ValueError):
            var_t2_res = pd.DataFrame()
    elif vr == 'precip':
        try:
            var_c_res = precip_res(var_c, tres)
        except (TypeError, NameError, ValueError):
            var_c_res = pd.DataFrame()
        try:
            var_e_res = precip_res(var_e, tres)
        except (TypeError, NameError, ValueError):
            var_e_res = pd.DataFrame()
        try:
            var_l_res = precip_res(var_l, tres)
        except:
            var_l_res = pd.DataFrame()
        try:
            var_t_res = precip_res(var_t, tres)
        except (TypeError, NameError, ValueError):
            var_t_res = pd.DataFrame()
        try:
            var_t1_res = precip_res(var_t1, tres)
        except (TypeError, NameError, ValueError):
            var_t1_res = pd.DataFrame()
        try:
            var_t2_res = precip_res(var_t2, tres)
        except (TypeError, NameError, ValueError):
            var_t2_res = pd.DataFrame()
    else:
        try:
            var_c_res = var_c.resample(tres).mean()
        except TypeError:
            var_c_res = pd.DataFrame()
        try:
            var_e_res = var_e.resample(tres).mean()
        except TypeError:
            var_e_res = pd.DataFrame()
        try:
            var_l_res = var_l.resample(tres).mean()
        except:
            var_l_res = pd.DataFrame()
        try:
            var_t_res = var_t.resample(tres).mean()
        except:
            var_t_res = pd.DataFrame()
        try:
            var_t1_res = var_t1.resample(tres).mean()
        except (TypeError, NameError):

            var_t1_res = pd.DataFrame()
        try:
            var_t2_res = var_t2.resample(tres).mean()
        except (TypeError, NameError):

            var_t2_res = pd.DataFrame()

    return var_c_res, var_e_res, var_l_res, var_t_res, var_t1_res, var_t2_res


if __name__ == "__main__":

    for var in var_list:
        print(var)
        var_c = pd.DataFrame()
        var_e = pd.DataFrame()
        var_l = pd.DataFrame()
        var_t = pd.DataFrame()
        var_t1 = pd.DataFrame()
        var_t2 = pd.DataFrame()
        if var in ['iwv', 'lwp']:
            [var_c, var_e, var_l, var_t, var_t1] = read(var)
        elif var in ['temp', 'msl_pres', 'surf_pres', 'rh', 'winds', 'windd', 'precip']:
            [var_c, var_e, var_l, var_t, var_t1, var_t2] = read(var)
        else:
            [var_c, var_e, var_l, var_t] = read(var)

        # time RESAMPLING (specific for windd-->using wind components, and precip--> cumulative)
        var_c_res, var_e_res, var_l_res, var_t_res, var_t1_res, var_t2_res = data_resampling(
                var, var_c, var_e, var_l, var_t, var_t1, var_t2)

        plot_ts('all')
        plot_residuals('all')
        for seas in seass:
            plot_scatter(seas)
