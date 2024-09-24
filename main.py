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
from typing import Any, Tuple

import julian
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_temp():
    c = pd.DataFrame()
    e = pd.DataFrame()
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
    c[2].name = var

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
    e[2].name = var

    # THAAO
    fn = 'Meteo'
    try:
        t = pd.read_table(
                os.path.join(basefol_t, 'thule_phaao_meteo', fn + '.txt'), skiprows=3, header=[2], engine='python')
        print('OK: ' + fn + '.txt')
    except FileNotFoundError:
        print('NOT FOUND: ' + fn + '.txt')
    t.index = pd.to_datetime(t['DateTime________UTC'], format='%Y-%m-%d %H:%M:%S')
    t.drop(columns=['DateTime________UTC', 'P (hPa)', 'U (%)'], inplace=True)
    t['T (K)'] = t.values - 273.15

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
    # "BP_mbar", "AirTC", "RH", "WS_aws", "WD_aws"
    return [c, e, t, t1, t2]


def read_rh():
    c = pd.DataFrame()
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
    c[2].name = var

    # THAAO
    fn = 'Meteo'
    try:
        t = pd.read_table(
                os.path.join(basefol_t, 'thule_phaao_meteo', fn + '.txt'), skiprows=3, header=[2], engine='python')
        print('OK: ' + fn + '.txt')
    except FileNotFoundError:
        print('NOT FOUND: ' + fn + '.txt')
    t.index = pd.to_datetime(t['DateTime________UTC'], format='%Y-%m-%d %H:%M:%S')
    t.drop(columns=['DateTime________UTC', 'P (hPa)', 'T (K)'], inplace=True)

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
    c[2].name = var

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
    # # "BP_mbar", "AirTC", "RH", "WS_aws", "WD_aws"
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
    c[2].name = var

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
    e[2].name = var

    # THAAO
    fn = 'Meteo'
    try:
        t = pd.read_table(
                os.path.join(basefol_t, 'thule_phaao_meteo', fn + '.txt'), skiprows=3, header=[2], engine='python')
        print('OK: ' + fn + '.txt')
    except FileNotFoundError:
        print('NOT FOUND: ' + fn + '.txt')
    t.index = pd.to_datetime(t['DateTime________UTC'], format='%Y-%m-%d %H:%M:%S')
    t.drop(columns=['DateTime________UTC', 'T (K)', 'U (%)'], inplace=True)

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
    c[2].name = var

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
    e[2].name = var

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
    # e[2].name = var

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
    c[2].name = var

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
    e[2].name = var

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
    t[2].name = var

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
    c[2].name = var

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
    c[2].name = var

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
    c[2].name = var

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
    c[2].name = var

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
    e[2].name = var

    return [c, e, t]


def read_precip():
    c = pd.DataFrame()
    e = pd.DataFrame()
    t = pd.DataFrame()

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
    e[2].name = var

    # # AWS ECAPAC
    # fn = 'DISDRO_THAAO_'
    # for i in pd.date_range(start=dt.datetime(2023, 4, 1), end=dt.datetime(2024, 6, 30), freq='1D'):
    #     try:
    #         file = os.path.join(
    #                 basefol_t, 'thule_phaao_ecapac_disdro_precip',
    #                 fn + i.strftime('%Y_%m_%d') + '_00_00.dat')
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
    # # "BP_mbar", "AirTC", "RH", "WS_aws", "WD_aws"
    # return [c, e, t, t1, t2]

    return [c, e, t]


def read_lwp():
    c = pd.DataFrame()
    e = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()

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
    e[2].name = var

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


c_col = 'red'
e_col = 'blue'
t_col = 'black'
t1_col = 'green'
t2_col = 'purple'

c_col_ori = 'orange'
e_col_ori = 'cyan'
t_col_ori = 'grey'
t1_col_ori = 'lightgreen'
t2_col_ori = 'violet'

tres = '12h'
var_list = ['surf_pres', 'temp', 'alb', 'rh', 'windd', 'winds', 'iwv', 'lwp', 'precip', 'cbh', 'tcc']  # , 'msl_pres']

years = np.arange(2016, 2025, 1)
# years = np.arange(2019, 2022, 1)  # zoom

basefol_c = os.path.join('H:\\Shared drives', 'Dati_elab_docs', 'reanalysis', 'carra', 'thaao')
basefol_e = os.path.join('H:\\Shared drives', 'Dati_elab_docs', 'reanalysis', 'era5', 'thaao')
basefol_t = os.path.join('H:\\Shared drives', 'Dati')

basefol_out = os.path.join('H:\\Shared drives', 'Dati_elab_docs', 'thaao_comparisons')

myFmt = mdates.DateFormatter('%b')

for var in var_list:
    print(var)
    var_c = pd.DataFrame()
    var_e = pd.DataFrame()
    var_t = pd.DataFrame()
    var_t1 = pd.DataFrame()
    var_t2 = pd.DataFrame()
    if var in ['iwv', 'lwp']:
        [var_c, var_e, var_t, var_t1] = read(var)
    elif var in ['temp', 'msl_pres', 'surf_pres', 'rh', 'winds', 'windd']:
        [var_c, var_e, var_t, va_t1, var_t2] = read(var)
    else:
        [var_c, var_e, var_t] = read(var)

    try:
        var_c_res = var_c.resample(tres).mean()
    except TypeError:
        var_c_res = pd.DataFrame()
    try:
        var_e_res = var_e.resample(tres).mean()
    except TypeError:
        var_e_res = pd.DataFrame()
    try:
        var_t_res = var_t.resample(tres).mean()
    except TypeError:
        var_t_res = pd.DataFrame()
    try:
        var_t1_res = var_t1.resample(tres).mean()
    except (TypeError, NameError):
        var_t1 = pd.DataFrame()
        var_t1_res = pd.DataFrame()
    try:
        var_t2_res = var_t2.resample(tres).mean()
    except (TypeError, NameError):
        var_t2 = pd.DataFrame()
        var_t2_res = pd.DataFrame()

    fig, ax = plt.subplots(len(years), 1, figsize=(12, 17), dpi=300)
    for [yy, year] in enumerate(years):
        print('plotting ' + str(year))
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
                    var_t2[var_t2.index.year == year], color=t2_col_ori, label='AWS ECAPAC ori', alpha=0.02, lw=0,
                    marker='.', ms=1)
        except AttributeError:
            pass

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
        if var in ['iwv']:
            ax[yy].set_ylim(0, 30)
        elif var in ['lwp']:
            ax[yy].set_ylim(0, 50)
        elif var in ['temp']:
            ax[yy].set_ylim(-40, 20)
        elif var in ['surf_pres', 'msl_pres']:
            ax[yy].set_ylim(925, 1013)
        elif var in ['winds']:
            ax[yy].set_ylim(0, 30)
        elif var in ['windd']:
            ax[yy].set_ylim(0, 360)
        elif var in ['tcc', 'rh']:
            ax[yy].set_ylim(0, 100)
        elif var in ['alb']:
            ax[yy].set_ylim(0, 1)
            range1 = pd.date_range(dt.datetime(year, 1, 1), dt.datetime(year, 2, 15), freq=tres)
            range2 = pd.date_range(dt.datetime(year, 11, 1), dt.datetime(year, 12, 31), freq=tres)
            ax[yy].vlines(range1.values, 0, 1, color='grey', alpha=0.3)
            ax[yy].vlines(range2.values, 0, 1, color='grey', alpha=0.3)
        else:
            ax[yy].set_ylim(0)
        ax[yy].text(0.45, 0.85, year, transform=ax[yy].transAxes)
        ax[yy].xaxis.set_major_formatter(myFmt)
        ax[yy].set_xlim(dt.datetime(year, 1, 1), dt.datetime(year, 12, 31))

    plt.xlabel('Time')
    plt.legend()
    fig.suptitle(var)
    plt.tight_layout()
    plt.savefig(os.path.join(basefol_out, tres + '_' + f'{var}.png'))
    plt.close('all')
