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

import julian
import pandas as pd

from inputs import *


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
    l = pd.DataFrame()
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

    return [c, e, l, t, t1, t2]


def read_msl_pres():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
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

    return [c, e, l, t, t1, t2]


def read_surf_pres():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
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
    return [c, e, l, t, t1, t2]


def read_alb():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
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

    return [c, e, l, t]


def read_iwv():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
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

    return [c, e, l, t, t1]


def read_winds():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
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

    return [c, e, l, t, t1, t2]


def read_windd():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
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

    return [c, e, l, t, t1, t2]


def read_tcc():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
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

    return [c, e, l, t]


def read_cbh():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
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

    return [c, e, l, t]


def read_precip():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
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

    return [c, e, l, t, t1, t2]


def read_lwp():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()

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

    return [c, e, l, t, t1]


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
