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
import xarray as xr
from metpy.calc import relative_humidity_from_dewpoint, wind_direction, wind_speed
from metpy.units import units

from inputs import *


def read_temp():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()
    vr = 'temp'

    # CARRA
    fn = 'thaao_carra_2m_temperature_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c[2] = c.values - 273.15
    c.columns = [vr]

    # ERA5
    fn = 'thaao_era5_2m_temperature_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)

            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e[2] = e[2].values - 273.15
    e.columns = [vr]

    # ERA5-L
    fn = 'thaao_era5-land_2m_temperature_'
    for yy, year in enumerate(years):
        try:
            l_tmp = pd.read_table(
                    os.path.join(basefol_l, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            l_tmp[l_tmp == -32767.0] = np.nan
            l = pd.concat([l, l_tmp], axis=0)

            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    l.index = pd.to_datetime(l[0] + ' ' + l[1], format='%Y-%m-%d %H:%M:%S')
    l.drop(columns=[0, 1], inplace=True)
    l[2] = l[2].values - 273.15
    l.columns = [vr]

    # THAAO
    fn = 'Meteo_weekly_all'
    try:
        t = xr.open_dataset(os.path.join(basefol_t, 'thaao_meteo', f'{fn}.nc'), engine='netcdf4').to_dataframe()
        print(f'OK: {fn}.nc')
    except FileNotFoundError:
        print(f'NOT FOUND: {fn}.nc')
    t.drop(columns=['BP_hPa', 'RH_%'], inplace=True)
    t['Air_K'] = t.values - 273.15
    t.columns = [vr]

    # AWS ECAPAC
    fn = 'AWS_THAAO_'
    for i in aws_ecapac_daterange:
        i_fmt = i.strftime('%Y_%m_%d')
        try:
            file = os.path.join(
                    basefol_t, 'thaao_ecapac_aws_snow', 'AWS_ECAPAC', f'{fn}{i_fmt}_00_00.dat')
            t2_tmp = pd.read_csv(
                    file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
                    index_col='TIMESTAMP').iloc[1:, :]
            t2 = pd.concat([t2, t2_tmp], axis=0)
            print(f'OK: {fn}{i_fmt}_00_00.dat')
        except (FileNotFoundError, pd.errors.EmptyDataError):
            print(f'NOT_FOUND: {fn}{i_fmt}_00_00.dat')
    t2.index = pd.DatetimeIndex(t2.index)
    t2.index.name = 'datetime'
    t2 = t2.iloc[:, :].filter(["AirTC"]).astype(float)
    t2.columns = [vr]

    return [c, e, l, t, t1, t2]


def read_rh():
    c = pd.DataFrame()
    e_t = pd.DataFrame()
    e_td = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()
    vr = 'rh'

    # CARRA
    fn = 'thaao_carra_2m_relative_humidity_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [vr]

    # ERA5
    fn1 = 'thaao_era5_2m_dewpoint_temperature_'
    fn3 = 'thaao_era5_2m_temperature_'
    for yy, year in enumerate(years):
        try:
            e_t_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn3}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_t_tmp[e_t_tmp == -32767.0] = np.nan
            e_t = pd.concat([e_t, e_t_tmp], axis=0)

            print(f'OK: {fn3}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn3}{year}.txt')
    e_t.index = pd.to_datetime(e_t[0] + ' ' + e_t[1], format='%Y-%m-%d %H:%M:%S')
    e_t.drop(columns=[0, 1], inplace=True)
    e_t[2].name = vr
    e_t.columns = ['e_t']

    for yy, year in enumerate(years):
        try:
            e_td_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn1}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_td[e_td_tmp == -32767.0] = np.nan
            e_td = pd.concat([e_td, e_td_tmp], axis=0)

            print(f'OK: {fn1}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn1}{year}.txt')
    e_td.index = pd.to_datetime(e_td[0] + ' ' + e_td[1], format='%Y-%m-%d %H:%M:%S')
    e_td.drop(columns=[0, 1], inplace=True)
    e_td[2].name = vr
    e_td.columns = ['e_td']

    e = pd.concat([e_td, e_t], axis=1)

    e['rh'] = relative_humidity_from_dewpoint(e['e_t'].values * units.K, e['e_td'].values * units.K).to('percent')
    e.drop(columns=['e_t', 'e_td'], inplace=True)
    e.columns = [vr]

    # THAAO
    fn = 'Meteo_weekly_all'
    try:
        t = xr.open_dataset(os.path.join(basefol_t, 'thaao_meteo', f'{fn}.nc'), engine='netcdf4').to_dataframe()
        print(f'OK: {fn}.nc')
    except FileNotFoundError:
        print(f'NOT FOUND: {fn}.nc')
    t.drop(columns=['BP_hPa', 'Air_K'], inplace=True)
    t.columns = [vr]
    #    t.drop(columns=['BP_hPa','Air_K', 'RH_%'], inplace=True)

    # AWS ECAPAC
    fn = 'AWS_THAAO_'
    for i in aws_ecapac_daterange:
        i_fmt = i.strftime('%Y_%m_%d')
        try:
            file = os.path.join(
                    basefol_t, 'thaao_ecapac_aws_snow', 'AWS_ECAPAC', f'{fn}{i_fmt}_00_00.dat')
            t2_tmp = pd.read_csv(
                    file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
                    index_col='TIMESTAMP').iloc[1:, :]
            t2 = pd.concat([t2, t2_tmp], axis=0)
            print(f'OK: {fn}{i_fmt}_00_00.dat')
        except (FileNotFoundError, pd.errors.EmptyDataError):
            print(f'NOT_FOUND: {fn}{i_fmt}_00_00.dat')
    t2.index = pd.DatetimeIndex(t2.index)
    t2.index.name = 'datetime'
    t2 = t2.iloc[:, :].filter(["RH"]).astype(float)
    t2.columns = [vr]

    return [c, e, l, t, t1, t2]


def read_msl_pres():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()
    vr = 'msl_pres'

    # CARRA
    fn = 'thaao_carra_mean_sea_level_pressure_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [vr]

    # # AWS ECAPAC
    # fn = 'AWS_THAAO_'
    # for i in aws_ecapac_daterange:
    # i_tmp = i.strftime('%Y_%m_%d')
    #     try:
    #         file = os.path.join(
    #                 basefol_t, 'thaao_ecapac_aws_snow', 'AWS_ECAPAC', f'{fn}{i_tmp}_00_00.dat')
    #         t2_tmp = pd.read_csv(
    #                 file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
    #                 index_col='TIMESTAMP').iloc[1:, :]
    #         t2 = pd.concat([t2, t2_tmp], axis=0)
    #         print(f'OK: {fn}{i_tmp}_00_00.dat')
    #     except FileNotFoundError:
    #         print(f'NOT_FOUND: {fn}{i_tmp}_00_00.dat')
    # t2.index = pd.DatetimeIndex(t2.index)
    # t2.index.name = 'datetime'
    # t2 = t2.iloc[:, :].filter(["AirTC"]).astype(float)
    # t2.columns = [vr]

    return [c, e, l, t, t1, t2]


def read_surf_pres():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()
    vr = 'surf_pres'
    # CARRA
    fn = 'thaao_carra_surface_pressure_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
            c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
            c.drop(columns=[0, 1], inplace=True)
            c[2] = c.values / 100.
            c.columns = [vr]

            # ERA5
            fn = 'thaao_era5_surface_pressure_'
            for yy, year in enumerate(years):
                try:
                    e_tmp = pd.read_table(
                            os.path.join(basefol_e, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None,
                            skiprows=1, engine='python')
                    e_tmp[e_tmp == -32767.0] = np.nan
                    e = pd.concat([e, e_tmp], axis=0)
                    print(f'OK: {fn}{year}.txt')
                except FileNotFoundError:
                    print(f'NOT FOUND: {fn}{year}.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e[2] = e.values / 100.
    e.columns = [vr]

    # THAAO
    fn = 'Meteo_weekly_all'
    try:
        t = xr.open_dataset(os.path.join(basefol_t, 'thaao_meteo', f'{fn}.nc'), engine='netcdf4').to_dataframe()
        print(f'OK: {fn}.nc')
    except FileNotFoundError:
        print(f'NOT FOUND: {fn}.nc')
    t.drop(columns=['Air_K', 'RH_%'], inplace=True)
    t.columns = [vr]

    # AWS ECAPAC
    fn = 'AWS_THAAO_'
    for i in aws_ecapac_daterange:
        i_fmt = i.strftime('%Y_%m_%d')
        try:
            file = os.path.join(
                    basefol_t, 'thaao_ecapac_aws_snow', 'AWS_ECAPAC', f'{fn}{i_fmt}_00_00.dat')
            t2_tmp = pd.read_csv(
                    file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
                    index_col='TIMESTAMP').iloc[1:, :]
            t2 = pd.concat([t2, t2_tmp], axis=0)
            print(f'OK: {fn}{i_fmt}_00_00.dat')
        except (FileNotFoundError, pd.errors.EmptyDataError):
            print(f'NOT_FOUND: {fn}{i_fmt}_00_00.dat')
    t2.index = pd.DatetimeIndex(t2.index)
    t2.index.name = 'datetime'
    t2 = t2.iloc[:, :].filter(["BP_mbar"]).astype(float)
    t2.columns = [vr]

    return [c, e, l, t, t1, t2]


def read_alb():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()
    vr = 'alb'

    # CARRA
    fn = 'thaao_carra_albedo_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c[2] = c.values / 100.
    c.columns = [vr]
    c[c <= 0.1] = np.nan

    # ERA5
    fn = 'thaao_era5_forecast_albedo_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e.columns = [vr]
    e[e <= 0.1] = np.nan

    # ERA5
    fn = 'thaao_era5_snow_albedo_'
    for yy, year in enumerate(years):
        try:
            t2_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            t2_tmp[t2_tmp == -32767.0] = np.nan
            t2 = pd.concat([t2, t2_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    t2.index = pd.to_datetime(t2[0] + ' ' + t2[1], format='%Y-%m-%d %H:%M:%S')
    t2.drop(columns=[0, 1], inplace=True)
    t2.columns = [vr]
    t2[t2 <= 0.1] = np.nan

    # THAAO
    # TODO: sostituire con questo blocco che prende direttamente dal file MERGED_SW_LW_UP_DW_METEO_YYYY.dat

    # fn = 'MERGED_SW_LW_UP_DW_METEO_'
    # for yy, year in enumerate(years):
    #     try:
    #         t_tmp = pd.read_table(
    #                 os.path.join(basefol_t, 'thaao_rad', f'{fn}{year}_5MIN.dat'), engine='python',
    #                 skiprows=None, header=0, decimal='.', sep='\s+')
    #         tmp = np.empty(t_tmp['JDAY_UT'].shape, dtype=dt.datetime)
    #         for ii, el in enumerate(t_tmp['JDAY_UT']):
    #             new_jd_ass = el + julian.to_jd(dt.datetime(year - 1, 12, 31, 0, 0), fmt='jd')
    #             tmp[ii] = julian.from_jd(new_jd_ass, fmt='jd')
    #             tmp[ii] = tmp[ii].replace(microsecond=0)
    #         t_tmp.index = pd.DatetimeIndex(tmp)
    #         t_tmp.drop(
    #                 ['JDAY_UT', 'JDAY_LOC', 'SZA', 'SW_DOWN', 'SW_UP','PAR_DOWN', 'PAR_UP', 'LW_DOWN', 'LW_UP', 'TBP',
    #                  'ALBEDO_LW', 'ALBEDO_PAR', 'P', 'T', 'RH', 'PE', 'RR2'], axis=1, inplace=True)
    #         t = pd.concat([t, t_tmp], axis=0)
    #         print(f'OK: {fn}{year}.txt')
    #     except FileNotFoundError:
    #         print(f'NOT FOUND: {fn}{year}.txt')
    # t.columns = [vr]

    fn = 'ALBEDO_SW_'
    for yy, year in enumerate(years):
        try:
            t_tmp = pd.read_table(
                    os.path.join(basefol_t, 'thaao_rad', f'{fn}{year}_5MIN.DAT'), engine='python', skiprows=None,
                    header=0, decimal='.', sep='\s+')
            tmp = np.empty(t_tmp['JDAY_UT'].shape, dtype=dt.datetime)
            for ii, el in enumerate(t_tmp['JDAY_UT']):
                new_jd_ass = el + julian.to_jd(dt.datetime(year - 1, 12, 31, 0, 0), fmt='jd')
                tmp[ii] = julian.from_jd(new_jd_ass, fmt='jd')
                tmp[ii] = tmp[ii].replace(microsecond=0)
            t_tmp.index = pd.DatetimeIndex(tmp)
            t_tmp.drop(['JDAY_UT', 'JDAY_LOC', 'SZA', 'SW_DOWN', 'SW_UP'], axis=1, inplace=True)
            t = pd.concat([t, t_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    t.columns = [vr]
    t[t <= 0.1] = np.nan

    return [c, e, l, t, t1, t2]


def read_iwv():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    vr = 'iwv'

    # CARRA
    fn = 'thaao_carra_total_column_integrated_water_vapour_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [vr]

    # ERA5
    fn = 'thaao_era5_total_column_water_vapour_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e.columns = [vr]

    # THAAO (vespa)
    fn = 'Vapor_20160712_20221130'
    try:
        t = pd.read_table(
                os.path.join(basefol_t, 'thaao_vespa', f'{fn}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                engine='python')
        print(f'OK: {fn}.txt')
    except FileNotFoundError:
        print(f'NOT FOUND: {fn}{year}.txt')
    t.index = pd.to_datetime(t[0] + ' ' + t[1], format='%Y-%m-%d %H:%M:%S')
    t.drop(columns=[0, 1], inplace=True)
    t.columns = [vr]

    # THAAO (hatpro)
    fn = 'QC_IWV_15_min_'
    for yy, year in enumerate(years):
        try:
            t1_tmp = pd.read_table(
                    os.path.join(
                            basefol_t, 'thaao_hatpro', 'definitivi_da_giando', f'{fn}{year}', f'{fn}{year}.DAT'),
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
            print(f'OK: {fn}{year}.DAT')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.DAT')
    t1['IWV'] = t1['IWV'].values
    t1.columns = [vr]
    # cleaning HATPRO DATA
    t1[t1 < 0] = np.nan

    return [c, e, l, t, t1]


def read_winds():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()
    vr = 'winds'

    # CARRA
    fn = 'thaao_carra_10m_wind_speed_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [vr]

    # ERA5
    fn_u = 'thaao_era5_10m_u_component_of_wind_'
    fn_v = 'thaao_era5_10m_v_component_of_wind_'
    e_u = pd.DataFrame()
    e_v = pd.DataFrame()
    for yy, year in enumerate(years):
        try:
            e_u_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn_u}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_u = pd.concat([e_u, e_u_tmp], axis=0)
            print(f'OK: {fn_u}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn_u}{year}.txt')
    e_u.drop(columns=[0, 1], inplace=True)

    for yy, year in enumerate(years):
        try:
            e_v_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn_v}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_v = pd.concat([e_v, e_v_tmp], axis=0)
            print(f'OK: {fn_v}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn_v}{year}.txt')
    e_v.index = pd.to_datetime(e_v[0] + ' ' + e_v[1], format='%Y-%m-%d %H:%M:%S')
    e_v.drop(columns=[0, 1], inplace=True)

    e_ws = wind_speed(e_u.values * units('m/s'), e_v.values * units('m/s'))

    e.index = e_v.index
    e[vr] = e_ws.magnitude
    e.columns = [vr]

    # AWS ECAPAC
    fn = 'AWS_THAAO_'
    for i in aws_ecapac_daterange:
        i_fmt = i.strftime('%Y_%m_%d')
        try:
            file = os.path.join(
                    basefol_t, 'thaao_ecapac_aws_snow', 'AWS_ECAPAC', f'{fn}{i_fmt}_00_00.dat')
            t2_tmp = pd.read_csv(
                    file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
                    index_col='TIMESTAMP').iloc[1:, :]
            t2 = pd.concat([t2, t2_tmp], axis=0)
            print(f'OK: {fn}{i_fmt}_00_00.dat')
        except (FileNotFoundError, pd.errors.EmptyDataError):
            print(f'NOT_FOUND: {fn}{i_fmt}_00_00.dat')
    t2.index = pd.DatetimeIndex(t2.index)
    t2.index.name = 'datetime'
    t2 = t2.iloc[:, :].filter(["WS_aws"]).astype(float)
    t2.columns = [vr]

    return [c, e, l, t, t1, t2]


def read_windd():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()
    vr = 'windd'

    # CARRA
    fn = 'thaao_carra_10m_wind_direction_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [vr]

    # ERA5
    fn_u = 'thaao_era5_10m_u_component_of_wind_'
    fn_v = 'thaao_era5_10m_v_component_of_wind_'
    e_u = pd.DataFrame()
    e_v = pd.DataFrame()
    for yy, year in enumerate(years):
        try:
            e_u_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn_u}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_u = pd.concat([e_u, e_u_tmp], axis=0)
            print(f'OK: {fn_u}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn_u}{year}.txt')
    e_u.drop(columns=[0, 1], inplace=True)

    for yy, year in enumerate(years):
        try:
            e_v_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn_v}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_v = pd.concat([e_v, e_v_tmp], axis=0)
            print(f'OK: {fn_v}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn_v}{year}.txt')
    e_v.index = pd.to_datetime(e_v[0] + ' ' + e_v[1], format='%Y-%m-%d %H:%M:%S')
    e_v.drop(columns=[0, 1], inplace=True)

    e_wd = wind_direction(e_u.values * units('m/s'), e_v.values * units('m/s'))
    e.index = e_v.index
    e[vr] = e_wd.magnitude
    e.columns = [vr]

    # AWS ECAPAC
    fn = 'AWS_THAAO_'
    for i in aws_ecapac_daterange:
        i_fmt = i.strftime('%Y_%m_%d')
        try:
            file = os.path.join(
                    basefol_t, 'thaao_ecapac_aws_snow', 'AWS_ECAPAC', f'{fn}{i_fmt}_00_00.dat')
            t2_tmp = pd.read_csv(
                    file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
                    index_col='TIMESTAMP').iloc[1:, :]
            t2 = pd.concat([t2, t2_tmp], axis=0)
            print(f'OK: {fn}{i_fmt}_00_00.dat')
        except (FileNotFoundError, pd.errors.EmptyDataError):
            print(f'NOT_FOUND: {fn}{i_fmt}_00_00.dat')
    t2.index = pd.DatetimeIndex(t2.index)
    t2.index.name = 'datetime'
    t2 = t2.iloc[:, :].filter(["WD_aws"]).astype(float)
    t2.columns = [vr]

    return [c, e, l, t, t1, t2]


def read_tcc():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    vr = 'tcc'

    # CARRA
    fn = 'thaao_carra_total_cloud_cover_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [vr]

    # ERA5
    fn = 'thaao_era5_total_cloud_cover_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e[2] = e.values * 100.
    e.columns = [vr]

    # THAAO (ceilometer)
    fn = '_Thule_CHM190147_000_0060cloud'
    for i in ceilometer_daterange:
        i_fmt = i.strftime('%Y%m%d')
        try:
            t_tmp = pd.read_table(
                    os.path.join(basefol_t_elab, 'thaao_ceilometer_elab', 'medie_tat_rianalisi', f'{i_fmt}{fn}.txt'),
                    skipfooter=0, sep='\s+', header=0, skiprows=9, engine='python')
            t_tmp[t_tmp == -9999.9] = np.nan
            t_tmp = pd.concat([t, t_tmp], axis=0)
            print(f'OK: {i_fmt}{fn}.txt')
        except (FileNotFoundError, pd.errors.EmptyDataError):
            print(f'NOT FOUND: {i_fmt}{fn}.txt')
    t.index = pd.to_datetime(t['#'] + ' ' + t['date[y-m-d]time[h:m:s]'], format='%Y-%m-%d %H:%M:%S')
    t.index.name = 'datetime'
    t = t.iloc[:, :].filter(['TCC[okt]']).astype(float)
    t.columns = [vr]

    return [c, e, l, t]


def read_cbh():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    vr = 'cbh'

    # CARRA
    fn = 'thaao_carra_cloud_base_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [vr]

    # ERA5
    fn = 'thaao_era5_cloud_base_height_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e.columns = [vr]

    # THAAO (ceilometer)
    fn = '_Thule_CHM190147_000_0060cloud'
    for i in ceilometer_daterange:
        i_fmt = i.strftime('%Y%m%d')
        try:
            t_tmp = pd.read_table(
                    os.path.join(basefol_t_elab, 'thaao_ceilometer_elab', 'medie_tat_rianalisi', f'{i_fmt}{fn}.txt'),
                    skipfooter=0, sep='\s+', header=0, skiprows=9, engine='python')
            t_tmp[t_tmp == -9999.9] = np.nan
            t = pd.concat([t, t_tmp], axis=0)
            print(f'OK: {i_fmt}{fn}.txt')
        except (FileNotFoundError, pd.errors.EmptyDataError):
            print(f'NOT FOUND: {i_fmt}{fn}.txt')
    t.index = pd.to_datetime(t['#'] + ' ' + t['date[y-m-d]time[h:m:s]'], format='%Y-%m-%d %H:%M:%S')
    t.index.name = 'datetime'
    t = t.iloc[:, :].filter(['CBH_L1[m]']).astype(float)
    t.columns = [vr]

    return [c, e, l, t]


def read_precip():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    t2 = pd.DataFrame()
    vr = 'precip'

    # CARRA
    fn = 'thaao_carra_total_precipitation_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c.columns = [vr]

    # ERA5
    fn = 'thaao_era5_total_precipitation_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e[2] = e.values * 1000.
    e.columns = [vr]

    # AWS ECAPAC
    fn = 'AWS_THAAO_'
    for i in aws_ecapac_daterange:
        i_fmt = i.strftime('%Y_%m_%d')
        try:
            file = os.path.join(
                    basefol_t, 'thaao_ecapac_aws_snow', 'AWS_ECAPAC', f'{fn}{i_fmt}_00_00.dat')
            t2_tmp = pd.read_csv(
                    file, skiprows=[0, 3], header=0, decimal='.', delimiter=',', engine='python',
                    index_col='TIMESTAMP').iloc[1:, :]
            t2 = pd.concat([t2, t2_tmp], axis=0)
            print(f'OK: {fn}{i_fmt}_00_00.dat')
        except (FileNotFoundError, pd.errors.EmptyDataError):
            print(f'NOT_FOUND: {fn}{i_fmt}_00_00.dat')
    t2.index = pd.DatetimeIndex(t2.index)
    t2.index.name = 'datetime'
    t2 = t2.iloc[:, :].filter(["PR"]).astype(float)
    t2.columns = [vr]
    t2 = t2

    return [c, e, l, t, t1, t2]


def read_lwp():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    t1 = pd.DataFrame()
    vr = 'lwp'

    # CARRA
    fn = 'thaao_carra_total_column_cloud_liquid_water_'
    for yy, year in enumerate(years):
        try:
            c_tmp = pd.read_table(
                    os.path.join(basefol_c, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            c = pd.concat([c, c_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    c.drop(columns=[0, 1], inplace=True)
    c[2] = c.values * 1000000
    c.columns = [vr]
    c[c < 0] = np.nan
    c[c < 15] = 0

    # ERA5
    fn = 'thaao_era5_total_column_cloud_liquid_water_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e[2] = e.values * 1000
    e.columns = [vr]
    e[c < 0] = np.nan
    e[e < 15] = 0

    # THAAO (hatpro)
    fn = 'LWP_15_min_'
    for yy, year in enumerate(years):
        try:
            t1_tmp = pd.read_table(
                    os.path.join(
                            basefol_t, 'thaao_hatpro', 'definitivi_da_giando', f'{fn}{year}_SITO',
                            f'{fn}{year}_SITO.dat'), sep='\s+', engine='python')
            tmp = np.empty(t1_tmp['JD_rif'].shape, dtype=dt.datetime)
            for ii, el in enumerate(t1_tmp['JD_rif']):
                new_jd_ass = el + julian.to_jd(dt.datetime(year - 1, 12, 31, 0, 0), fmt='jd')
                tmp[ii] = julian.from_jd(new_jd_ass, fmt='jd')
                tmp[ii] = tmp[ii].replace(microsecond=0)
            t1_tmp.index = pd.DatetimeIndex(tmp)
            t1_tmp['LWP_gm-2'] = t1_tmp['LWP_gm-2'].values
            t1_tmp.drop(columns=['JD_rif', 'RF', 'N', 'STD_LWP'], axis=1, inplace=True)
            t1 = pd.concat([t1, t1_tmp], axis=0)

            print(f'OK: {fn}{year}.dat')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.dat')
    t1.columns = [vr]
    # cleaning HATPRO DATA
    t1[t1 < 0] = np.nan
    t1[t1 < 15] = 0

    return [c, e, l, t, t1]


def read_lw_down():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    vr = 'lw_down'

    # # CARRA
    # fn = 'thaao_carra_thermal_surface_radiation_downwards_'

    # ERA5
    fn = 'thaao_era5_surface_thermal_radiation_downwards_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e[2] = e.values / 3600.  # originele in J*m-2
    e.columns = [vr]

    # THAAO
    fn = 'MERGED_SW_LW_UP_DW_METEO_'
    for yy, year in enumerate(years):
        try:
            t_tmp = pd.read_table(
                    os.path.join(basefol_t, 'thaao_rad', f'{fn}{year}_5MIN.dat'), engine='python', skiprows=None,
                    header=0, decimal='.', sep='\s+')
            tmp = np.empty(t_tmp['JDAY_UT'].shape, dtype=dt.datetime)
            for ii, el in enumerate(t_tmp['JDAY_UT']):
                new_jd_ass = el + julian.to_jd(dt.datetime(year - 1, 12, 31, 0, 0), fmt='jd')
                tmp[ii] = julian.from_jd(new_jd_ass, fmt='jd')
                tmp[ii] = tmp[ii].replace(microsecond=0)
            t_tmp.index = pd.DatetimeIndex(tmp)
            t_tmp.drop(
                    ['JDAY_UT', 'JDAY_LOC', 'SZA', 'SW_DOWN', 'PAR_DOWN', 'PAR_UP', 'SW_UP', 'LW_UP', 'TBP',
                     'ALBEDO_SW', 'ALBEDO_LW', 'ALBEDO_PAR', 'P', 'T', 'RH', 'PE', 'RR2'], axis=1, inplace=True)
            t = pd.concat([t, t_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    t.columns = [vr]

    return [c, e, l, t]


def read_lw_up():
    c = pd.DataFrame()
    e_n = pd.DataFrame()
    e_d = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    vr = 'lw_up'

    # # CARRA
    # fn1 = 'thaao_carra_surface_net_solar_radiation_'
    # fn2 = 'thaao_carra_surface_solarl_radiation_downwards_'

    # ERA5
    fn1 = 'thaao_era5_surface_net_thermal_radiation_'
    fn2 = 'thaao_era5_surface_thermal_radiation_downwards_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn1}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e_n = pd.concat([e_n, e_tmp], axis=0)
            print(f'OK: {fn1}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn1}{year}.txt')
    e_n.index = pd.to_datetime(e_n[0] + ' ' + e_n[1], format='%Y-%m-%d %H:%M:%S')
    e_n.drop(columns=[0, 1], inplace=True)
    e_n[2] = e_n.values / 3600.  # originele in J*m-2
    e_n.columns = ['surface_net_thermal_radiation']

    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn2}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e_d = pd.concat([e_d, e_tmp], axis=0)
            print(f'OK: {fn2}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn2}{year}.txt')
    e_d.index = pd.to_datetime(e_d[0] + ' ' + e_d[1], format='%Y-%m-%d %H:%M:%S')
    e_d.drop(columns=[0, 1], inplace=True)
    e_d[2] = e_d.values / 3600.  # originele in J*m-2
    e_d.columns = ['surface_thermal_radiation_downwards']

    e = pd.concat([e_n, e_d], axis=1)

    e['surface_thermal_radiation_upwards'] = e['surface_thermal_radiation_downwards'] - e[
        'surface_net_thermal_radiation']
    e.drop(columns=['surface_net_thermal_radiation', 'surface_thermal_radiation_downwards'], inplace=True)
    e.columns = [vr]

    # THAAO
    fn = 'MERGED_SW_LW_UP_DW_METEO_'
    for yy, year in enumerate(years):
        try:
            t_tmp = pd.read_table(
                    os.path.join(basefol_t, 'thaao_rad', f'{fn}{year}_5MIN.dat'), engine='python', skiprows=None,
                    header=0, decimal='.', sep='\s+')
            tmp = np.empty(t_tmp['JDAY_UT'].shape, dtype=dt.datetime)
            for ii, el in enumerate(t_tmp['JDAY_UT']):
                new_jd_ass = el + julian.to_jd(dt.datetime(year - 1, 12, 31, 0, 0), fmt='jd')
                tmp[ii] = julian.from_jd(new_jd_ass, fmt='jd')
                tmp[ii] = tmp[ii].replace(microsecond=0)
            t_tmp.index = pd.DatetimeIndex(tmp)
            t_tmp.drop(
                    ['JDAY_UT', 'JDAY_LOC', 'SZA', 'SW_DOWN', 'PAR_DOWN', 'PAR_UP', 'LW_DOWN', 'SW_UP', 'TBP',
                     'ALBEDO_SW', 'ALBEDO_LW', 'ALBEDO_PAR', 'P', 'T', 'RH', 'PE', 'RR2'], axis=1, inplace=True)
            t = pd.concat([t, t_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    t.columns = [vr]

    return [c, e, l, t]


def read_sw_down():
    c = pd.DataFrame()
    e = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    vr = 'sw_down'

    # # CARRA
    # fn = 'thaao_carra_surface_solar_radiation_downwards_'
    # for yy, year in enumerate(years):
    #     try:
    #         c_tmp = pd.read_table(
    #                 os.path.join(basefol_c, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
    #                 engine='python')
    #         c = pd.concat([c, c_tmp], axis=0)
    #         print(f'OK: {fn}{year}.txt')
    #     except FileNotFoundError:
    #         print(f'NOT FOUND: {fn}{year}.txt')
    # c.index = pd.to_datetime(c[0] + ' ' + c[1], format='%Y-%m-%d %H:%M:%S')
    # c.drop(columns=[0, 1], inplace=True)
    # c[2] = c.values / 10800.
    # c.columns = [vr]

    # ERA5
    fn = 'thaao_era5_surface_solar_radiation_downwards_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e = pd.concat([e, e_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    e.index = pd.to_datetime(e[0] + ' ' + e[1], format='%Y-%m-%d %H:%M:%S')
    e.drop(columns=[0, 1], inplace=True)
    e[2] = e.values / 3600.  # originale in J*m-2
    e.columns = [vr]

    # THAAO
    fn = 'MERGED_SW_LW_UP_DW_METEO_'
    for yy, year in enumerate(years):
        try:
            t_tmp = pd.read_table(
                    os.path.join(basefol_t, 'thaao_rad', f'{fn}{year}_5MIN.dat'), engine='python', skiprows=None,
                    header=0, decimal='.', sep='\s+')
            tmp = np.empty(t_tmp['JDAY_UT'].shape, dtype=dt.datetime)
            for ii, el in enumerate(t_tmp['JDAY_UT']):
                new_jd_ass = el + julian.to_jd(dt.datetime(year - 1, 12, 31, 0, 0), fmt='jd')
                tmp[ii] = julian.from_jd(new_jd_ass, fmt='jd')
                tmp[ii] = tmp[ii].replace(microsecond=0)
            t_tmp.index = pd.DatetimeIndex(tmp)
            t_tmp.drop(
                    ['JDAY_UT', 'JDAY_LOC', 'SZA', 'SW_UP', 'PAR_DOWN', 'PAR_UP', 'LW_DOWN', 'LW_UP', 'TBP',
                     'ALBEDO_SW', 'ALBEDO_LW', 'ALBEDO_PAR', 'P', 'T', 'RH', 'PE', 'RR2'], axis=1, inplace=True)
            t = pd.concat([t, t_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    t.columns = [vr]

    return [c, e, l, t]


def read_sw_up():
    c = pd.DataFrame()
    e_n = pd.DataFrame()
    e_d = pd.DataFrame()
    l = pd.DataFrame()
    t = pd.DataFrame()
    vr = 'sw_up'

    # # CARRA
    # fn1 = 'thaao_carra_surface_net_solar_radiation_'
    # fn2 = 'thaao_carra_surface_solar_radiation_downwards_'

    # ERA5
    fn1 = 'thaao_era5_surface_net_solar_radiation_'
    fn2 = 'thaao_era5_surface_solar_radiation_downwards_'
    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn1}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e_n = pd.concat([e_n, e_tmp], axis=0)
            print(f'OK: {fn1}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn1}{year}.txt')
    e_n.index = pd.to_datetime(e_n[0] + ' ' + e_n[1], format='%Y-%m-%d %H:%M:%S')
    e_n.drop(columns=[0, 1], inplace=True)
    e_n[2] = e_n.values / 3600.  # originele in J*m-2
    e_n.columns = ['surface_net_solar_radiation']

    for yy, year in enumerate(years):
        try:
            e_tmp = pd.read_table(
                    os.path.join(basefol_e, f'{fn2}{year}.txt'), skipfooter=1, sep='\s+', header=None, skiprows=1,
                    engine='python')
            e_tmp[e_tmp == -32767.0] = np.nan
            e_d = pd.concat([e_d, e_tmp], axis=0)
            print(f'OK: {fn2}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn2}{year}.txt')
    e_d.index = pd.to_datetime(e_d[0] + ' ' + e_d[1], format='%Y-%m-%d %H:%M:%S')
    e_d.drop(columns=[0, 1], inplace=True)
    e_d[2] = e_d.values / 3600.  # originele in J*m-2
    e_d.columns = ['surface_solar_radiation_downwards']

    e = pd.concat([e_n, e_d], axis=1)

    e['surface_solar_radiation_upwards'] = e['surface_solar_radiation_downwards'] - e['surface_net_solar_radiation']
    e.drop(columns=['surface_net_solar_radiation', 'surface_solar_radiation_downwards'], inplace=True)
    e.columns = [vr]

    # THAAO
    fn = 'MERGED_SW_LW_UP_DW_METEO_'
    for yy, year in enumerate(years):
        try:
            t_tmp = pd.read_table(
                    os.path.join(basefol_t, 'thaao_rad', f'{fn}{year}_5MIN.dat'), engine='python', skiprows=None,
                    header=0, decimal='.', sep='\s+')
            tmp = np.empty(t_tmp['JDAY_UT'].shape, dtype=dt.datetime)
            for ii, el in enumerate(t_tmp['JDAY_UT']):
                new_jd_ass = el + julian.to_jd(dt.datetime(year - 1, 12, 31, 0, 0), fmt='jd')
                tmp[ii] = julian.from_jd(new_jd_ass, fmt='jd')
                tmp[ii] = tmp[ii].replace(microsecond=0)
            t_tmp.index = pd.DatetimeIndex(tmp)
            t_tmp.drop(
                    ['JDAY_UT', 'JDAY_LOC', 'SZA', 'SW_DOWN', 'PAR_DOWN', 'PAR_UP', 'LW_DOWN', 'LW_UP', 'TBP',
                     'ALBEDO_SW', 'ALBEDO_LW', 'ALBEDO_PAR', 'P', 'T', 'RH', 'PE', 'RR2'], axis=1, inplace=True)
            # 'JDAY_UT', 'JDAY_LOC', 'SZA', 'SW_DOWN', 'SW_UP', 'PAR_DOWN', 'PAR_UP', 'LW_DOWN', 'LW_UP', 'TBP', 'ALBEDO_SW', 'ALBEDO_LW', 'ALBEDO_PAR', 'P', 'T', 'RH', 'PE', 'RR2'
            t = pd.concat([t, t_tmp], axis=0)
            print(f'OK: {fn}{year}.txt')
        except FileNotFoundError:
            print(f'NOT FOUND: {fn}{year}.txt')
    t.columns = [vr]

    return [c, e, l, t]


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
    if var == 'lwp':
        return read_lwp()
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
    if var == 'lw_down':
        return read_lw_down()
    if var == 'lw_up':
        return read_lw_up()
    if var == 'sw_down':
        return read_sw_down()
    if var == 'sw_up':
        return read_sw_up()
