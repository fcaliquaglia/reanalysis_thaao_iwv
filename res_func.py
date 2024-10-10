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

import pandas as pd
from metpy.calc import wind_components, wind_direction
from metpy.units import units

from read_func import read


def precip_res(pr, timeres):
    pr_res = pr.resample(timeres).sum()

    return pr_res


def lwp_res(lwp, timeres):
    lwp_res = lwp.resample(timeres).sum()

    return lwp_res

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


def data_resampling(vr, tres, var_c, var_e, var_l, var_t, var_t1, var_t2):
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
    elif vr == 'lwp':
        try:
            var_c_res = lwp_res(var_c, tres)
        except (TypeError, NameError, ValueError):
            var_c_res = pd.DataFrame()
        try:
            var_e_res = lwp_res(var_e, tres)
        except (TypeError, NameError, ValueError):
            var_e_res = pd.DataFrame()
        try:
            var_l_res = lwp_res(var_l, tres)
        except:
            var_l_res = pd.DataFrame()
        try:
            var_t_res = lwp_res(var_t, tres)
        except (TypeError, NameError, ValueError):
            var_t_res = pd.DataFrame()
        try:
            var_t1_res = lwp_res(var_t1, tres)
        except (TypeError, NameError, ValueError):
            var_t1_res = pd.DataFrame()
        try:
            var_t2_res = lwp_res(var_t2, tres)
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
