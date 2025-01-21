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

from plot_func import plot_residuals, plot_scatter, plot_scatter_cum, plot_ts
from read_func import *
from res_func import data_resampling

if __name__ == "__main__":

    var = 'iwv'
    print(var)
    [var_c, var_e, var_l, var_t, var_t1, var_t2] = read_iwv()

    # time RESAMPLING
    var_c_res, var_e_res, var_l_res, var_t_res, var_t1_res, var_t2_res = data_resampling(
            var, var_c, var_e, var_l, var_t, var_t1, var_t2)

    all_var = [var_c, var_e, var_l, var_t, var_t1, var_t2, var_c_res, var_e_res, var_l_res, var_t_res, var_t1_res,
               var_t2_res]
    plot_ts(var, all_var, 'all')
    plot_residuals(var, all_var, 'all')
    plot_scatter_cum(var, all_var)
    for seas in seass:
        plot_scatter(var, all_var, seas)
        # plot_ba(var, all_var, seas)
