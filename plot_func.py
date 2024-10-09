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

import matplotlib.pyplot as plt

from inputs import *


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
