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

import datetime as dt
import string

import matplotlib.pyplot as plt
import numpy.ma as ma
import pandas as pd

from inputs import *

letters = list(string.ascii_lowercase)


def plot_ts(vr, avar, period_label):
    """

    :param vr:
    :param avar:
    :param period_label:
    :return:
    """
    print('TIMESERIES')
    [vr_c, vr_e, vr_l, vr_t, vr_t1, vr_t2, vr_c_res, vr_e_res, vr_l_res, vr_t_res, vr_t1_res, vr_t2_res] = avar
    fig, ax = plt.subplots(len(years), 1, figsize=(12, 17), dpi=300)
    fig.suptitle(f'{vr.upper()} all {tres}', fontweight='bold')

    label_t2_ori = 'RS'
    label_t2 = 'RS'

    for [yy, year] in enumerate(years):
        print(f'plotting {year}')

        # original resolution
        try:
            ax[yy].plot(
                    vr_c[vr_c.index.year == year], color=c_col_ori, label='CARRA 3h', alpha=0.2, lw=0, marker='.', ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_e[vr_e.index.year == year], color=e_col_ori, label='ERA5 1h', alpha=0.2, lw=0, marker='.', ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_l[vr_l.index.year == year], color=l_col_ori, label='ERA5 1h', alpha=0.2, lw=0, marker='.', ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_t[vr_t.index.year == year], color=t_col_ori, label='THAAO ori', alpha=0.2, lw=0, marker='.',
                    ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_t1[vr_t1.index.year == year], color=t1_col_ori, label='HATPRO ori', alpha=0.2, lw=0, marker='.',
                    ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_t2[vr_t2.index.year == year], color=t2_col_ori, label=label_t2_ori, alpha=0.02, lw=0, marker='.',
                    ms=1)
        except AttributeError:
            pass

        # resampled resolution
        try:
            ax[yy].plot(
                    vr_c_res[vr_c_res.index.year == year], color=c_col, label='CARRA', lw=0, marker='.', ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_e_res[vr_e_res.index.year == year], color=e_col, label='ERA5', lw=0, marker='.', ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_l_res[vr_l_res.index.year == year], color=l_col, label='ERA5-L', lw=0, marker='.', ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_t_res[vr_t_res.index.year == year], color=t_col, label='THAAO', lw=0, marker='.', ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_t1_res[vr_t1_res.index.year == year], color=t1_col, label='HATPRO', lw=0, marker='.', ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_t2_res[vr_t2_res.index.year == year], color=t2_col, label=label_t2, lw=0, marker='.', ms=2)
        except AttributeError:
            pass

        ax[yy].set_ylim(extr[vr]['min'], extr[vr]['max'])
        ax[yy].text(0.45, 0.85, year, transform=ax[yy].transAxes)
        ax[yy].set_xticklabels([])
        ax[yy].set_xlim(dt.datetime(year, 1, 1), dt.datetime(year, 12, 31))
        ax[yy].text(0.1, 0.8, letters[yy] + ')', transform=ax[yy].transAxes)
    ax[-1].xaxis.set_major_formatter(myFmt)
    ax[-1].set_xlabel('Time')
    plt.legend(ncol=2)
    plt.savefig(os.path.join(basefol_out, tres, f'{tres}_{period_label}_{vr}_only.png'))
    plt.close('all')


def plot_residuals(vr, avar, period_label):
    """

    :param vr:
    :param avar:
    :param period_label:
    :return:
    """
    print('RESIDUALS')
    [vr_c, vr_e, vr_l, vr_t, vr_t1, vr_t2, vr_c_res, vr_e_res, vr_l_res, vr_t_res, vr_t1_res, vr_t2_res] = avar
    fig, ax = plt.subplots(len(years), 1, figsize=(12, 17), dpi=300)
    fig.suptitle(f'residuals {vr.upper()} all {tres}', fontweight='bold')
    for [yy, year] in enumerate(years):
        print(f'plotting {year}')

        if vr != 'iwv':
            label_t2 = 'AWS ECAPAC'
        else:
            label_t2 = 'RS'

        if vr == 'lwp':
            vr_ref = vr_t1_res
        elif vr in ['precip', 'windd', 'winds']:
            vr_ref = vr_t2_res
        elif vr in ['iwv']:
            vr_ref = vr_t_res.resample(tres).mean()
        else:
            vr_ref = vr_t_res
        # resampled resolution
        daterange = pd.date_range(dt.datetime(year, 1, 1), dt.datetime(year, 12, 31))
        ax[yy].plot(daterange, np.repeat(0, len(daterange)), color='black', lw=2, ls='--')
        try:
            ax[yy].plot(
                    (vr_c_res[vr_c_res.index.year == year] - vr_ref[vr_ref.index.year == year]), color=c_col,
                    label='CARRA', lw=1, marker='.', ms=0)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    (vr_e_res[vr_e_res.index.year == year] - vr_ref[vr_ref.index.year == year]), color=e_col,
                    label='ERA5', lw=1, marker='.', ms=0)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    (vr_l_res[vr_l_res.index.year == year] - vr_ref[vr_ref.index.year == year]), color=l_col,
                    label='ERA5-L', lw=1, marker='.', ms=0)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    (vr_t1_res[vr_t1_res.index.year == year] - vr_ref[vr_ref.index.year == year]), color=t1_col,
                    label='HATPRO', lw=1, marker='.', ms=0)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    (vr_t2_res[vr_t2_res.index.year == year] - vr_ref[vr_ref.index.year == year]), color=t2_col,
                    label=label_t2, lw=1, marker='.', ms=0)
        except AttributeError:
            pass

        ax[yy].set_ylim(extr[vr]['res_min'], extr[vr]['res_max'])
        ax[yy].text(0.45, 0.85, year, transform=ax[yy].transAxes)
        ax[yy].set_xlim(dt.datetime(year, 1, 1), dt.datetime(year, 12, 31))
        # panel letters
        ax[yy].set_xticklabels([])
        ax[yy].text(0.1, 0.8, letters[yy] + ')', transform=ax[yy].transAxes)
    ax[-1].xaxis.set_major_formatter(myFmt)
    ax[-1].set_xlabel('Time')
    plt.legend(ncol=2)
    # plt.tight_layout()
    plt.savefig(os.path.join(basefol_out, tres, f'{tres}_{period_label}_residuals_{vr}_only.png'))
    plt.close('all')


def plot_scatter(vr, avar, period_label):
    """

    :param vr:
    :param avar:
    :param period_label:
    :return:
    """
    print('SCATTERPLOTS')
    [vr_c, vr_e, vr_l, vr_t, vr_t1, vr_t2, vr_c_res, vr_e_res, vr_l_res, vr_t_res, vr_t1_res, vr_t2_res] = avar
    seas_name = seass[period_label]['name']
    fig, ax = plt.subplots(2, 2, figsize=(12, 12), dpi=300)
    axs = ax.ravel()



    comps = ['c', 'e', 't1', 't2']
    x = vr_t_res[vr]
    xlabel = 'VESPA'

    for i, comp in enumerate(comps):
        axs[i].set_xlabel(xlabel)
        if comp == 'c':
            label = 'CARRA'
            axs[i].set_ylabel(label)
            try:
                y = vr_c_res[vr]
            except KeyError:
                print(f'error with {label}')
                continue
        if comp == 'e':
            label = 'ERA5'
            axs[i].set_ylabel(label)
            try:
                y = vr_e_res[vr]
            except KeyError:
                print(f'error with {label}')
                continue
        if comp == 'l':
            label = 'ERA5-L'
            axs[i].set_ylabel(label)
            try:
                y = vr_l_res[vr]
            except KeyError:
                print(f'error with {label}')
                continue
        if comp == 't':
            label = 'THAAO'
            axs[i].set_ylabel(label)
            try:
                y = vr_t_res[vr]
            except KeyError:
                print(f'error with {label}')
                continue
        if comp == 't1':
            label = 'HATPRO'
            axs[i].set_ylabel(label)
            try:
                y = vr_t1_res[vr]
            except KeyError:
                print(f'error with {label}')
                continue
        if comp == 't2':
            label = 'RS'
            axs[i].set_ylabel(label)
            try:
                y = vr_t2_res[vr]
            except KeyError:
                print(f'error with {label}')
                continue
        try:
            print(f'plotting scatter THAAO-{label}')

            fig.suptitle(f'{vr.upper()} {seas_name} {tres}', fontweight='bold')
            axs[i].set_title(label)

            time_list = pd.date_range(start=dt.datetime(years[0], 1, 1), end=dt.datetime(years[-1], 12, 31), freq=tres)
            time_list_rs = pd.date_range(
                start=dt.datetime(years[0], 1, 1), end=dt.datetime(years[-1], 12, 31), freq=tres_rs)

            x_all = x.reindex(time_list)
            x_s = x_all.loc[(x_all.index.month.isin(seass[period_label]['months']))]
            y_all = y.reindex(time_list)
            y_s = y_all.loc[(y_all.index.month.isin(seass[period_label]['months']))]
            idx = ~(np.isnan(x_s) | np.isnan(y_s))

            if seas_name != 'all':
                axs[i].scatter(x_s[idx], y_s[idx], color=seass[period_label]['col'])
            else:
                if label == 'RS':
                    x_all_rs = x.reindex(time_list_rs)
                    x_s = x_all_rs.loc[(x_all_rs.index.month.isin(seass[period_label]['months']))]
                    y_all_rs = y.reindex(time_list_rs)
                    y_s = y_all_rs.loc[(y_all_rs.index.month.isin(seass[period_label]['months']))]
                    idx_rs = ~(np.isnan(x_s) | np.isnan(y_s))
                    axs[i].scatter(x_s[idx_rs], y_s[idx_rs], color=seass[period_label]['col'])
                else:
                    h = axs[i].hist2d(x_s[idx], y_s[idx], bins=(250, 250), cmap=plt.cm.jet, cmin=1, vmin=1)
                    fig.colorbar(h[3], ax=axs[i], extend='both')

            b, a = np.polyfit(x_s[idx], y_s[idx], deg=1)
            xseq = np.linspace(extr[vr]['min'], extr[vr]['max'], num=1000)
            axs[i].plot(xseq, a + b * xseq, color='red', lw=2.5, ls='--', alpha=0.5)
            axs[i].plot(
                    [extr[vr]['min'], extr[vr]['max']], [extr[vr]['min'], extr[vr]['max']], color='black', lw=1.5,
                    ls='-')
            corcoef = ma.corrcoef(x_s[idx], y_s[idx])

            N = x_s[idx].shape[0]
            rmse = np.sqrt(np.nanmean((x_s[idx] - y_s[idx]) ** 2))
            mae = np.nanmean(np.abs(x_s[idx] - y_s[idx]))
            axs[i].text(
                    0.60, 0.15, f'R={corcoef[0, 1]:1.3}\nrmse={rmse:1.3}\nN={N}\nmae={mae:1.3}', fontsize=14,
                    transform=axs[i].transAxes)
            axs[i].set_xlim(extr[vr]['min'], extr[vr]['max'])
            axs[i].set_ylim(extr[vr]['min'], extr[vr]['max'])
            axs[i].text(0.1, 0.8, letters[i] + ')', transform=axs[i].transAxes)
        except:
            print(f'error with {label}')

    plt.savefig(os.path.join(basefol_out, tres, f'{tres}_scatter_{seas_name}_{vr}_only.png'))
    plt.close('all')
