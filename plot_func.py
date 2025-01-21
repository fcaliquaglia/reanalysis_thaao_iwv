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
uom = ' [kg/m2]'

def plot_ts(avar, period_label):
    """

    :param avar:
    :param period_label:
    :return:
    """
    print('TIMESERIES')
    [vr_c, vr_e, vr_l, vr_t, vr_t1, vr_t2, vr_c_res, vr_e_res, vr_l_res, vr_t_res, vr_t1_res, vr_t2_res] = avar
    fig, ax = plt.subplots(len(years), 1, figsize=(12, 17), dpi=300)
    fig.suptitle(f'{var_name_u} all {tres}', fontweight='bold')

    label_t2_ori = 'RS'
    label_t2 = 'RS'

    for [yy, year] in enumerate(years):
        print(f'plotting {year}')

        # original resolution
        try:
            ax[yy].plot(
                    vr_c[vr_c.index.year == year], color=c_col_ori, label=var_name_u + ' CARRA 3h', alpha=0.2, lw=0, marker='.', ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_e[vr_e.index.year == year], color=e_col_ori, label=var_name_u + ' ERA5 1h', alpha=0.2, lw=0, marker='.', ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_l[vr_l.index.year == year], color=l_col_ori, label=var_name_u + ' ERA5 1h', alpha=0.2, lw=0, marker='.', ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_t[vr_t.index.year == year], color=t_col_ori, label=var_name_u + ' THAAO ori', alpha=0.2, lw=0, marker='.',
                    ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_t1[vr_t1.index.year == year], color=t1_col_ori, label=var_name_u + ' HATPRO ori', alpha=0.2, lw=0, marker='.',
                    ms=1)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_t2[vr_t2.index.year == year], color=t2_col_ori, label=var_name_u + ' ' + label_t2_ori, alpha=0.02, lw=0, marker='.',
                    ms=1)
        except AttributeError:
            pass

        # resampled resolution
        try:
            ax[yy].plot(
                    vr_c_res[vr_c_res.index.year == year], color=c_col, label=var_name_u + ' CARRA', lw=0, marker='.', ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_e_res[vr_e_res.index.year == year], color=e_col, label=var_name_u + ' ERA5', lw=0, marker='.', ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_l_res[vr_l_res.index.year == year], color=l_col, label=var_name_u + ' ERA5-L', lw=0, marker='.', ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_t_res[vr_t_res.index.year == year], color=t_col, label=var_name_u + ' THAAO', lw=0, marker='.', ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_t1_res[vr_t1_res.index.year == year], color=t1_col, label=var_name_u + ' HATPRO' , lw=0, marker='.', ms=2)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    vr_t2_res[vr_t2_res.index.year == year], color=t2_col, label=var_name_u + ' ' + label_t2, lw=0, marker='.', ms=2)
        except AttributeError:
            pass

        ax[yy].set_ylim(extr[var_name]['min'], extr[var_name]['max'])
        ax[yy].text(0.45, 0.85, year, transform=ax[yy].transAxes)
        ax[yy].set_xticklabels([])
        ax[yy].set_xlim(dt.datetime(year, 1, 1), dt.datetime(year, 12, 31))
        ax[yy].text(0.01, 0.90, letters[yy] + ')', transform=ax[yy].transAxes)
    ax[-1].xaxis.set_major_formatter(myFmt)
    ax[-1].set_xlabel('Time')
    plt.legend(ncol=2)
    plt.savefig(os.path.join(basefol_out, tres, f'{tres}_{period_label}_{var_name}_only.png'))
    plt.close('all')


def plot_residuals(avar, period_label):
    """

    :param avar:
    :param period_label:
    :return:
    """
    print('RESIDUALS')
    [vr_c, vr_e, vr_l, vr_t, vr_t1, vr_t2, vr_c_res, vr_e_res, vr_l_res, vr_t_res, vr_t1_res, vr_t2_res] = avar
    fig, ax = plt.subplots(len(years), 1, figsize=(12, 17), dpi=300)
    fig.suptitle(f'residuals {var_name_u} all {tres}', fontweight='bold')
    for [yy, year] in enumerate(years):
        print(f'plotting {year}')

        label_t2 = var_name_u + ' RS'

        vr_ref = vr_t_res.resample(tres).mean()
        # resampled resolution
        daterange = pd.date_range(dt.datetime(year, 1, 1), dt.datetime(year, 12, 31))
        ax[yy].plot(daterange, np.repeat(0, len(daterange)), color='black', lw=2, ls='--')
        try:
            ax[yy].plot(
                    (vr_c_res[vr_c_res.index.year == year] - vr_ref[vr_ref.index.year == year]), color=c_col,
                    label=var_name_u + ' CARRA', lw=1, marker='.', ms=0)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    (vr_e_res[vr_e_res.index.year == year] - vr_ref[vr_ref.index.year == year]), color=e_col,
                    label=var_name_u + ' ERA5', lw=1, marker='.', ms=0)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    (vr_l_res[vr_l_res.index.year == year] - vr_ref[vr_ref.index.year == year]), color=l_col,
                    label=var_name_u + ' ERA5-L', lw=1, marker='.', ms=0)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    (vr_t1_res[vr_t1_res.index.year == year] - vr_ref[vr_ref.index.year == year]), color=t1_col,
                    label=var_name_u + ' HATPRO', lw=1, marker='.', ms=0)
        except AttributeError:
            pass
        try:
            ax[yy].plot(
                    (vr_t2_res[vr_t2_res.index.year == year] - vr_ref[vr_ref.index.year == year]), color=t2_col,
                    label=var_name_u + ' ' + label_t2, lw=1, marker='.', ms=0)
        except AttributeError:
            pass

        ax[yy].set_ylim(extr[var_name]['res_min'], extr[var_name]['res_max'])
        ax[yy].text(0.45, 0.85, year, transform=ax[yy].transAxes)
        ax[yy].set_xlim(dt.datetime(year, 1, 1), dt.datetime(year, 12, 31))
        # panel letters
        ax[yy].set_xticklabels([])
        ax[yy].text(0.01, 0.90, letters[yy] + ')', transform=ax[yy].transAxes)
    ax[-1].xaxis.set_major_formatter(myFmt)
    ax[-1].set_xlabel('Time')
    plt.legend(ncol=2)
    # plt.tight_layout()
    plt.savefig(os.path.join(basefol_out, tres, f'{tres}_{period_label}_residuals_{var_name}_only.png'))
    plt.close('all')


def plot_scatter(avar, period_label):
    """

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
    x = vr_t_res[var_name].resample(tres).mean()
    xlabel = 'VESPA'

    for i, comp in enumerate(comps):
        axs[i].set_xlabel(xlabel)
        if comp == 'c':
            label = 'CARRA'
            axs[i].set_ylabel(label)
            try:
                y = vr_c_res[var_name]
            except KeyError:
                print(f'error with {label}')
                continue
        if comp == 'e':
            label = 'ERA5'
            axs[i].set_ylabel(label)
            try:
                y = vr_e_res[var_name]
            except KeyError:
                print(f'error with {label}')
                continue
        if comp == 'l':
            label = 'ERA5-L'
            axs[i].set_ylabel(label)
            try:
                y = vr_l_res[var_name]
            except KeyError:
                print(f'error with {label}')
                continue
        if comp == 't':
            label = 'THAAO'
            axs[i].set_ylabel(label)
            try:
                y = vr_t_res[var_name]
            except KeyError:
                print(f'error with {label}')
                continue
        if comp == 't1':
            label = 'HATPRO'
            axs[i].set_ylabel(label)
            try:
                y = vr_t1_res[var_name]
            except KeyError:
                print(f'error with {label}')
                continue
        if comp == 't2':
            label = 'RS'
            axs[i].set_ylabel(label)
            try:
                y = vr_t2_res[var_name].dropna()
            except KeyError:
                print(f'error with {label}')
                continue
        try:
            print(f'plotting scatter THAAO-{label}')

            fig.suptitle(f'{var_name_u} {seas_name} {tres}', fontweight='bold')
            axs[i].set_title(label)

            time_list = pd.date_range(start=dt.datetime(years[0], 1, 1), end=dt.datetime(years[-1], 12, 31), freq=tres)

            x_all = x.reindex(time_list).fillna(np.nan)
            x_s = x_all.loc[(x_all.index.month.isin(seass[period_label]['months']))]
            y_all = y.reindex(time_list).fillna(np.nan)
            y_s = y_all.loc[(y_all.index.month.isin(seass[period_label]['months']))]
            idx = ~(np.isnan(x_s) | np.isnan(y_s))

            if seas_name != 'all':
                if label == 'RS':
                    y_s = y.loc[(y.index.month.isin(seass[period_label]['months']))]
                    x_s = pd.Series(vr_t_res.reindex(y_s.index)['iwv'])

                    idx = ~(np.isnan(x_s) | np.isnan(y_s))
                    axs[i].scatter(
                            x_s[idx].values, y_s[idx].values, s=50, facecolor='none', color=seass[period_label]['col'],
                            label=period_label)
                else:
                    axs[i].scatter(
                            x_s[idx], y_s[idx], s=5, color=seass[period_label]['col'], facecolor='none', alpha=0.5,
                            label=period_label)
            else:
                if label == 'RS':
                    y_s = y.loc[(y.index.month.isin(seass[period_label]['months']))]
                    x_s = pd.Series(vr_t_res.reindex(y_s.index)['iwv'])

                    idx = ~(np.isnan(x_s) | np.isnan(y_s))
                    axs[i].scatter(x_s[idx], y_s[idx], facecolor='none', s=50, color=seass[period_label]['col'])
                else:
                    bin_nr= 200
                    bin_size=extr[var_name]['max']/bin_nr
                    h = axs[i].hist2d(x_s[idx], y_s[idx], bins=bin_nr, cmap=plt.cm.jet, cmin=1, vmin=1)
                    axs[i].text(0.30, 0.60, f'bin_size={bin_size} kg/m3')
                    #fig.colorbar(h[3], ax=axs[i], extend='both')

            if len(x_s[idx]) < 2 | len(y_s[idx]) < 2:
                print('EEEEEEEEEEEEEEEEEEE')
            else:
                b, a = np.polyfit(x_s[idx], y_s[idx], deg=1)
                xseq = np.linspace(extr[var_name]['min'], extr[var_name]['max'], num=1000)
                axs[i].plot(xseq, a + b * xseq, color='red', lw=2.5, ls='--', alpha=0.5)
                axs[i].plot(
                        [extr[var_name]['min'], extr[var_name]['max']], [extr[var_name]['min'], extr[var_name]['max']], color='black', lw=1.5,
                        ls='-')
                corcoef = ma.corrcoef(x_s[idx], y_s[idx])

                N = len(y_s[idx])
                rmse = np.sqrt(np.nanmean((x_s[idx] - y_s[idx]) ** 2))
                mae = np.nanmean(np.abs(x_s[idx] - y_s[idx]))
                axs[i].text(
                        0.60, 0.15, f'R={corcoef[0, 1]:1.3}\nrmse={rmse:1.3}\nN={N}\nmae={mae:1.3}', fontsize=14,
                        transform=axs[i].transAxes)
                axs[i].set_xlim(extr[var_name]['min'], extr[var_name]['max'])
                axs[i].set_ylim(extr[var_name]['min'], extr[var_name]['max'])
                axs[i].text(0.05, 0.95, letters[i] + ')', transform=axs[i].transAxes)

        except:
            print(f'error with {label}')

    plt.savefig(os.path.join(basefol_out, tres, f'{tres}_scatter_{seas_name}_{var_name}_only.png'))
    plt.close('all')


def plot_scatter_cum(avar):
    """

    :param avar:
    :return:
    """
    import copy as cp
    fig, ax = plt.subplots(2, 2, figsize=(12, 12), dpi=300)
    seass_new = cp.copy(seass)
    seass_new.pop('all')
    for period_label in seass_new:
        print('SCATTERPLOTS')
        [vr_c, vr_e, vr_l, vr_t, vr_t1, vr_t2, vr_c_res, vr_e_res, vr_l_res, vr_t_res, vr_t1_res, vr_t2_res] = avar
        seas_name = seass[period_label]['name']
        axs = ax.ravel()

        comps = ['c', 'e', 't1', 't2']
        x = vr_t_res[vr].resample(tres).mean()
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
                    y = vr_t2_res[vr].dropna()
                except KeyError:
                    print(f'error with {label}')
                    continue
            try:
                print(f'plotting scatter THAAO-{label}')

                fig.suptitle(f'{var_name_u} cumulative plot', fontweight='bold')
                axs[i].set_title(label)

                time_list = pd.date_range(
                        start=dt.datetime(years[0], 1, 1), end=dt.datetime(years[-1], 12, 31), freq=tres)

                x_all = x.reindex(time_list).fillna(np.nan)
                x_s = x_all.loc[(x_all.index.month.isin(seass[period_label]['months']))]
                y_all = y.reindex(time_list).fillna(np.nan)
                y_s = y_all.loc[(y_all.index.month.isin(seass[period_label]['months']))]
                idx = ~(np.isnan(x_s) | np.isnan(y_s))

                if seas_name != 'all':
                    if label == 'RS':
                        y_s = y.loc[(y.index.month.isin(seass[period_label]['months']))]
                        x_s = pd.Series(vr_t_res.reindex(y_s.index)['iwv'])

                        idx = ~(np.isnan(x_s) | np.isnan(y_s))
                        axs[i].scatter(
                                x_s[idx].values, y_s[idx].values, s=50, facecolor='none',
                                color=seass[period_label]['col'], label=period_label)

                    else:
                        axs[i].scatter(
                                x_s[idx], y_s[idx], s=5, color=seass[period_label]['col'], edgecolors='none', alpha=0.5,
                                label=period_label)

                if len(x_s[idx]) < 2 | len(y_s[idx]) < 2:
                    print('ERROR, ERROR, NO DATA ENOUGH FOR PROPER FIT (i.e. only 1 point available)')
                else:
                    b, a = np.polyfit(x_s[idx], y_s[idx], deg=1)
                    xseq = np.linspace(extr[vr]['min'], extr[vr]['max'], num=1000)
                    axs[i].plot(xseq, a + b * xseq, color=seass[period_label]['col'], lw=2.5, ls='--', alpha=0.5)
                    axs[i].plot(
                            [extr[vr]['min'], extr[vr]['max']], [extr[vr]['min'], extr[vr]['max']], color='black',
                            lw=1.5, ls='-')

                    axs[i].set_xlim(extr[vr]['min'], extr[vr]['max'])
                    axs[i].set_ylim(extr[vr]['min'], extr[vr]['max'])
                    axs[i].text(0.05, 0.95, letters[i] + ')', transform=axs[i].transAxes)
                    axs[i].legend()
            except:
                print(f'error with {label}')

    plt.savefig(os.path.join(basefol_out, tres, f'{tres}_scatter_cum_{vr}_only.png'))
    plt.close('all')
