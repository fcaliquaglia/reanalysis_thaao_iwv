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
    kwargs_ori = {'alpha': 0.02, 'lw': 0, 'marker': '.', 'ms': 1}
    kwargs = {'lw': 0, 'marker': '.', 'ms': 2}

    for [yy, year] in enumerate(years):
        print(f'plotting {year}')

        # original resolution
        for (vr, vr_n) in zip([vr_c, vr_e, vr_l, vr_t, vr_t1, vr_t2], var_names):
            try:
                data = vr[vr.index.year == year]
                ax[yy].plot(data, color=var_dict[vr_n]['col_ori'], **kwargs_ori)
            except AttributeError:
                pass

        # resampled resolution
        for (vr, vr_n) in zip([vr_c_res, vr_e_res, vr_l_res, vr_t_res, vr_t1_res, vr_t2_res], var_names):
            try:
                data = vr[vr.index.year == year]
                ax[yy].plot(data, color=var_dict[vr_n]['col'], label=var_dict[vr_n]['label_uom'], **kwargs)
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
    kwargs = {'lw': 1, 'marker': '.', 'ms': 0}

    for [yy, year] in enumerate(years):
        print(f'plotting {year}')

        daterange = pd.date_range(dt.datetime(year, 1, 1), dt.datetime(year, 12, 31))
        ax[yy].plot(daterange, np.repeat(0, len(daterange)), color='black', lw=2, ls='--')

        # resampled resolution
        vr_ref = vr_t_res.resample(tres).mean()
        for (vr, vr_n) in zip([vr_c_res, vr_e_res, vr_l_res, vr_t_res, vr_t1_res, vr_t2_res], var_names):
            try:
                data = vr[vr.index.year == year]
                ax[yy].plot(
                        (data - vr_ref[vr_ref.index.year == year]), color=var_dict[vr_n]['col'],
                        label=var_dict[vr_n]['label_uom'], **kwargs)
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
    plt.savefig(os.path.join(basefol_out, tres, f'{tres}_{period_label}_residuals_{var_name}_only.png'))
    plt.close('all')


def plot_scatter(avar, period_label):
    """

    :param avar:
    :param period_label:
    :return:
    """
    print('SCATTERPLOTS')
    seas_name = seass[period_label]['name']
    fig, ax = plt.subplots(2, 2, figsize=(12, 12), dpi=300)
    axs = ax.ravel()

    for i, comp in enumerate(comps):
        x, y, vr_t_res = var_selection(avar, comp=comp)

        axs[i].set_ylabel(var_dict[comp]['label_uom'])
        try:
            print(f'plotting scatter VESPA-{var_dict[comp]['label']}')

            fig.suptitle(f'{var_name_u} {seas_name} {tres}', fontweight='bold')
            axs[i].set_title(var_dict[comp]['label'])

            time_list = pd.date_range(start=dt.datetime(years[0], 1, 1), end=dt.datetime(years[-1], 12, 31), freq=tres)

            x_all = x.reindex(time_list).fillna(np.nan)
            x_s = x_all.loc[(x_all.index.month.isin(seass[period_label]['months']))]
            y_all = y.reindex(time_list).fillna(np.nan)
            y_s = y_all.loc[(y_all.index.month.isin(seass[period_label]['months']))]
            idx = ~(np.isnan(x_s) | np.isnan(y_s))

            if seas_name != 'all':
                if var_dict[comp]['label'] == 'RS':
                    y_s = y.loc[(y.index.month.isin(seass[period_label]['months']))]
                    x_s = pd.Series(vr_t_res.reindex(y_s.index)[var_name])

                    idx = ~(np.isnan(x_s) | np.isnan(y_s))
                    axs[i].scatter(
                            x_s[idx].values, y_s[idx].values, s=50, facecolor='none', color=seass[period_label]['col'],
                            label=period_label)
                    axs[i].set_xlabel(var_dict[comp]['label_uom'])
                else:
                    axs[i].scatter(
                            x_s[idx], y_s[idx], s=5, color=seass[period_label]['col'], facecolor='none', alpha=0.5,
                            label=period_label)
                    axs[i].set_xlabel(var_dict['vr_t']['label_uom'])
            else:
                if var_dict[comp]['label'] == 'RS':
                    y_s = y.loc[(y.index.month.isin(seass[period_label]['months']))]
                    x_s = pd.Series(vr_t_res.reindex(y_s.index)[var_name])

                    idx = ~(np.isnan(x_s) | np.isnan(y_s))
                    axs[i].scatter(x_s[idx], y_s[idx], facecolor='none', s=50, color=seass[period_label]['col'])
                    axs[i].set_xlabel(var_dict[comp]['label_uom'])
                else:
                    bin_size = extr[var_name]['max'] / bin_nr
                    h = axs[i].hist2d(x_s[idx], y_s[idx], bins=bin_nr, cmap=plt.cm.jet, cmin=1, vmin=1)
                    axs[i].text(
                            0.10, 0.80, f'bin_size={bin_size} kg/m3',
                            transform=axs[i].transAxes)  # fig.colorbar(h[3], ax=axs[i], extend='both')
                    axs[i].set_xlabel(var_dict['vr_t']['label_uom'])

            if len(x_s[idx]) < 2 | len(y_s[idx]) < 2:
                print('ERROR, ERROR, NO DATA ENOUGH FOR PROPER FIT (i.e. only 1 point available)')
            else:
                b, a = np.polyfit(x_s[idx], y_s[idx], deg=1)
                xseq = np.linspace(extr[var_name]['min'], extr[var_name]['max'], num=1000)
                axs[i].plot(xseq, a + b * xseq, color='red', lw=2.5, ls='--', alpha=0.5)
                axs[i].plot(
                        [extr[var_name]['min'], extr[var_name]['max']], [extr[var_name]['min'], extr[var_name]['max']],
                        color='black', lw=1.5, ls='-')

                corcoef = ma.corrcoef(x_s[idx], y_s[idx])
                N = len(y_s[idx])
                rmse = np.sqrt(np.nanmean((y_s[idx] - x_s[idx]) ** 2))
                mbe = np.nanmean(y_s[idx] - x_s[idx])
                axs[i].text(
                        0.70, 0.80,
                        '$LWP_{WH}$: R=' + f'{corcoef[0, 1]:.2f}' + f' N={N:}' + f'\n y={b:+.2f}x{a:+.2f}' + f'\n MBE={mbe:.2f}' + f' RMSE={rmse:.2f}',
                        transform=ax.transAxes, fontsize=14, color='black', ha='left', va='center',
                        bbox=dict(facecolor='white', edgecolor='white'))

                axs[i].set_xlim(extr[var_name]['min'], extr[var_name]['max'])
                axs[i].set_ylim(extr[var_name]['min'], extr[var_name]['max'])
                axs[i].text(0.05, 0.95, letters[i] + ')', transform=axs[i].transAxes)

        except:
            print(f'error with {var_dict[comp]['label']}')

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
        seas_name = seass[period_label]['name']
        axs = ax.ravel()

        for i, comp in enumerate(comps):

            axs[i].set_ylabel(var_dict[comp]['label_uom'])
            x, y, vr_t_res = var_selection(avar, comp=comp)

            try:
                print(f'plotting scatter VESPA-{var_dict[comp]['label']}')

                fig.suptitle(f'{var_name_u} cumulative plot', fontweight='bold')
                axs[i].set_title(var_dict[comp]['label'])

                time_list = pd.date_range(
                        start=dt.datetime(years[0], 1, 1), end=dt.datetime(years[-1], 12, 31), freq=tres)

                x_all = x.reindex(time_list).fillna(np.nan)
                x_s = x_all.loc[(x_all.index.month.isin(seass[period_label]['months']))]
                y_all = y.reindex(time_list).fillna(np.nan)
                y_s = y_all.loc[(y_all.index.month.isin(seass[period_label]['months']))]
                idx = ~(np.isnan(x_s) | np.isnan(y_s))

                if seas_name != 'all':
                    if var_dict[comp]['label'] == 'RS':
                        y_s = y.loc[(y.index.month.isin(seass[period_label]['months']))]
                        x_s = pd.Series(vr_t_res.reindex(y_s.index)[var_name])

                        idx = ~(np.isnan(x_s) | np.isnan(y_s))
                        axs[i].scatter(
                                x_s[idx].values, y_s[idx].values, s=50, facecolor='none',
                                color=seass[period_label]['col'], label=period_label)
                        axs[i].set_xlabel(var_dict[comp]['label_uom'])

                    else:
                        axs[i].scatter(
                                x_s[idx], y_s[idx], s=5, color=seass[period_label]['col'], edgecolors='none', alpha=0.5,
                                label=period_label)
                        axs[i].set_xlabel(var_dict['vr_t']['label_uom'])

                if len(x_s[idx]) < 2 | len(y_s[idx]) < 2:
                    print('ERROR, ERROR, NO DATA ENOUGH FOR PROPER FIT (i.e. only 1 point available)')
                else:
                    b, a = np.polyfit(x_s[idx], y_s[idx], deg=1)
                    xseq = np.linspace(extr[var_name]['min'], extr[var_name]['max'], num=1000)
                    axs[i].plot(xseq, a + b * xseq, color=seass[period_label]['col'], lw=2.5, ls='--', alpha=0.5)
                    axs[i].plot(
                            [extr[var_name]['min'], extr[var_name]['max']],
                            [extr[var_name]['min'], extr[var_name]['max']], color='black', lw=1.5, ls='-')

                    axs[i].set_xlim(extr[var_name]['min'], extr[var_name]['max'])
                    axs[i].set_ylim(extr[var_name]['min'], extr[var_name]['max'])
                    axs[i].text(0.05, 0.95, letters[i] + ')', transform=axs[i].transAxes)
                    axs[i].legend()
            except:
                print(f'error with {var_dict[comp]['label']}')

    plt.savefig(os.path.join(basefol_out, tres, f'{tres}_scatter_cum_{var_name}_only.png'))
    plt.close('all')


def var_selection(avar, comp):
    [vr_c, vr_e, vr_l, vr_t, vr_t1, vr_t2, vr_c_res, vr_e_res, vr_l_res, vr_t_res, vr_t1_res, vr_t2_res] = avar

    x = vr_t_res[var_name].resample(tres).mean()

    if comp == 'vr_c':
        try:
            y = vr_c_res[var_name]
        except KeyError:
            print(f'error with {var_dict[comp]['label']}')
    if comp == 'vr_e':
        try:
            y = vr_e_res[var_name]
        except KeyError:
            print(f'error with {var_dict[comp]['label']}')
    if comp == 'vr_l':
        try:
            y = vr_l_res[var_name]
        except KeyError:
            print(f'error with {var_dict[comp]['label']}')
    if comp == 'vr_t':
        try:
            y = vr_t_res[var_name]
        except KeyError:
            print(f'error with {var_dict[comp]['label']}')
    if comp == 'vr_t1':
        try:
            y = vr_t1_res[var_name]
        except KeyError:
            print(f'error with {var_dict[comp]['label']}')
    if comp == 'vr_t2':
        try:
            y = vr_t2_res[var_name].dropna()
        except KeyError:
            print(f'error with {var_dict[comp]['label']}')

    return x, y, vr_t_res
