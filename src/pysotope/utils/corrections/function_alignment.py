#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 14:42:43 2025

@author: gerard
"""

import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess as _lowess
from . .queries import *

def _choose_x_star_overlap(df, area_col="area", chain_col="chain"):
    mins = df.groupby(chain_col)[area_col].min()
    maxs = df.groupby(chain_col)[area_col].max()
    xmin_star = mins.max()
    xmax_star = maxs.min()
    if not np.isfinite(xmin_star) or not np.isfinite(xmax_star) or xmin_star >= xmax_star:
        return None
    return 0.5 * (xmin_star + xmax_star)  # midpoint of overlap

def _level_at_xstar(x, y, x_star):
    """Estimate y(x_star) using LOWESS if available; otherwise monotone interp."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3:
        return None
    order = np.argsort(x)
    x, y = x[order], y[order]
    sm = _lowess(y, x, frac=0.6, it=0, return_sorted=True)
    xs, ys = sm[:, 0], sm[:, 1]
    return float(np.interp(x_star, xs, ys))

def build_norm_with_fallback(lin_std, y_col, area_col="area", chain_col="chain", name_col ="Identifier 1", log_file_path=None):
    """
    Try 1) Anchor alignment at a common area x*: for each chain g, estimate y_g(x*)
           and shift by a constant so all chains coincide at x*.
       If that fails, 2) Moment/centering fallback: per-chain mean-centering.

    Returns
    -------
    norm : DataFrame with [area_col, chain_col, y_col] where y_col is shifted
    meta : dict with method info and diagnostics
    """
    df = lin_std[[name_col, area_col, chain_col, y_col]].copy()
    x_star = _choose_x_star_overlap(df, area_col=area_col, chain_col=chain_col)
    if x_star is not None:
        levels = {}
        ok_all = True
        for g, sub in df.groupby(chain_col):
            lvl = _level_at_xstar(sub[area_col], sub[y_col], x_star)
            if lvl is None or not np.isfinite(lvl):
                ok_all = False
                break
            levels[g] = float(lvl)
        if ok_all and len(levels) == df[chain_col].nunique():
            R0 = float(np.median(list(levels.values())))
            norm = df.copy()
            norm[y_col] = norm.apply(lambda r: r[y_col] - levels[r[chain_col]] + R0, axis=1)
            append_to_log(log_file_path, f"- Linearity normalization: anchor alignment at x*={x_star:.4g} (R0=median levels)")
            return norm, {"method": "anchor", "x_star": x_star, "R0": R0, "levels_at_xstar": levels}

    chain_means = df.groupby(chain_col)[y_col].mean().to_dict()
    global_center = float(np.median(list(chain_means.values())))
    norm = df.copy()
    norm[y_col] = norm.apply(lambda r: r[y_col] - chain_means[r[chain_col]] + global_center, axis=1)
    append_to_log(log_file_path, "- Linearity normalization: fallback to per-chain mean centering (moment matching)")
    print(norm)
    norm[y_col] = norm[y_col]+norm[y_col].min()+1
    return norm, {"method": "moment_center", "global_center": global_center, "chain_means": chain_means}