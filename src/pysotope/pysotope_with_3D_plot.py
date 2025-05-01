# src/OSIBL_correction
import os 

from .utils.corrections.drift import *
from .utils.corrections.linearity import *
from .utils.corrections.methanol import *
from .utils.corrections.vsmow import *
from .utils.outliers.outliers import *
from .utils.queries import *
from .utils.regression import *
from .utils.uncertainty_and_output import *
from .utils.figures import *
from .utils.base_functions import *
from .utils.config import CorrectionConfig




def iso_process(pame=False, user_linearity_conditions = False):
    cfg = CorrectionConfig()
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    import os
    from IPython.display import clear_output
    from scipy.stats import linregress
    from matplotlib.dates import date2num
    from datetime import datetime, timedelta
    import time
    from scipy.stats import zscore
    from sklearn.linear_model import HuberRegressor
    import scipy.stats as stats

    # Query isotope system
    isotope = isotope_type()
    
    # Setup output folder
    folder_path, fig_path, results_path, loc, log_file_path = create_folder(isotope)
    
    # Set standards
    standards_df = load_standards(isotope)#query_stds(alt_stds, isotope)
    append_to_log(log_file_path, standards_df)
    
    # Import data
    lin_std, drift_std, samples, correction_log, pame = import_data(loc, folder_path, log_file_path, isotope, standards_df)
    uncorrected_samples = samples.copy()
    
    from sklearn.linear_model import LinearRegression
    from mpl_toolkits.mplot3d import Axes3D
    import plotly.express as px
    import plotly.graph_objects as go
    
    lin_std["dD_norm"] = (lin_std["dD"]/lin_std.groupby(["Identifier 1", "chain"])["dD"].transform("mean"))

    drift_std["dD_norm"] = (drift_std["dD"]/drift_std.groupby(["Identifier 1", "chain"])["dD"].transform("mean"))
    
    y  = np.concatenate([lin_std['dD_norm'].values,     drift_std['dD_norm'].values])
    x1 = np.concatenate([lin_std['time_rel'].values,
                         drift_std['time_rel'].values])
    x2 = np.concatenate([lin_std['area'].values,
                         drift_std['area'].values])
    X = np.column_stack((x1, x2))          # shape (n_samples, 2)
    model = LinearRegression().fit(X, y)
    print("β0 (intercept):", model.intercept_)
    print("β1, β2 (slopes):", model.coef_)
    
    # ── 3.  BUILD A GRID OVER (x1, x2) FOR THE REGRESSION PLANE  ─────────────
    x1_lin = np.linspace(x1.min(), x1.max(), 30)
    x2_lin = np.linspace(x2.min(), x2.max(), 30)
    x1_grid, x2_grid = np.meshgrid(x1_lin, x2_lin)           # 30×30 grid
    X_grid = np.column_stack((x1_grid.ravel(), x2_grid.ravel()))
    y_pred_grid = model.predict(X_grid).reshape(x1_grid.shape)
    
    df = pd.DataFrame({'x1': x1, 'x2': x2, 'y': y})
    fig = px.scatter_3d(df, x='x1', y='x2', z='y',
                        opacity=0.7, title='Multiple-Linear Regression fit')
    
    # 2) regression plane
    fig.add_trace(
        go.Surface(
            x=x1_grid,          # 2-D arrays
            y=x2_grid,
            z=y_pred_grid,
            showscale=False,
            opacity=0.4,
            name='Fitted plane'
        )
    )
    
    # 3) axis labels, view angle (optional)
    fig.update_layout(
        scene=dict(
            xaxis_title='x1',
            yaxis_title='x2',
            zaxis_title='y'
        ),
        legend=dict(yanchor='top', y=0.99, xanchor='left', x=0.01)
    )
    
    fig.show()
    # print(list(lin_std))
    # lin_std["dD_norm"] = (lin_std["dD"]/lin_std.groupby(["Identifier 1", "chain"])["dD"].transform("mean"))

    # drift_std["dD_norm"] = (drift_std["dD"]/drift_std.groupby(["Identifier 1", "chain"])["dD"].transform("mean"))
    
    # fig = go.Figure()                       # one figure for everything
    # colors = ['royalblue', 'crimson']       # optional colour list
    
    # for i, df_x in enumerate([lin_std, drift_std]):
    #     # -- 1. grab data -----------------------------------------------------
    #     y  = df_x['dD_norm'].values
    #     x1 = df_x['time_rel'].values
    #     x2 = df_x['area'].values
    
    #     # -- 2. fit regression ------------------------------------------------
    #     X = np.column_stack((x1, x2))
    #     model = LinearRegression().fit(X, y)
    
    #     # -- 3. make prediction surface --------------------------------------
    #     x1_lin = np.linspace(x1.min(), x1.max(), 30)
    #     x2_lin = np.linspace(x2.min(), x2.max(), 30)
    #     x1_grid, x2_grid = np.meshgrid(x1_lin, x2_lin)
    #     y_pred_grid = model.predict(
    #         np.column_stack((x1_grid.ravel(), x2_grid.ravel()))
    #     ).reshape(x1_grid.shape)
    
    #     # -- 4. add scatter trace --------------------------------------------
    #     fig.add_trace(
    #         go.Scatter3d(
    #             x=x1, y=x2, z=y,
    #             mode='markers',
    #             marker=dict(size=4, color=colors[i]),
    #             name=f'Data set {i+1}'
    #         )
    #     )
    
    #     # -- 5. add surface trace --------------------------------------------
    #     fig.add_trace(
    #         go.Surface(
    #             x=x1_grid, y=x2_grid, z=y_pred_grid,
    #             showscale=False,
    #             opacity=0.4,
    #             colorscale=[[0, colors[i]], [1, colors[i]]],
    #             name=f'Plane {i+1}'
    #         )
    #     )
    
    # # -- 6. tidy up axes / title ---------------------------------------------
    # fig.update_layout(
    #     title='Multiple-linear regressions for two data sets',
    #     scene=dict(
    #         xaxis_title='time_rel',
    #         yaxis_title='area',
    #         zaxis_title='dD_norm'
    #     )
    # )
    
    # fig.show()
    
    
    
    # # Run standard plots for area
    # std_plot(lin_std, drift_std, folder_path=folder_path, fig_path=fig_path,isotope=isotope, dD=isotope)
    
    # # Drift Correction
    # samples, lin_std, drift_std, dD_temp, correction_log = process_drift_correction(cfg, samples, lin_std, drift_std, correction_log, log_file_path=log_file_path, fig_path=fig_path,isotope=isotope)
    
    # # # Show plots again
    # # std_plot(lin_std, drift_std, folder_path=folder_path, fig_path=fig_path, dD=dD_temp,isotope=isotope)

    # # Linearity (area) correction
    # drift_std, correction_log, lin_std, samples = process_linearity_correction(cfg, samples, drift_std, lin_std, dD_temp, correction_log, folder_path, fig_path, isotope, user_linearity_conditions, log_file_path=log_file_path)
 
    # # VSMOW correction
    # samples, standards = vsmow_correction(cfg, samples, lin_std, drift_std, correction_log, folder_path, fig_path, log_file_path, isotope, standards_df)

    # # Methylation Correction
    # if isotope =="dD":
    #     samples, standards = q_methylation(samples, standards, log_file_path);

    # # PAME
    # if pame:
    #     samples, pame_unknown = calculate_methanol_dD(samples, isotope, log_file_path)
        
    # # Remove outliers
    # samples, excluded_samples = outlier_removal(samples, fig_path, log_file_path)
    # raw_samples = samples

    # # Calculate mean values of replicate analyses
    # samples = mean_values_with_uncertainty(samples, cfg, sample_name_header="Identifier 1", chain_header="chain", iso=isotope)
    # if pame:
    #     pame_unknown = mean_values_with_uncertainty(pame_unknown, sample_name_header="Identifier 1", chain_header="chain", iso=isotope)
    # else:
    #     pame_unknown = None
    
    # # Final Data Correction and Plot
    # output_results(raw_samples, samples, standards, pame_unknown, folder_path, fig_path, results_path, isotope, pame)
    
    
    
    
    
    