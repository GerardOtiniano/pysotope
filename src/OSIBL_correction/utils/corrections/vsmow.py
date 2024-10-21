import pandas as pdimport numpy as npimport osimport matplotlib.pyplot as pltfrom . .regression import wls_regressionimport statsmodels.api as smfrom statsmodels.sandbox.regression.predstd import wls_prediction_std# VSMOWdef vsmow_correction(unknown, lin, drift, correction_log, folder_path,fig_path, log_file_path, isotope, fig=True):    print("Applying VSMOW correction")    if isotope == "dD":         label="δD"        vsmow_data = {        "id"      : ["C18", "C20Z2", "C20", "C28", "C24"],        "dD" : [-206.2, -4.9, -166.7, -89.28190994, -179.3],        "std"   : [1.7, 1, 0.3, 1.062744893, 1.7],        "n"       : [5, 6, 3, 924, 5]        }        standards_iso = ["C18", "C20", "C28"]    else:         label = "δC"        vsmow_data = {        "id"      : ["C18", "C20", "C24"],        "dD" : [-23.24, -30.68, -26.57],        "std"   : [0.01, 0.02, 0.02]        }        standards_iso = ["C18", "C20", "C24"]    vsmow = pd.DataFrame(vsmow_data)    stds  = pd.concat([lin, drift])        for i in standards_iso:        mask = stds.chain==i        stds.loc[mask, 'VSMOW_dD_actual']    = vsmow.loc[vsmow.id == i, 'dD'].iloc[0]        stds.loc[mask, 'vsmow_error_actual'] = vsmow.loc[vsmow.id == i, 'std'].iloc[0]    # Need to revise to pass the proper δD correction from previous runs    # Determine which column to use based on whether linearity correction was applied    dD_id = "linearity_corrected_dD" if correction_log.loc["Linearity", 'samples'] == 1 else "drift_corrected_dD"    stds=stds[~stds[dD_id].isna()]    mask = stds.chain.isin(standards_iso)    m_d, b_d, r_squared, p_value, std_error, model = wls_regression(stds.loc[mask, dD_id],stds.loc[mask,'VSMOW_dD_actual'], log_file_path)    # Predict VSMOW values - standards    """     Predicted error is: standard error of the prediction, from the scipy package wls_prediction_std.     This is the correct prediction error assocaited with the prediction of weighted least squares (wls)    regression.    """    new_x = sm.add_constant(np.array(stds[dD_id]))    stds['VSMOW_dD'] = model.predict(new_x)    prstd, std_pred_ci_lower, std_pred_ci_upper = wls_prediction_std(model, exog = new_x)#, alpha=0.05)    stds['VSMOW_error'] = prstd    # Predict VSMOW values - samples    new_x = sm.add_constant(np.array(unknown[dD_id]))    unknown['VSMOW_dD'] = model.predict(new_x)    prstd, unknown_pred_ci_lower, unknown_pred_ci_upper = wls_prediction_std(model, exog = new_x, alpha=0.05)    unknown['VSMOW_error'] = prstd    relative_standard = -179.3    if fig:        mask = stds.chain.isin(["C18","C20","C28"])        plt.scatter(stds.loc[mask,dD_id], stds.loc[mask,'VSMOW_dD_actual'], c='red', ec='k', label='C18, 20, 28 standards')        plt.plot(stds.loc[mask, dD_id], stds.loc[mask, dD_id]*m_d + b_d, c='k', linestyle='-', alpha=0.6)        if isotope=='dD':            mask = stds.chain=="C24"            plt.scatter(stds.loc[mask,dD_id], stds.loc[mask,'VSMOW_dD'], c='green', ec='k', label='Predicted C24')            plt.axhline(relative_standard,c='green',linestyle='--',zorder=0,label="Actual C24")        # plt.plot(new_x, unknown_pred_ci_lower, c='red')        # plt.plot(new_x, unknown_pred_ci_upper, c='red')        plt.scatter(unknown[isotope], unknown['VSMOW_dD'], c='blue', marker='x', alpha=0.6, label = 'corrected samples')        plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), frameon=False, ncol=2)                plt.xlabel("Measured "+str(label)+" (‰)")        plt.ylabel("Predicted VSMOW "+str(label)+" (‰)")        plt.savefig(os.path.join(fig_path, 'VSMOW correction.png'), dpi=300, bbox_inches='tight')        plt.show()            print("\nRegression statistics VSMOW standards:")    print("Linear equation: "+str(label)+f" = ({m_d:.2f})(time) + {b_d:.2f}")    print(f"Adjusted R²: {r_squared:.2f}")    print(f"P-value: {p_value:.2f}")    print(f"Standard Error: {std_error:.2f}")    print("\nAccuracy of C24 prediction")    mask = stds.chain=="C24"    test = stds.loc[mask,'VSMOW_dD']    mae = (test - relative_standard).abs().mean()    rmse = np.sqrt(((test - (relative_standard)) ** 2).mean())    print(f"MAE: {mae:.2f}")    print(f"RMSE: {rmse:.2f}")    # Uncertainty for samples    #stds['vsmow_error'] = std_error    if 'drift_error' not in stds.columns:        stds['drift_error'] = np.nan    if 'linearity_error' not in stds.columns:        stds['linearity_error'] = np.nan    return unknown, stds