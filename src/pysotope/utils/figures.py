import matplotlib.pyplot as plt
import numpy as np
import os
from .regression import *
from .curve_fitting import *


# def std_plot(lin, drift, folder_path, fig_path, isotope, cutoff_line=None, regress=False, dD = "dD"):
#     """
#     Function to plot linearity and drift standards.
#     ~GAO~12/1/2023
#     """
#     fig, ax = plt.subplots(2, 2, figsize=[6, 4], sharex=False)
#     for i in [0, 1]:
#         drift.chain.unique()
#         temp = drift[drift.chain == drift.chain.unique()[i]]
#         ax[i,0].scatter(temp["date-time_true"], temp[dD], alpha=0.4, ec='k', s=80, c='blue')
#         ax[i,0].text(0.9, 0.9, drift.chain.unique()[i],
#                    horizontalalignment='center', verticalalignment='center',
#                    transform=ax[i,0].transAxes)

#         temp = lin[lin.chain == lin.chain.unique()[i]]
#         ax[i,1].scatter(temp["area"], temp[dD], alpha=0.2, ec='k', s=80, c='orange')
#         ax[i,1].text(0.9, 0.9, lin.chain.unique()[i],
#                       horizontalalignment='center', verticalalignment='center',
#                       transform=ax[i,1].transAxes)
#     # Set x-axis labels for the third subplot (ax2)
#     ax[1,0].set_xlabel('Date (mm-dd-yyyy)')
#     ax[1,0].set_xticks(ax[1,0].get_xticks())
#     ax[1,0].set_xticklabels(labels=ax[1,0].get_xticklabels(), rotation=45)
#     ax[0,0].set_xticks([]);#ax[1].set_xticks([])

#     ax[1,1].set_xlabel('Peak Area (mVs)')
#     ax[0,0].set_title("Drift Standards")
#     ax[0,1].set_title("Linearity Standards")
#     if isotope == "dD": label = "dD"
#     else: label = "dC"
#     fig.supylabel('Normalized '+str(label)+' (‰)')

#     # Plot user-defined cutoff line
#     if cutoff_line is not None:
#         ax[0,1].axvline(cutoff_line[0], c='red', linestyle='--')
#         ax[1,1].axvline(cutoff_line[1], c='red', linestyle='--')
#     # plt.show(block=True)
#     plt.tight_layout()
#     plt.savefig(os.path.join(fig_path, 'Standards Raw.png'), dpi=300, bbox_inches='tight')
#     plt.show()
def std_plot(lin, drift, folder_path, fig_path, isotope, cutoff_line=None, regress=False, dD="dD"):
    """
    Plot drift (left) and linearity (right) standards.
    Robust to 1 or 2+ chains; shows up to 2 per type in a 2x2 layout.
    """

    fig, ax = plt.subplots(2, 2, figsize=[6, 4], sharex=False)

    drift_chains = list(drift["chain"].dropna().unique())
    lin_chains   = list(lin["chain"].dropna().unique())

    # always attempt up to 2 rows
    for i in range(2):

        # --- Drift column (0) ---
        if i < len(drift_chains):
            ch = drift_chains[i]
            temp = drift[drift["chain"] == ch]
            ax[i, 0].scatter(temp["date-time_true"], temp[dD], alpha=0.4, ec="k", s=80, c='blue')
            ax[i, 0].text(
                0.9, 0.9, str(ch),
                horizontalalignment="center", verticalalignment="center",
                transform=ax[i, 0].transAxes
            )
        else:
            ax[i, 0].axis("off")

        # --- Linearity column (1) ---
        if i < len(lin_chains):
            ch = lin_chains[i]
            temp = lin[lin["chain"] == ch]
            ax[i, 1].scatter(temp["area"], temp[dD], alpha=0.2, ec="k", s=80, c='orange')
            ax[i, 1].text(
                0.9, 0.9, str(ch),
                horizontalalignment="center", verticalalignment="center",
                transform=ax[i, 1].transAxes
            )
        else:
            ax[i, 1].axis("off")

    # Labels/titles (only for axes that exist)
    if len(drift_chains) > 0:
        ax[min(1, len(drift_chains)-1), 0].set_xlabel("Date (mm-dd-yyyy)")
        ax[min(1, len(drift_chains)-1), 0].set_xticks(ax[min(1, len(drift_chains)-1), 0].get_xticks())
        ax[min(1, len(drift_chains)-1), 0].set_xticklabels(
            labels=ax[min(1, len(drift_chains)-1), 0].get_xticklabels(), rotation=45
        )

    if len(lin_chains) > 0:
        ax[min(1, len(lin_chains)-1), 1].set_xlabel("Peak Area (mVs)")

    if ax[0, 0].has_data():
        ax[0, 0].set_title("Drift Standards")
    if ax[0, 1].has_data():
        ax[0, 1].set_title("Linearity Standards")
    ax[0,0].set_xticks([]);#ax[1].set_xticks([])
    label = "dD" if isotope == "dD" else "dC"
    fig.supylabel("Normalized " + str(label) + " (‰)")

    # Cutoff lines: draw on whichever linearity rows exist
    if cutoff_line is not None and len(lin_chains) > 0:
        if len(cutoff_line) >= 1 and ax[0, 1].has_data():
            ax[0, 1].axvline(cutoff_line[0], c="red", linestyle="--")
        if len(cutoff_line) >= 2 and ax[1, 1].has_data():
            ax[1, 1].axvline(cutoff_line[1], c="red", linestyle="--")

    plt.tight_layout()
    plt.savefig(os.path.join(fig_path, "Standards Raw.png"), dpi=300, bbox_inches="tight")
    plt.show()

def verify_lin_plot(lin, samples, fig_path, dD_id, log_file_path,cutoff_line, isotope, regress=False, include_parabolic=False):
    """
    Function to plot linearity standards with color differentiation based on cutoff,
    but split into one subplot per unique 'chain' in the `samples` DataFrame.
    ~GAO~12/4/2023 (modified)
    """
    # Ensure cutoff_line is a float
    cutoff_line = float(cutoff_line)

    # Find all unique chains in samples
    chains = samples["chain"].unique()
    n_chains = len(chains)
    if n_chains == 0:
        raise ValueError("`samples.chain.unique()` returned no chains. Nothing to plot.")

    # Create a figure with n_chains subplots, in a single column
    fig, axes = plt.subplots(n_chains, 1, figsize=(5, 3 * n_chains), sharex=False)
    # If there's only one chain, axes won't be an array; force it into an array
    if n_chains == 1:
        axes = np.array([axes])

    # Precompute all points above/below cutoff (these are the same for every chain)
    above_cutoff = lin[lin["area"] >= cutoff_line]
    below_cutoff = lin[lin["area"] <  cutoff_line]

    # Fit model once (on the "above_cutoff" points); we'll re‐plot the same best‐fit curve in each subplot
    xdata = above_cutoff["area"].values
    ydata = above_cutoff[dD_id].values
    best_model, popt, sse, pcov = fit_and_select_best(xdata, ydata, include_parabolic)
    x_fit = np.linspace(xdata.min(), xdata.max(), 200)

    if best_model == "linear":
        y_fit = linear_func(x_fit, *popt)
        model_label = "Linear Fit"
        parameter_text = f"y = {popt[0]:.3g} x + {popt[1]:.3g}"
    elif best_model == "decay":
        y_fit = exp_decay(x_fit, *popt)
        model_label = "Exponential Decay"
        parameter_text = f"y = {popt[0]:.3g} e^(−{popt[1]:.3g} x + {popt[2]:.3g})"
    elif best_model == "growth":
        y_fit = exp_growth(x_fit, *popt)
        model_label = "Exponential Growth"
        parameter_text = f"y = {popt[0]:.3g} (1 − e^(−{popt[1]:.3g} x) + {popt[2]:.3g})"
    elif best_model == "parabolic":
        y_fit = parabolic_func(x_fit, *popt)
        model_label = "Parabolic"
        parameter_text = f"y = {popt[0]:.6g} · x^2 + {popt[1]:.6g} · x + {popt[2]:.6g}"
    else:
        raise RuntimeError("No model converged in fit_and_select_best().")

    # Compute R²
    tss = np.sum((ydata - ydata.mean()) ** 2)
    r_squared = 1.0 if tss == 0 else 1.0 - (sse / tss)

    # Now loop over each chain, make a subplot
    for idx, chain in enumerate(chains):
        ax = axes[idx]

        # 1) Scatter points: above cutoff (orange) and below cutoff (grey)
        ax.scatter(above_cutoff["area"], above_cutoff[dD_id],  alpha=0.4,  ec="k", s=80, c="orange", label="Above peak area threshold")
        ax.scatter(below_cutoff["area"], below_cutoff[dD_id], alpha=0.4, ec="k", s=80, c="grey",label="Below peak area threshold")
        # Vertical cutoff line
        ax.axvline(cutoff_line, color="red", linestyle="--")

        # 2) Vertical lines for this specific chain
        subset = samples[samples["chain"] == chain]
        for i, row in subset.iterrows():
            ax.axvline(row["area"], color="k", alpha=0.3, zorder=0,label="Samples" if i == subset.index[0] else None)

        # 3) Plot the best-fit curve (same for all subplots)
        ax.plot(x_fit, y_fit, "red", label=model_label)
        ax.text(0.01,1.1, str(chain), transform=ax.transAxes, ha='left', va='top', fontsize=10)

        # 5) Axis labels (only bottom‐most subplot gets the X label; all get a Y label)
        if isotope == "dD":
            y_label = "Normalized dD (‰)"
        else:
            y_label = "Normalized dC (‰)"
        ax.set_ylabel(y_label)

        if idx == n_chains - 1:
            ax.set_xlabel("Peak Area (mVs)")

        # 6) Minor formatting
        ax.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.2),
            frameon=False,
            fancybox=False,
            shadow=False,
            ncol=2
        )
        ax.grid(alpha=0)

    # Adjust layout and save
    plt.tight_layout()
    out_fname = os.path.join(fig_path, "Linearity_by_chain.png")
    plt.savefig(out_fname, bbox_inches="tight", dpi=300)
    plt.show()

    # Finally, print fit diagnostics once
    print(f"Chosen Model: {model_label}")
    print(f"Parameters: {parameter_text}")
    print(f"R² = {r_squared:.3f} | SSE = {sse:.3f}")

def total_dD_correction_plot(uncorrected_unknown, unknown , folder_path, fig_path, isotope):
    unique_chains = unknown['Chain Length'].unique()
    num_chains    = len(unique_chains)
    if isotope == "dD": label = "dD"
    else: label = "dC"
    if num_chains>1:
        fig, axes     = plt.subplots(num_chains, 1, figsize=(5,3 * num_chains))
        for i, chain in enumerate(unique_chains):
            chain_unknown             = unknown[unknown['Chain Length'] == chain]
            uncorrected_chain_unknown = uncorrected_unknown[uncorrected_unknown.Component == chain]
            if "Raw "+str(isotope) in uncorrected_chain_unknown:
                axes[i].scatter(uncorrected_chain_unknown['Peak area'], uncorrected_chain_unknown[f'Raw {isotope}'], label=f'Original {isotope}', marker = 'x', alpha=0.6, s=60, c='k')

            if 'Final - Methanol Corrected '+str(isotope) in chain_unknown:
                axes[i].errorbar(chain_unknown["Mean Area"], chain_unknown['Final - Methanol Corrected '+str(isotope)],
                                 yerr=chain_unknown['Total Uncertainty'],
                                 linestyle="", fmt='', ecolor='red', alpha=0.5)
                axes[i].scatter(chain_unknown["Mean Area"], chain_unknown['Final - Methanol Corrected '+str(isotope)], label='Corrected '+str(label),alpha=0.6, edgecolor='k', s=60, color='red')
            axes[i].set_title(f'Chain: {chain}')
            axes[i].set_xlabel('Peak Area (mVs)')
            axes[i].set_ylabel('Normalized '+str(label)+' (‰)')
            if i == len(unique_chains):
                axes[i].legend(loc='upper left', bbox_to_anchor=(1, 1))
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', ncol=2, bbox_to_anchor=(0.5, 0))
    else:
        plt.scatter(uncorrected_unknown['Peak area'], uncorrected_unknown[f'Raw {isotope}'], label='Original '+str(label), marker = 'x', alpha=0.6, s=60, c='k')
        plt.errorbar(unknown["Mean Area"], unknown['Final - Methanol Corrected '+str(isotope)],
                                 yerr=unknown['Total Uncertainty'],
                                 linestyle="", fmt='', ecolor='red', alpha=0.5)
        plt.scatter(unknown["Mean Area"], unknown['Final - Methanol Corrected '+str(isotope)], label='Corrected '+str(label),alpha=0.6, edgecolor='k', s=60, color='red')
        plt.xlabel('Peak Area (mVs)')
        plt.ylabel(str(label)+' (‰)')
        plt.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.5, 0))
    plt.subplots_adjust(hspace=1)
    plt.savefig(os.path.join(fig_path, 'isotope_corrections.png'), bbox_inches='tight')
    plt.close()

def drift_std_corr(norm, isotope, drift_std, t_mean, intercept, slope, fig_path):
        fig, ax = plt.subplots(figsize=(6, 4))
        for i, magic_markers in zip(norm.chain.unique(),['o','s']):
            temp = norm[norm.chain==i]
            ax.scatter(temp["time_rel"], temp[isotope], marker=magic_markers,
                         c="k", ec='k', alpha=0.5, label=f"{i} Drift Raw")
            ax.scatter(temp["time_rel"], temp['corrected_norm'], marker=magic_markers,
                       c="red", ec='k', label= f"{i} Drift Corrected")
        t_line = np.linspace(drift_std["time_rel"].min(), drift_std["time_rel"].max(), 50)
        ax.plot(t_line, (slope * (t_line - t_mean) + intercept), "--k", lw=1, label="WLS fit")
        ax.set_xlabel("Centered Time")
        ax.set_ylabel(f"Centered {isotope} (‰)")
        ax.legend()
        fig.tight_layout()
        plt.savefig(os.path.join(fig_path, 'Drift.png'), bbox_inches='tight')
        plt.show()


def standard_check_figures(cfg, stds, fig_path, label, vsmow, isotope):
    cl = vsmow[vsmow['VSMOW accuracy check']==True]['chain length'].values[0]
    cl_val = vsmow[vsmow['VSMOW accuracy check']==True]['isotope value'].values[0]
    fig = plt.figure()
    plt.title(f"Raw {label}")
    for j in vsmow['chain length'].unique():
        temp = stds[stds.chain==j]
        plt.scatter(temp['VSMOW_dD_actual'], temp[label], label = j)
    plt.ylabel(f"Measured {label}")
    plt.xlabel(f"VSMOW {label}")
    plt.scatter([cl_val]*len(stds[stds.chain==f"{cl}"]), stds[stds.chain==f"{cl}"]['label'], label = cl)
    plt.legend()
    plt.savefig(f"{fig_path}/Standards_RawVsVSMOW.png", dpi=300)
    plt.close()

    if cfg.drift_applied:
        fig = plt.figure()
        plt.title("Drift Corrected")
        for j in vsmow['chain length'].unique():
            temp = stds[stds.chain==j]
            plt.scatter(temp[f'VSMOW_{label}_actual'], temp[f'drift_corrected_{label}'], label =j)
        plt.scatter([cl_val]*len(stds[stds.chain==f"{cl}"]), stds[stds.chain==f"{cl}"][f'drift_corrected_{label}'], label= cl)
        plt.legend()
        plt.ylabel(f"Measured {label}")
        plt.xlabel(f"VSMOW {label}")
        plt.savefig(f"{fig_path}/Standards_DriftCorrVsVSMOW.png", dpi=300)
        plt.close()
    if cfg.linearity_applied:
        fig = plt.figure()
        plt.title("Linearity Corrected")
        for j in vsmow['chain length'].unique():
            temp = stds[stds.chain==j]
            plt.scatter(temp[f'VSMOW_{label}_actual'], temp[f'linearity_corrected_{label}'], label =j)
        plt.scatter([cl_val]*len(stds[stds.chain==f"{cl}"]), stds[stds.chain==f"{cl}"][f'linearity_corrected_{label}'], label = cl)
        plt.legend()
        plt.ylabel(f"Measured {label}")
        plt.xlabel(f"VSMOW {label}")
        plt.savefig(f"{fig_path}/Standards_LinearityCorrVsVSMOW.png", dpi=300)
        plt.close()




