import matplotlib.pyplot as plt
import numpy as np
import os
from .regression import *
from .curve_fitting import *


def std_plot(lin, drift, folder_path, fig_path, isotope, cutoff_line=None, regress=False, dD = "dD"):
    """
    Function to plot linearity and drift standards.
    ~GAO~12/1/2023
    """
    fig, axes = plt.subplots(2, 2, figsize=[6, 4], sharex=False)
    chain_length = ["C18", "C20", "C24", "C28", "C30"]
    ax = axes.flatten()

    for i in [0, 2]:
        temp = drift[drift.chain == chain_length[i]]
        ax[i].scatter(temp["date-time_true"], temp[dD], alpha=0.4, ec='k', s=80, c='blue')
        ax[i].text(0.9, 0.9, chain_length[i],
                   horizontalalignment='center', verticalalignment='center',
                   transform=ax[i].transAxes)

        temp = lin[lin.chain == chain_length[i + 1]]
        ax[i + 1].scatter(temp["area"], temp[dD], alpha=0.2, ec='k', s=80, c='orange')
        ax[i + 1].text(0.9, 0.9, chain_length[i + 1],
                      horizontalalignment='center', verticalalignment='center',
                      transform=ax[i + 1].transAxes)
    # Set x-axis labels for the third subplot (ax2)
    ax[2].set_xlabel('Date (mm-dd-yyyy)')
    ax[2].set_xticks(ax[2].get_xticks())
    ax[2].set_xticklabels(labels=ax[2].get_xticklabels(), rotation=45)
    ax[0].set_xticks([]);#ax[1].set_xticks([])

    ax[3].set_xlabel('Peak Area (mVs)')
    ax[0].set_title("Drift Standards")
    ax[1].set_title("Linearity Standards")
    if isotope == "dD": label = "dD"
    else: label = "dC"
    fig.supylabel('Normalized '+str(label)+' (‰)')

    # Plot user-defined cutoff line
    if cutoff_line is not None:
        ax[1].axvline(cutoff_line[0], c='red', linestyle='--')
        ax[3].axvline(cutoff_line[1], c='red', linestyle='--')
    # plt.show(block=True)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_path, 'Standards Raw.png'), dpi=300, bbox_inches='tight')
    plt.show()
    
def verify_lin_plot(lin, fig_path, dD_id, log_file_path,cutoff_line, isotope, regress=False):
    """
    Function to plot linearity and drift standards with color differentiation based on cutoff.
    ~GAO~12/4/2023
    """
    cutoff_line=float(cutoff_line)
    fig = plt.figure(figsize=[5, 3])
    above_cutoff = lin[lin["area"] >= cutoff_line]
    plt.scatter(above_cutoff["area"], above_cutoff[dD_id], alpha=0.4, ec='k', s=80, c='orange', label = "Above peak area threshold.")
    below_cutoff = lin[lin["area"] < cutoff_line]
    plt.axvline(cutoff_line,color='red',linestyle="--")
    plt.scatter(below_cutoff["area"], below_cutoff[dD_id], alpha=0.4, ec='k', s=80, c='grey', label = "Below peak area threshold.")
    if isotope == "dD": label = "dD"
    else: label = "dC"
    plt.ylabel("Normalized "+str(label)+" (‰)")
    plt.xlabel('Peak Area (mVs)')
    temp = lin[lin.area > cutoff_line]   
    # Fit both exponential and log
    xdata = above_cutoff["area"]
    ydata = above_cutoff[dD_id]
    best_model, popt, sse, pcov = fit_and_select_best(xdata, ydata)
    # Generate smooth x for plotting
    x_fit = np.linspace(xdata.min(), xdata.max(), 200)
    if best_model == "linear":
        y_fit = linear_func(x_fit, *popt)
        model_label = "Linear Fit"
        parameter_text = f"y = {popt[0]}x + {popt[1]}"
    elif best_model == "decay":
        y_fit = exp_decay(x_fit, *popt)
        model_label = "Exponential Decay"
        parameter_text = f"y = {popt[0]} e^(-{popt[1]}x + {popt[2]})"
    elif best_model == "growth": # "growth"
        y_fit = exp_growth(x_fit, *popt)
        model_label="Exponential Growth"
        parameter_text = f"y = {popt[0]} (1-e^(-{popt[1]})x + {popt[2]})"
    else:
        print("Fatal error: no model converged")
    tss = np.sum((ydata - ydata.mean()) ** 2)
    if tss == 0:
        r_squared = 1.0
    else:
        r_squared = 1 - (sse / tss)
    
    # Plot the chosen best-fit curve
    plt.plot(x_fit, y_fit, 'k--', label=model_label)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), frameon=False, fancybox=False, shadow=True, ncol=2)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_path, 'Linearity.png'), bbox_inches='tight')
    plt.show()
    print(f"Chosen Model: {model_label}")
    print(f"Parameters: {parameter_text}")
    print(f" R²: {r_squared:.3f} | SSE: {sse:.3f}")

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
                axes[i].scatter(uncorrected_chain_unknown['Peak area'], uncorrected_chain_unknown['Raw dD'], label='Original dD', marker = 'x', alpha=0.6, s=60, c='k')
                
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
        plt.scatter(uncorrected_unknown['Peak area'], uncorrected_unknown['Raw dD'], label='Original '+str(label), marker = 'x', alpha=0.6, s=60, c='k')
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
        
        
        
        
        
        
        