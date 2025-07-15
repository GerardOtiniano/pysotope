#src/pysotope/EA/utils/
import matplotlib.pyplot as plt

def apply_drift_correction_plot(df):
    """
    Filter standards with name 'SORGHUM' and plot Rt vs isotope values
    for both Nitrogen and Carbon.

    Parameters:
        df (pd.DataFrame): Full EA dataset with labeled element type.
    """
    # Filter rows with Identifier 1 = "SORGHUM" (case insensitive)
    sorghum = df[df['Identifier 1'].str.lower() == 'sorghum'].copy()

    if sorghum.empty:
        print("No 'SORGHUM' standards found.")
        return

    # Plot for Nitrogen
    sorghum_N = sorghum[sorghum['Element Type'] == 'Nitrogen']
    # samp_n = df[df['Element Type']=='Nitrogen']
    if not sorghum_N.empty:
        # plt.figure(figsize=(8,4))
        # # plt.scatter(samp_n['Rt'], samp_n['d 15N/14N'], alpha=0.6, c = 'k', marker='x')
        # plt.scatter(sorghum_N['Rt'], sorghum_N['d 15N/14N'], c=sorghum_N['Seconds Since Start'],
        #             ec = 'k', label='SORGHUM - N', s=200, alpha=0.6)
        # plt.xlabel("Rt")
        # plt.ylabel("δ15N/14N")
        # plt.legend()
        # plt.show()
        
        plt.figure(figsize=(8,4))
        # # plt.scatter(samp_n['Rt'], samp_n['d 15N/14N'], alpha=0.6, c = 'k', marker='x')
        # plt.scatter(sorghum_N['Rt'], sorghum_N['d 15N/14N'], c=sorghum_N['Area All'],
        #             ec = 'k', label='SORGHUM - N', s=200, alpha=0.6)
        # plt.xlabel("Rt")
        # plt.ylabel("δ15N/14N")
        # plt.legend()
        # plt.show()
        for x in ["Area All", "Ampl. 28", "Ampl. 29", "Area 28", "Area 29"]:
            plt.figure(figsize=(8,4))
            # plt.scatter(samp_n['Rt'], samp_n['d 15N/14N'], alpha=0.6, c = 'k', marker='x')
            plt.scatter(sorghum_N['Rt'], sorghum_N['d 15N/14N'], c=sorghum_N[x],
                        ec = 'k', label=f'SORGHUM - N {x}', s=200, alpha=0.6)
            plt.xlabel("Rt")
            plt.ylabel("δ15N/14N")
            plt.legend()
            plt.show()

    # # Plot for Carbon
    # sorghum_C = sorghum[sorghum['Element Type'] == 'Carbon']
    # # samp_n = df[df['Element Type']=='Carbon']
    # if not sorghum_C.empty:
    #     plt.figure(figsize=(8,4))
    #     # plt.scatter(samp_n['Rt'], samp_n['d 15N/14N'], alpha=0.6, c = 'k', marker='x')
    #     plt.scatter(sorghum_C['Rt'], sorghum_C['d 13C/12C'], c=sorghum_C['Seconds Since Start'],
    #                 ec = 'k', label='SORGHUM - C')
    #     plt.xlabel("Rt")
    #     plt.ylabel("δ13C/12C")
    #     plt.legend()
    #     plt.show()

    return sorghum