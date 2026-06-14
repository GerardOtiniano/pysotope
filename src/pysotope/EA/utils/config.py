#src/pysotope/EA/utils/
from dataclasses import dataclass

@dataclass
class CorrectionConfig:
    """
    Configuration object storing which EA isotope corrections have been applied.

    Attributes
    ----------
    drift_C_applied : bool, default=False
        True if drift correction for Carbon has been applied; otherwise False.
    drift_N_applied : bool, default=False
        True if drift correction for Nitrogen has been applied; otherwise False.
    linearity_C_applied : bool, default=False
        True if linearity correction for Carbon has been applied; otherwise False.
    linearity_N_applied : bool, default=False
        True if linearity correction for Nitrogen has been applied; otherwise False.

    Properties
    ----------
    dN_col : str
        Returns the current Nitrogen isotope-value column and a description.
        The order is linearity-corrected, drift-corrected, then raw.

    dC_col : str
        Returns the current Carbon isotope-value column and a description.
        The order is linearity-corrected, drift-corrected, then raw.
    """

    drift_C_applied: bool = False
    drift_N_applied: bool = False
    linearity_C_applied: bool = False
    linearity_N_applied: bool = False

    @property
    def dN_col(self) -> str:
        """
        Return the appropriate column name for Nitrogen isotope data depending
        on which correction has been applied.

        Returns
        -------
        tuple[str, str]
            The current Nitrogen isotope-value column and a short description.
        """

        if self.linearity_N_applied:
            return "d 15N/14N_lin_corr","Linearity corrected"
        if self.drift_N_applied:
            return "d 15N/14N_corr", "Drift corrected"
        return "d 15N/14N", "Raw"

    @property
    def dC_col(self) -> str:
        """
        Return the appropriate column name for Carbon isotope data depending
        on which correction has been applied.

        Returns
        -------
        tuple[str, str]
            The current Carbon isotope-value column and a short description.
        """

        if self.linearity_C_applied:
            return "d 13C/12C_lin_corr","Linearity corrected"
        if self.drift_C_applied:
            return "d 13C/12C_corr", "Drift corrected"
        return "d 13C/12C", "Raw"
