#src/pysotope/EA/utils/
from dataclasses import dataclass

@dataclass
class CorrectionConfig:
    """
    Configuration object storing whether isotopic drift corrections have been applied.

    Attributes
    ----------
    drift_C_applied : bool, default=False
        True if drift correction for Carbon has been applied; otherwise False.
    drift_N_applied : bool, default=False
        True if drift correction for Nitrogen has been applied; otherwise False.

    Properties
    ----------
    dN_col : str
        Returns the appropriate column name for Nitrogen isotope data depending
        on whether drift correction has been applied.
        - Returns "d 15N/14N_corr" if Nitrogen drift correction applied.
        - Returns "d 15N/14N" otherwise.

    dC_col : str
        Returns the appropriate column name for Carbon isotope data depending
        on whether drift correction has been applied.
        - Returns "d 13C/12C_corr" if Carbon drift correction applied.
        - Returns "d 13C/12C" otherwise.
    """

    drift_C_applied: bool = False
    drift_N_applied: bool = False

    @property
    def dN_col(self) -> str:
        """
        Return the appropriate column name for Nitrogen isotope data depending
        on whether drift correction has been applied.

        Returns
        -------
        str
            "d 15N/14N_corr" if Nitrogen drift correction applied; otherwise "d 15N/14N".
        """

        if self.drift_N_applied:
            return "d 15N/14N_corr"
        return "d 15N/14N"

    @property
    def dC_col(self) -> str:
        """
        Return the appropriate column name for Carbon isotope data depending
        on whether drift correction has been applied.

        Returns
        -------
        str
            "d 13C/12C_corr" if Carbon drift correction applied; otherwise "d 13C/12C".
        """

        if self.drift_C_applied:
            return "d 13C/12C_corr"
        return "d 13C/12C"