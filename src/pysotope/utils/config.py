# utils/config.py
from dataclasses import dataclass

@dataclass
class CorrectionConfig:
    drift_applied:     bool = False
    linearity_applied: bool = False

    @property
    def dD_col(self) -> str:
        # which column holds the “current” δD
        if self.linearity_applied:
            return "linearity_corrected_dD"
        if self.drift_applied:
            return "drift_corrected_dD"
        return "dD"