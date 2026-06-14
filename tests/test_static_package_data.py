import csv
import subprocess
import sys
import tomllib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class StaticPackageDataTests(unittest.TestCase):
    def test_pyproject_includes_runtime_data_files(self):
        pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text())
        includes = set(pyproject["tool"]["poetry"].get("include", []))

        self.assertIn("src/pysotope/standards_manager/data/*.csv", includes)
        self.assertIn("src/pysotope/utils/chains/chains.json", includes)

    def test_ea_standard_files_have_required_selector_columns(self):
        for filename in ("ea_RS_dN.csv", "ea_RS_dC.csv"):
            with (ROOT / "src" / "pysotope" / "standards_manager" / "data" / filename).open(
                newline=""
            ) as handle:
                header = next(csv.reader(handle))

            self.assertIn("ID", header)
            self.assertIn("Linearity Standard", header)
            self.assertIn("RS accuracy check", header)
            self.assertIn("drift", header)

    def test_gcirms_standard_files_have_required_selector_columns(self):
        for filename in ("gcirms_RS_dD.csv", "gcirms_RS_dC.csv"):
            with (ROOT / "src" / "pysotope" / "standards_manager" / "data" / filename).open(
                newline=""
            ) as handle:
                header = next(csv.reader(handle))

            self.assertIn("ID", header)
            self.assertIn("type", header)
            self.assertIn("chain length", header)
            self.assertIn("RS accuracy check", header)
            self.assertIn("Use as Standard", header)

    def test_standard_boolean_parsers_accept_csv_round_trips(self):
        sources = [
            ROOT / "src" / "pysotope" / "standards_manager" / "editor.py",
            ROOT / "src" / "pysotope" / "corrections" / "drift_correction.py",
            ROOT / "src" / "pysotope" / "corrections" / "linearity_correction.py",
            ROOT / "src" / "pysotope" / "EA" / "utils" / "VPDB_correction.py",
        ]
        for source in sources:
            with self.subTest(source=source):
                text = source.read_text()
                self.assertIn('"1.0"', text)

    def test_gcirms_standards_are_coerced_before_filtering(self):
        source = ROOT / "src" / "pysotope" / "utils" / "base_functions.py"
        text = source.read_text()

        coerce_pos = text.index('df[col] = _coerce_bool_series(df[col])')
        filter_pos = text.index("df = df[df['Use as Standard']].copy()")
        self.assertLess(coerce_pos, filter_pos)

    def test_ea_element_classification_prefers_component(self):
        source = ROOT / "src" / "pysotope" / "import" / "ea_importer.py"
        text = source.read_text()

        component_pos = text.index('component == "N2"')
        fallback_pos = text.index("unlabeled & out['d 15N/14N'].notna()")
        self.assertLess(component_pos, fallback_pos)
        self.assertIn('component == "CO2"', text)

    def test_ea_drift_diagnostics_are_built_before_user_confirmation(self):
        source = ROOT / "src" / "pysotope" / "EA" / "eaAnalyze.py"
        text = source.read_text()

        prepare_n_pos = text.index('drift_n = prepare_drift_diagnostic("N"')
        prepare_c_pos = text.index('drift_c = prepare_drift_diagnostic("C"')
        confirm_pos = text.index("apply_n, apply_c = drift_confirm")
        apply_helper_pos = text.index("def apply_selected_drift")

        self.assertLess(prepare_n_pos, confirm_pos)
        self.assertLess(prepare_c_pos, confirm_pos)
        self.assertLess(confirm_pos, apply_helper_pos)

        prepare_block = text[text.index("def prepare_drift_diagnostic"):confirm_pos]
        self.assertIn("get_ea_drift_standards(", prepare_block)
        self.assertIn("build_drift_model(", prepare_block)
        self.assertIn("plot_drift_diagnostics(", prepare_block)

    def test_ea_linearity_uses_cfg_chain_before_reference_calibration(self):
        source = ROOT / "src" / "pysotope" / "EA" / "eaAnalyze.py"
        text = source.read_text()

        linearity_pos = text.index("# Linearity Correction")
        reference_pos = text.index("# Reference standard calibration")
        self.assertLess(linearity_pos, reference_pos)

        linearity_block = text[linearity_pos:reference_pos]
        self.assertIn('input_col, input_desc = cfg.dN_col if tag == "N" else cfg.dC_col', linearity_block)
        self.assertIn("target_column=input_col", linearity_block)
        self.assertIn("input_column=input_col", linearity_block)
        self.assertIn("cfg.linearity_N_applied = True", linearity_block)
        self.assertIn("cfg.linearity_C_applied = True", linearity_block)

        n_apply_pos = linearity_block.index('apply_linearity_model(')
        n_cfg_pos = linearity_block.index("cfg.linearity_N_applied = True")
        self.assertLess(n_apply_pos, n_cfg_pos)

    def test_ea_config_prefers_linearity_then_drift_then_raw(self):
        source = ROOT / "src" / "pysotope" / "EA" / "utils" / "config.py"
        text = source.read_text()

        n_linearity_pos = text.index("if self.linearity_N_applied:")
        n_drift_pos = text.index("if self.drift_N_applied:")
        n_raw_pos = text.index('return "d 15N/14N", "Raw"')
        self.assertLess(n_linearity_pos, n_drift_pos)
        self.assertLess(n_drift_pos, n_raw_pos)

        c_linearity_pos = text.index("if self.linearity_C_applied:")
        c_drift_pos = text.index("if self.drift_C_applied:")
        c_raw_pos = text.index('return "d 13C/12C", "Raw"')
        self.assertLess(c_linearity_pos, c_drift_pos)
        self.assertLess(c_drift_pos, c_raw_pos)

    @unittest.skipUnless((ROOT / ".git").exists(), "requires a git checkout")
    def test_no_python_cache_files_are_tracked(self):
        result = subprocess.run(
            ["git", "-C", str(ROOT), "ls-files", "*__pycache__*", "*.pyc"],
            check=True,
            capture_output=True,
            text=True,
        )

        self.assertEqual("", result.stdout.strip())


if __name__ == "__main__":
    unittest.main()
