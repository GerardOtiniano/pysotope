from importlib import import_module

_ea_importer = import_module('...import.ea_importer', __package__)

label_element_type = _ea_importer.label_element_type
load_ea_standards = _ea_importer.load_ea_standards
add_seconds_since_start = _ea_importer.add_seconds_since_start
import_EA_data = _ea_importer.import_EA_data

__all__ = [
    'add_seconds_since_start',
    'import_EA_data',
    'label_element_type',
    'load_ea_standards',
]
