import json
from pathlib import Path

CONFIG_PATH = Path(__file__).resolve().parent.parent / "chains" / "chains.json"

def load_chain_config():
    with open(CONFIG_PATH, "r") as f:
        return json.load(f)

def get_selected_chains():
    cfg = load_chain_config()
    selected_key = cfg["selected"]
    return cfg["sets"][selected_key]

def save_chain_config(cfg):
    with open(CONFIG_PATH, "w") as f:
        json.dump(cfg, f, indent=4)