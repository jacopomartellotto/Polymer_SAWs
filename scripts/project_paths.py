"""Canonical paths for the clean thesis simulation pipeline."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "simulation_protocol.json"
OUTPUTS = ROOT / "outputs"
DATA = OUTPUTS / "data"
SIMPLE = DATA / "simple_sampling"
PIVOT = DATA / "pivot"
PIVOT_RUNS = PIVOT / "runs"
PIVOT_SUMMARY = PIVOT / "summary.csv"
DIAGNOSTICS = OUTPUTS / "diagnostics"
THERMALIZATION = DIAGNOSTICS / "thermalization"
ANALYSIS = OUTPUTS / "analysis"
FIGURES = OUTPUTS / "figures"


def object_sha256(value) -> str:
    """Stable hash for one protocol section, independent of unrelated settings."""
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def create_output_tree() -> None:
    for path in (SIMPLE, PIVOT_RUNS, THERMALIZATION,
                 ANALYSIS, FIGURES):
        path.mkdir(parents=True, exist_ok=True)
