"""Load pooled measurements and check their sampling scheme against the protocol."""
import json
import pandas as pd
from project_paths import CONFIG, SIMPLE


def load_grid():
    refined = SIMPLE / "refined" / "grid_refined.csv"
    data = pd.read_csv(refined if refined.exists() else SIMPLE / "base.csv")
    expected = json.loads(CONFIG.read_text(encoding="utf-8"))["simple_sampling"]["proposal_scheme"]
    if ("proposal_scheme" not in data
            or not data.proposal_scheme.eq(expected).all()):
        raise RuntimeError("Simple-sampling data do not match the configured proposal scheme")
    return data
