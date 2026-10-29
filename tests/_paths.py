"""Put scripts/ and src/ on sys.path and load 04_evaluate.py as a module."""
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
for d in (ROOT / "scripts", ROOT / "src"):
    if str(d) not in sys.path:
        sys.path.insert(0, str(d))
FIXTURES = Path(__file__).resolve().parent / "fixtures"


def load_evaluator():
    spec = importlib.util.spec_from_file_location("evaluate04", ROOT / "scripts" / "04_evaluate.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
