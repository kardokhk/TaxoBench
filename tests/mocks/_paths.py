"""Put the mock builders (mocks/scripts/ and its lib/) on sys.path."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
for d in (ROOT / "mocks" / "scripts", ROOT / "mocks" / "scripts" / "lib"):
    if str(d) not in sys.path:
        sys.path.insert(0, str(d))
