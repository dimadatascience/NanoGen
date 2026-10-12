# Put bin/long_reads/ on sys.path so unit tests can import the pipeline scripts
# as modules (e.g. `import count_consensus`), and tests/fixtures/ for the
# synthetic fixture generator (`import make_fixtures`).

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BIN_DIR = REPO / "bin" / "long_reads"
FIXTURES_DIR = REPO / "tests" / "fixtures"
sys.path.insert(0, str(BIN_DIR))
sys.path.insert(0, str(FIXTURES_DIR))
