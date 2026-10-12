# Put bin/long_reads/ on sys.path so unit tests can import the pipeline scripts
# as modules (e.g. `import count_consensus`).

import sys
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parents[2] / "bin" / "long_reads"
sys.path.insert(0, str(BIN_DIR))
