# Every process needs a stub: block, or the -stub nf-test cannot run the
# workflow end to end.

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PROCESS = re.compile(r"^process\s+(\w+)\s*\{", re.MULTILINE)


def _processes():
    for nf in sorted(REPO.glob("subworkflows/**/*.nf")):
        text = nf.read_text()
        starts = [(m.group(1), m.start()) for m in PROCESS.finditer(text)]
        ends = [start for _, start in starts[1:]] + [len(text)]
        for (name, start), end in zip(starts, ends):
            yield pytest.param(text[start:end], id=f"{nf.name}:{name}")


@pytest.mark.parametrize("body", list(_processes()))
def test_process_has_stub(body):
    assert re.search(r"^\s*stub:", body, re.MULTILINE)
