# nextflow_schema.json is the single source of truth for parameters. These
# tests keep nextflow.config and the README Parameters table consistent with
# it, and make sure retired parameter names don't come back.
#
# The schema is also loaded by everit at runtime: the -stub nf-test runs with
# validate_params on, so a schema everit rejects fails there too.

import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCHEMA = json.loads((REPO / "nextflow_schema.json").read_text())

# Retired names and their migration messages are declared once, in WorkflowMain
RETIRED = re.findall(
    r"^\s+(\w+)\s*:\s*'",
    re.search(
        r"RETIRED_PARAMS = \[\n(.*?)\n\s*\]",
        (REPO / "lib" / "WorkflowMain.groovy").read_text(),
        re.DOTALL,
    ).group(1),
    re.MULTILINE,
)
GROUPS = ["UMI grouping", "Consensus", "Genotyping", "BLAZE"]
DRAFT_07_TYPES = {"string", "number", "integer", "boolean", "object", "array", "null"}
EMPTY_CELL = "—"


def _config_params():
    # Parse the literal assignments of the top-level params { } block.
    text = (REPO / "nextflow.config").read_text()
    block = re.search(r"^params\s*\{\n(.*?)^\}", text, re.MULTILINE | re.DOTALL).group(1)
    params = {}
    for line in block.splitlines():
        line = re.sub(r"(^|\s)//.*", "", line).strip()
        if not line:
            continue
        name, value = (s.strip() for s in line.split("=", 1))
        params[name] = _literal(value)
    return params


def _literal(value):
    if value == "null":
        return None
    if value in ("true", "false"):
        return value == "true"
    if value[0] in "'\"":
        assert value[-1] == value[0], f"unparsed config value: {value}"
        return value[1:-1]
    if re.fullmatch(r"-?\d+", value):
        return int(value)
    return float(value)


def _schema_params():
    return {
        name: prop
        for group in SCHEMA["definitions"].values()
        for name, prop in group["properties"].items()
    }


def _readme_table():
    # Rows of the table under "## ... Parameters": name -> (default, symbol, meaning)
    text = (REPO / "README.md").read_text()
    section = re.search(r"^## .*Parameters\n(.*?)(?=^## )", text, re.MULTILINE | re.DOTALL).group(1)
    rows = {}
    for line in section.splitlines():
        m = re.match(r"\|\s*`--(\w+)`\s*\|(.*)\|\s*$", line)
        if m:
            rows[m.group(1)] = tuple(c.strip() for c in m.group(2).split("|"))
    return rows


CONFIG = _config_params()
PARAMS = _schema_params()
VISIBLE = sorted(name for name, prop in PARAMS.items() if not prop.get("hidden"))


def test_schema_declares_draft_07():
    assert SCHEMA["$schema"].rstrip("#") == "http://json-schema.org/draft-07/schema"
    assert SCHEMA["type"] == "object"


def test_schema_refs_resolve_to_every_group():
    refs = [entry["$ref"] for entry in SCHEMA["allOf"]]
    assert sorted(refs) == sorted(f"#/definitions/{key}" for key in SCHEMA["definitions"])


@pytest.mark.parametrize("key", list(SCHEMA["definitions"]))
def test_schema_group_is_well_formed(key):
    group = SCHEMA["definitions"][key]
    assert group["type"] == "object"
    assert group["properties"]
    assert set(group.get("required", [])) <= set(group["properties"])


@pytest.mark.parametrize("name", sorted(PARAMS))
def test_schema_param_is_well_formed(name):
    prop = PARAMS[name]
    assert prop["type"] in DRAFT_07_TYPES
    assert prop["description"]
    for bound in ("minimum", "maximum", "exclusiveMinimum", "exclusiveMaximum"):
        if bound in prop:
            assert prop["type"] in ("number", "integer")
            assert type(prop[bound]) in (int, float)


@pytest.mark.parametrize("name", sorted(PARAMS))
def test_help_text_names_spec_symbol_and_units(name):
    # Generic options are pipeline plumbing, outside the model
    if name in SCHEMA["definitions"]["generic_options"]["properties"]:
        return
    help_text = PARAMS[name]["help_text"]
    assert "spec symbol" in help_text.lower() or "Spec:" in help_text
    if PARAMS[name]["type"] in ("number", "integer"):
        assert "Units:" in help_text


def test_schema_has_the_spec_groups():
    titles = [group["title"] for group in SCHEMA["definitions"].values()]
    for title in GROUPS:
        assert title in titles


def test_config_and_schema_have_equal_param_sets():
    assert sorted(CONFIG) == sorted(PARAMS)


@pytest.mark.parametrize("name", sorted(CONFIG))
def test_schema_default_equals_config(name):
    config_default = CONFIG[name]
    if config_default in (None, ""):
        assert "default" not in PARAMS.get(name, {})
    else:
        assert PARAMS[name]["default"] == config_default


@pytest.mark.parametrize("name", sorted(CONFIG))
def test_config_default_is_within_schema_bounds(name):
    prop, value = PARAMS[name], CONFIG[name]
    if value is None or prop["type"] not in ("number", "integer"):
        return
    assert value >= prop.get("minimum", value)
    assert value <= prop.get("maximum", value)
    assert value > prop.get("exclusiveMinimum", value - 1)
    assert value < prop.get("exclusiveMaximum", value + 1)


@pytest.mark.parametrize(
    "name, bounds",
    [
        ("min_fraction", {"exclusiveMinimum": 0.5, "maximum": 1}),
        ("error_floor", {"exclusiveMinimum": 0, "exclusiveMaximum": 1}),
        ("min_reads", {"minimum": 1}),
        ("snv_min_base_quality", {"minimum": 0}),
        ("min_mapping_quality", {"minimum": 0}),
        ("blaze_expect_cells", {"minimum": 1}),
        ("blaze_edit_distance", {"minimum": 0}),
    ],
)
def test_model_param_bounds(name, bounds):
    for bound, value in bounds.items():
        assert PARAMS[name][bound] == value


def test_min_fraction_default_is_0_6():
    assert CONFIG["min_fraction"] == 0.6


def test_readme_table_lists_the_visible_params():
    assert sorted(_readme_table()) == VISIBLE


@pytest.mark.parametrize("name", VISIBLE)
def test_readme_row_matches_schema(name):
    default, symbol, meaning = _readme_table()[name]
    prop = PARAMS[name]
    expected = f"`{prop['default']}`" if "default" in prop else EMPTY_CELL
    assert default == expected
    if symbol != EMPTY_CELL:
        assert symbol in prop["help_text"]
    assert meaning


def test_retired_names_are_read_from_workflow_main():
    assert RETIRED == ["window_width", "min_dist_umi_collapsing", "baseline_error", "min_base_quality"]


@pytest.mark.parametrize(
    "path",
    ["nextflow.config", "nextflow_schema.json", "main.nf", "subworkflows", "bin"],
)
def test_no_retired_names_remain(path):
    root = REPO / path
    files = [root] if root.is_file() else [
        f for f in root.rglob("*") if f.is_file() and "__pycache__" not in f.parts
    ]
    for f in files:
        for name in RETIRED:
            assert not re.search(rf"\b{name}\b", f.read_text()), f"{name} in {f}"


def test_no_retired_names_in_readme_table():
    assert not set(_readme_table()) & set(RETIRED)
