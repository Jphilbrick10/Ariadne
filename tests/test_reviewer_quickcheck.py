import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "reviewer_quickcheck.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("reviewer_quickcheck", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_reviewer_quickcheck_dry_run_writes_machine_readable_report(tmp_path):
    mod = _load_module()
    report = tmp_path / "reviewer_report.json"

    rc = mod.main(["--profile", "smoke", "--dry-run", "--report", str(report)])

    assert rc == 0
    data = json.loads(report.read_text(encoding="utf-8"))
    assert data["status"] == "pass"
    assert data["dry_run"] is True
    assert [gate["name"] for gate in data["gates"]] == ["import_smoke", "stage24_packaging"]
    assert all(gate["status"] == "dry_run" for gate in data["gates"])
