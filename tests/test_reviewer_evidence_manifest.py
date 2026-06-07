import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_reviewer_evidence_manifest.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("build_reviewer_evidence_manifest", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_reviewer_evidence_manifest_contains_core_release_artifacts(tmp_path):
    mod = _load_module()
    manifest = mod.build_manifest()
    outputs = mod.write_manifest(manifest, tmp_path)

    assert manifest["status"] == "complete"
    paths = {record["path"] for record in manifest["files"]}
    assert "docs/REVIEWER_GUIDE.md" in paths
    assert "data/benchmarks/real_corpus_mpc_500/corpus_manifest.json" in paths
    assert "data/benchmarks/closure/closure_report.json" in paths
    assert manifest["certificate_hash"]

    written = json.loads(Path(outputs["json"]).read_text(encoding="utf-8"))
    assert written["certificate_hash"] == manifest["certificate_hash"]
    assert Path(outputs["markdown"]).exists()
