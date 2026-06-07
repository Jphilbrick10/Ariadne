"""Build the frozen reviewer evidence manifest.

This manifest is deliberately smaller and more reviewer-facing than the full
benchmark artifact manifest. It hashes the public claims map, launch audit,
closure ledger, frozen labelled corpus artifacts, and key benchmark summaries
that a scientist is likely to inspect first.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "data" / "benchmarks" / "reviewer_evidence"


@dataclass(frozen=True)
class EvidenceFile:
    path: str
    role: str
    claim: str


EVIDENCE_FILES = [
    EvidenceFile("README.md", "public_entrypoint", "Top-level project claims and honest scope."),
    EvidenceFile("docs/REVIEWER_GUIDE.md", "review_map", "Claim-to-evidence review path."),
    EvidenceFile("docs/PUBLIC_LAUNCH_AUDIT.md", "release_audit", "Commands run for public launch."),
    EvidenceFile(
        "docs/WHITE_PAPER.md", "scientific_narrative", "Methods, validation, and limitations."
    ),
    EvidenceFile(
        "docs/VALIDATION_RESULTS.md", "validation_report", "Discovery pipeline validation results."
    ),
    EvidenceFile("docs/REAL_BENCHMARKS.md", "real_data_report", "Real-data benchmark summary."),
    EvidenceFile(
        "data/benchmarks/closure/closure_report.json", "closure_ledger", "Readiness ledger."
    ),
    EvidenceFile(
        "data/benchmarks/closure/closure_report.md",
        "closure_summary",
        "Human-readable closure summary.",
    ),
    EvidenceFile(
        "data/benchmarks/real_corpus_mpc_500/corpus_manifest.json",
        "frozen_labelled_corpus_manifest",
        "500-case MPC-derived labelled corpus manifest.",
    ),
    EvidenceFile(
        "data/benchmarks/real_corpus_mpc_500/labelled_cases.jsonl",
        "frozen_labelled_cases",
        "Frozen labelled cases used by the real-corpus benchmark.",
    ),
    EvidenceFile(
        "data/benchmarks/real_corpus_mpc_500/benchmark/metrics.json",
        "real_corpus_metrics",
        "Accuracy/calibration metrics for the frozen labelled corpus.",
    ),
    EvidenceFile(
        "data/benchmarks/real_corpus_mpc_500/benchmark/reliability_curve.csv",
        "real_corpus_reliability",
        "Reliability curve for calibration review.",
    ),
    EvidenceFile(
        "data/benchmarks/external_inference/metrics.json",
        "external_inference_metrics",
        "External inference benchmark metrics.",
    ),
    EvidenceFile(
        "data/benchmarks/inference_adversarial_builtin/metrics.json",
        "adversarial_metrics",
        "Adversarial false-positive benchmark metrics.",
    ),
    EvidenceFile(
        "data/benchmarks/artifact_integrity/artifact_manifest.json",
        "artifact_integrity_manifest",
        "Full benchmark artifact hash manifest.",
    ),
    EvidenceFile("CITATION.cff", "citation_metadata", "Citation and release metadata."),
    EvidenceFile("LICENSE", "license", "PolyForm Noncommercial license text."),
]


REPRO_COMMANDS = [
    {
        "name": "reviewer_quickcheck",
        "command": "PYTHONPATH=src python scripts/reviewer_quickcheck.py",
        "purpose": "Run production lint, fast tests, strict docs, and closure ledger.",
    },
    {
        "name": "closure_ledger",
        "command": "python scripts/build_closure_report.py --fail-on-critical",
        "purpose": "Regenerate closure evidence and fail on critical residuals.",
    },
    {
        "name": "stage24_packaging",
        "command": "PYTHONPATH=src python -m ariadne.validate.stage24",
        "purpose": "Verify installability, license metadata, CI, and core imports.",
    },
    {
        "name": "stage32_selector",
        "command": "PYTHONPATH=src python -m ariadne.validate.stage32",
        "purpose": "Verify ensemble backend faithfulness and selector calibration.",
    },
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def text_line_count(path: Path) -> int | None:
    try:
        return len(path.read_text(encoding="utf-8").splitlines())
    except UnicodeDecodeError:
        return None


def json_summary(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        return {}
    keys = [
        "schema",
        "status",
        "readiness_score",
        "critical_failures",
        "certificate_hash",
        "n_cases",
        "case_hash",
        "accuracy",
        "macro_f1",
        "safe_accuracy",
        "ece",
    ]
    return {k: data[k] for k in keys if k in data}


def build_manifest() -> dict:
    files = []
    missing = []
    for item in EVIDENCE_FILES:
        path = ROOT / item.path
        if not path.exists():
            missing.append(item.path)
            continue
        record = {
            "path": item.path,
            "role": item.role,
            "claim": item.claim,
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        line_count = text_line_count(path)
        if line_count is not None:
            record["line_count"] = line_count
        if path.suffix.lower() == ".json":
            summary = json_summary(path)
            if summary:
                record["summary"] = summary
        files.append(record)

    certificate_payload = json.dumps(
        [{"path": f["path"], "sha256": f["sha256"]} for f in files],
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    certificate_hash = hashlib.sha256(certificate_payload).hexdigest()
    return {
        "schema": "ariadne.reviewer_evidence_manifest.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete" if not missing else "partial",
        "missing": missing,
        "file_count": len(files),
        "certificate_hash": certificate_hash,
        "repro_commands": REPRO_COMMANDS,
        "files": files,
        "limits": [
            "Source-available under PolyForm Noncommercial, not OSI open-source.",
            "Frozen MPC corpus is useful evidence, not a blind NOIRLab/ZTF/LSST competition.",
            "Network-backed external benchmarks should be rerun before new scientific claims.",
            "Operational mission design still requires certified workflows such as GMAT/Monte/Copernicus.",
        ],
    }


def write_markdown(manifest: dict, path: Path) -> None:
    lines = [
        "# Reviewer Evidence Manifest",
        "",
        f"- status: `{manifest['status']}`",
        f"- file count: `{manifest['file_count']}`",
        f"- certificate hash: `{manifest['certificate_hash']}`",
        "",
        "## Reproduction Commands",
        "",
        "| Name | Command | Purpose |",
        "|---|---|---|",
    ]
    for cmd in manifest["repro_commands"]:
        lines.append(f"| `{cmd['name']}` | `{cmd['command']}` | {cmd['purpose']} |")
    lines.extend(
        [
            "",
            "## Frozen Evidence Files",
            "",
            "| Role | Path | Bytes | SHA-256 |",
            "|---|---|---:|---|",
        ]
    )
    for rec in manifest["files"]:
        lines.append(f"| `{rec['role']}` | `{rec['path']}` | {rec['bytes']} | `{rec['sha256']}` |")
    lines.extend(["", "## Limits", ""])
    for limit in manifest["limits"]:
        lines.append(f"- {limit}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_manifest(manifest: dict, out_dir: Path) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "reviewer_evidence_manifest.json"
    md_path = out_dir / "reviewer_evidence_manifest.md"
    json_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(manifest, md_path)
    return {"json": str(json_path), "markdown": str(md_path)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    manifest = build_manifest()
    outputs = write_manifest(manifest, args.out_dir)
    print(f"status={manifest['status']}")
    print(f"file_count={manifest['file_count']}")
    print(f"certificate_hash={manifest['certificate_hash']}")
    print(f"json={outputs['json']}")
    print(f"markdown={outputs['markdown']}")
    return 0 if manifest["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
