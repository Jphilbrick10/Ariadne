"""One-command reviewer evidence harness for Ariadne.

The public docs tell reviewers what to run; this script makes that path
executable and records the result as JSON. It intentionally runs only local,
read-only gates by default: no network catalog fetches and no writes outside the
chosen report path.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = ROOT / "results" / "reviewer_quickcheck_report.json"
TAIL_CHARS = 6000


@dataclass(frozen=True)
class Gate:
    name: str
    command: tuple[str, ...]
    proof: str


@dataclass
class GateResult:
    name: str
    command: list[str]
    proof: str
    returncode: int | None
    duration_s: float
    status: str
    stdout_tail: str = ""
    stderr_tail: str = ""


def _py(*args: str) -> tuple[str, ...]:
    return (sys.executable, *args)


def quick_gates() -> list[Gate]:
    """Local gates a scientific reviewer can run without network credentials."""
    return [
        Gate(
            "production_lint_sanity",
            _py("-m", "ruff", "check", "src/ariadne", "--select", "E9,F63,F7,F82,B023,B904"),
            "fatal Python diagnostics, loop-closure bugs, and exception-chain correctness",
        ),
        Gate(
            "fast_offline_tests",
            _py("-m", "pytest", "-m", "not slow", "-q"),
            "offline unit/integration behavior across the maintained test suite",
        ),
        Gate(
            "sphinx_docs_strict",
            _py(
                "-m",
                "sphinx",
                "-b",
                "html",
                "-W",
                "--keep-going",
                "docs/sphinx",
                "docs/sphinx/_build/html",
            ),
            "public docs build with warnings treated as errors",
        ),
        Gate(
            "closure_ledger",
            _py("scripts/build_closure_report.py", "--fail-on-critical"),
            "proof ledger reports complete readiness and zero critical failures",
        ),
    ]


def release_gates() -> list[Gate]:
    """Additional packaging checks for release reviewers."""
    wheel = "dist/ariadne_astro-1.0.0rc2-py3-none-any.whl"
    return [
        Gate(
            "package_build",
            _py("-m", "build"),
            "source distribution and wheel can be rebuilt from the checkout",
        ),
        Gate(
            "distribution_metadata",
            _py("-m", "twine", "check", "dist/*"),
            "built distributions pass metadata validation",
        ),
        Gate(
            "wheel_install_smoke",
            _py(
                "-c",
                (
                    "import subprocess, sys; "
                    f"subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--force-reinstall', '--no-deps', {wheel!r}]); "
                    "import ariadne; "
                    "assert ariadne.__version__ == '1.0.0rc2'; "
                    "assert ariadne.system('EARTH_MOON').mu > 0; "
                    "print('wheel smoke OK', ariadne.__version__)"
                ),
            ),
            "freshly built wheel installs and imports with the top-level API available",
        ),
    ]


def gates_for_profile(profile: str) -> list[Gate]:
    if profile == "quick":
        return quick_gates()
    if profile == "release":
        return quick_gates() + release_gates()
    if profile == "smoke":
        return [
            Gate(
                "import_smoke",
                _py("-c", "import ariadne; print('ariadne', ariadne.__version__)"),
                "package imports from the current environment",
            ),
            Gate(
                "stage24_packaging",
                _py("-m", "ariadne.validate.stage24"),
                "packaging/source-availability validation stage passes",
            ),
        ]
    raise ValueError(f"unknown profile: {profile}")


def _tail(text: str) -> str:
    return text[-TAIL_CHARS:] if len(text) > TAIL_CHARS else text


def run_gate(gate: Gate, *, env: dict[str, str], dry_run: bool) -> GateResult:
    start = time.perf_counter()
    if dry_run:
        return GateResult(
            name=gate.name,
            command=list(gate.command),
            proof=gate.proof,
            returncode=None,
            duration_s=0.0,
            status="dry_run",
        )

    proc = subprocess.run(
        gate.command,
        cwd=ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    duration = time.perf_counter() - start
    return GateResult(
        name=gate.name,
        command=list(gate.command),
        proof=gate.proof,
        returncode=proc.returncode,
        duration_s=round(duration, 3),
        status="pass" if proc.returncode == 0 else "fail",
        stdout_tail=_tail(proc.stdout),
        stderr_tail=_tail(proc.stderr),
    )


def build_report(profile: str, results: Iterable[GateResult], *, dry_run: bool) -> dict:
    rows = [asdict(r) for r in results]
    failed = [r["name"] for r in rows if r["status"] == "fail"]
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "profile": profile,
        "dry_run": dry_run,
        "status": "pass" if not failed else "fail",
        "failed_gates": failed,
        "python": sys.version,
        "platform": platform.platform(),
        "repo_root": str(ROOT),
        "gates": rows,
    }


def write_report(report: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile",
        choices=("quick", "release", "smoke"),
        default="quick",
        help="Gate set to run. quick is local scientific review; release adds package checks.",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=DEFAULT_REPORT,
        help="JSON report path. Defaults to ignored results/reviewer_quickcheck_report.json.",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Write planned gates without running commands."
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    env = os.environ.copy()
    src = str(ROOT / "src")
    env["PYTHONPATH"] = src + os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else src

    results = []
    for gate in gates_for_profile(args.profile):
        print(f"[{gate.name}] {' '.join(gate.command)}", flush=True)
        result = run_gate(gate, env=env, dry_run=args.dry_run)
        print(f"  -> {result.status} ({result.duration_s:.3f}s)", flush=True)
        results.append(result)

    report = build_report(args.profile, results, dry_run=args.dry_run)
    write_report(report, args.report)
    print(f"report={args.report}")
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
