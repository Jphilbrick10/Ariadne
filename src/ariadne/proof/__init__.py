"""System-wide proof, closure, and residual intelligence for Ariadne."""

from .artifact_manifest import build_artifact_manifest, write_artifact_manifest
from .closure import (
    ArtifactEvidence,
    ClosureLedger,
    ClosureReport,
    GateResult,
    ResidualSignal,
    SubsystemContract,
    build_closure_report,
    load_json_artifact,
    stable_hash,
    write_closure_report,
)
from .defaults import (
    audit_png_directory,
    build_default_ariadne_closure,
    collect_default_evidence,
    default_contracts,
    known_residuals,
)
from .dream import DreamExperiment, DreamRun, build_dream_run, write_dream_run
from .high_fidelity import (
    covariance_envelope_evidence,
    independent_crosscheck_evidence,
    nbody_replay_evidence,
)
from .promotion import (
    PromotionEvidence,
    PromotionReport,
    PromotionRung,
    PromotionThresholds,
    RoutePromotionCertificate,
    load_routes_from_navigator_report,
    promote_route,
    promote_routes,
    write_promotion_report,
)
from .visuals import navigator_visual_contract_evidence

__all__ = [
    "ArtifactEvidence",
    "ClosureLedger",
    "ClosureReport",
    "DreamExperiment",
    "DreamRun",
    "GateResult",
    "PromotionEvidence",
    "PromotionReport",
    "PromotionRung",
    "PromotionThresholds",
    "ResidualSignal",
    "RoutePromotionCertificate",
    "SubsystemContract",
    "audit_png_directory",
    "build_artifact_manifest",
    "build_closure_report",
    "build_default_ariadne_closure",
    "build_dream_run",
    "collect_default_evidence",
    "covariance_envelope_evidence",
    "default_contracts",
    "independent_crosscheck_evidence",
    "known_residuals",
    "load_json_artifact",
    "load_routes_from_navigator_report",
    "navigator_visual_contract_evidence",
    "nbody_replay_evidence",
    "promote_route",
    "promote_routes",
    "stable_hash",
    "write_artifact_manifest",
    "write_closure_report",
    "write_dream_run",
    "write_promotion_report",
]
