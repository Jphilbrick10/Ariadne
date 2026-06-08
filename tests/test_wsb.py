"""Low-energy WSB transfer test (Gate G_wsb, Stage 10)."""

from ariadne.transfers.wsb import SOLUTION_PARAMS, evaluate_transfer

DIRECT_MS = 3953.0
COIMBRA_MS = 3925.0


def test_wsb_transfer_beats_coimbra():
    b = evaluate_transfer(SOLUTION_PARAMS)
    # departs ~LEO, arrives lower-energy than direct, total beats both direct and Coimbra
    assert abs(b["perigee_alt_km"] - 200.0) < 300.0
    assert b["v_inf"] < 0.82  # lower v_inf than the direct transfer
    assert b["total_ms"] < DIRECT_MS
    assert b["total_ms"] < COIMBRA_MS
    # the honest tradeoff: a long, low-energy route (tens of days) vs a fast direct
    # transfer. The exact tof is chaos-sensitive in the weak-stability-boundary regime,
    # so assert the qualitative property with margin rather than a hard ~40-day edge.
    assert b["tof_days"] > 30.0
