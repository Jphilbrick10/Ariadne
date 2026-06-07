"""Interactive 3D map of the solar system with the candidate outer perturber placed
where the data says it should be.

Shows, in the heliocentric ecliptic frame (AU):
  * the Sun and the 8 planets (real J2000 elements; current positions by solving
    Kepler from J2000 mean anomalies to a given date)
  * the real extreme/detached TNOs from JPL (the actual anomaly a perturber would
    explain), with the famous ones labelled
  * "Planet Nine": the literature best-estimate orbit, with its most-likely current
    location marked near aphelion -- and, honestly, the large uncertainty on that

  python scripts/map_solar_system_3d.py
Outputs a self-contained HTML you can rotate/zoom in any browser.
"""

from __future__ import annotations

import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from ariadne.discovery.frontier.tno_clustering import extreme_population, fetch_tnos  # noqa: E402

# planet J2000 elements: a(AU), e, i, Omega, omega (deg), M0 (deg), period (yr), color
PLANETS = {
    "Mercury": (0.387, 0.2056, 7.00, 48.33, 29.12, 174.79, 0.2408, "#b1b1b1"),
    "Venus": (0.723, 0.0068, 3.39, 76.68, 54.88, 50.45, 0.6152, "#e8cda2"),
    "Earth": (1.000, 0.0167, 0.00, 348.74, 114.21, 357.52, 1.0000, "#4f8ff7"),
    "Mars": (1.524, 0.0934, 1.85, 49.56, 286.50, 19.41, 1.8808, "#e2733b"),
    "Jupiter": (5.203, 0.0484, 1.30, 100.49, 273.87, 19.65, 11.862, "#d9a066"),
    "Saturn": (9.537, 0.0539, 2.49, 113.64, 339.39, 317.51, 29.457, "#e6cd8a"),
    "Uranus": (19.19, 0.0473, 0.77, 74.00, 96.99, 142.27, 84.011, "#9fe3e3"),
    "Neptune": (30.07, 0.0086, 1.77, 131.78, 276.34, 256.23, 164.79, "#5b7cff"),
}
J2000 = datetime(2000, 1, 1, 12, tzinfo=timezone.utc)


def _rot(a, e, i, Om, om, nu):
    """Perifocal->ecliptic xyz (AU) at true anomaly nu (rad). Arrays ok."""
    i, Om, om = map(math.radians, (i, Om, om))
    r = a * (1 - e**2) / (1 + e * np.cos(nu))
    xp, yp = r * np.cos(nu), r * np.sin(nu)
    cO, sO, ci, si, co, so = (
        math.cos(Om),
        math.sin(Om),
        math.cos(i),
        math.sin(i),
        math.cos(om),
        math.sin(om),
    )
    x = (cO * co - sO * so * ci) * xp + (-cO * so - sO * co * ci) * yp
    y = (sO * co + cO * so * ci) * xp + (-sO * so + cO * co * ci) * yp
    z = (so * si) * xp + (co * si) * yp
    return x, y, z


def orbit(a, e, i, Om, om, n=400):
    nu = np.linspace(0, 2 * math.pi, n)
    return _rot(a, e, i, Om, om, nu)


def _kepler(M, e):
    M = math.radians(M % 360)
    E = M
    for _ in range(50):
        E -= (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
    return 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2), math.sqrt(1 - e) * math.cos(E / 2))


def planet_now(a, e, i, Om, om, M0, period, when):
    yrs = (when - J2000).total_seconds() / (365.25 * 86400)
    nu = _kepler(M0 + 360 * yrs / period, e)
    return _rot(a, e, i, Om, om, np.array([nu]))


def main():
    import plotly.graph_objects as go

    when = datetime.now(timezone.utc)
    fig = go.Figure()

    # Sun
    fig.add_trace(
        go.Scatter3d(
            x=[0],
            y=[0],
            z=[0],
            mode="markers+text",
            marker=dict(size=7, color="#ffd23f"),
            text=["Sun"],
            textposition="top center",
            name="Sun",
            textfont=dict(color="#ffd23f"),
        )
    )

    # planets: orbit + current position
    for name, (a, e, i, Om, om, M0, per, col) in PLANETS.items():
        ox, oy, oz = orbit(a, e, i, Om, om)
        fig.add_trace(
            go.Scatter3d(
                x=ox,
                y=oy,
                z=oz,
                mode="lines",
                line=dict(color=col, width=2),
                name=name,
                showlegend=False,
                hoverinfo="skip",
            )
        )
        px, py, pz = planet_now(a, e, i, Om, om, M0, per, when)
        fig.add_trace(
            go.Scatter3d(
                x=px,
                y=py,
                z=pz,
                mode="markers+text",
                marker=dict(size=4, color=col),
                text=[name],
                textposition="top center",
                name=name,
                textfont=dict(color=col, size=9),
            )
        )

    # real detached / extreme TNOs (the anomaly)
    tn = extreme_population(fetch_tnos(), a_min=150, q_min=30)
    famous = {
        "sedna": "Sedna",
        "leleakuhonua": "Leleakuhonua (Goblin)",
        "2012 vp113": "2012 VP113 (Biden)",
        "2015 tg387": "Leleakuhonua",
        "541132": "Leleakuhonua",
        "90377": "Sedna",
        "2017 of201": "2017 OF201",
    }
    first = True
    for t in tn:
        if t.a > 1500:
            continue  # keep the frame readable
        ox, oy, oz = orbit(t.a, t.e, t.inc, t.Omega, t.omega, n=300)
        fig.add_trace(
            go.Scatter3d(
                x=ox,
                y=oy,
                z=oz,
                mode="lines",
                line=dict(color="rgba(120,200,255,0.35)", width=1),
                name="detached TNOs" if first else None,
                legendgroup="tno",
                showlegend=first,
                hoverinfo="text",
                text=t.name.strip(),
            )
        )
        first = False
        lo = t.name.strip().lower()
        lab = next((v for k, v in famous.items() if k in lo), None)
        if lab:
            mx, my, mz = _rot(t.a, t.e, t.inc, t.Omega, t.omega, np.array([0.0]))  # perihelion
            fig.add_trace(
                go.Scatter3d(
                    x=mx,
                    y=my,
                    z=mz,
                    mode="markers+text",
                    marker=dict(size=3, color="#7fd3ff"),
                    text=[lab],
                    textposition="bottom center",
                    textfont=dict(color="#7fd3ff", size=8),
                    showlegend=False,
                    hoverinfo="text",
                )
            )

    # Planet Nine: literature best-estimate orbit + most-likely current location
    P9 = dict(a=500.0, e=0.25, i=16.0, Om=100.0, om=150.0)
    ox, oy, oz = orbit(P9["a"], P9["e"], P9["i"], P9["Om"], P9["om"])
    fig.add_trace(
        go.Scatter3d(
            x=ox,
            y=oy,
            z=oz,
            mode="lines",
            line=dict(color="#ff4d4d", width=4, dash="dash"),
            name="Planet Nine orbit (best estimate, uncertain)",
        )
    )
    ax, ay, az = _rot(
        P9["a"], P9["e"], P9["i"], P9["Om"], P9["om"], np.array([math.pi])
    )  # aphelion
    fig.add_trace(
        go.Scatter3d(
            x=ax,
            y=ay,
            z=az,
            mode="markers+text",
            marker=dict(size=8, color="#ff2d2d", symbol="diamond"),
            text=["Planet Nine?  (most likely here, near aphelion ~575 AU)"],
            textposition="top center",
            textfont=dict(color="#ff6b6b", size=11),
            name="Planet Nine (predicted position)",
        )
    )

    rng = 650
    fig.update_layout(
        title=dict(
            text="Solar System + the candidate outer perturber (heliocentric ecliptic, AU)<br>"
            "<sub>Planet Nine location is a best estimate with large uncertainty; "
            "blue = real detached TNOs (the anomaly). Drag to rotate, scroll to zoom.</sub>",
            font=dict(color="#e0e0e0"),
        ),
        paper_bgcolor="#06080f",
        scene=dict(
            xaxis=dict(
                title="x (AU)",
                range=[-rng, rng],
                backgroundcolor="#06080f",
                color="#888",
                gridcolor="#1a2030",
            ),
            yaxis=dict(
                title="y (AU)",
                range=[-rng, rng],
                backgroundcolor="#06080f",
                color="#888",
                gridcolor="#1a2030",
            ),
            zaxis=dict(
                title="z (AU)",
                range=[-rng / 2, rng / 2],
                backgroundcolor="#06080f",
                color="#888",
                gridcolor="#1a2030",
            ),
            aspectmode="manual",
            aspectratio=dict(x=1, y=1, z=0.5),
        ),
        legend=dict(font=dict(color="#ccc")),
        margin=dict(l=0, r=0, t=70, b=0),
    )

    out = ROOT / "solar_system_planet_nine_3d.html"
    fig.write_html(str(out), include_plotlyjs="cdn")
    print(f"wrote {out}")
    print(
        f"  {len([t for t in tn if t.a <= 1500])} real detached TNOs plotted; "
        f"Planet Nine placed at literature best-estimate aphelion (a={P9['a']}, e={P9['e']}, i={P9['i']} deg)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
