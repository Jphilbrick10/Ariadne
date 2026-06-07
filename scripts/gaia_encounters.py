"""Search Gaia DR3 for stars that pass close to the Sun -- the close stellar
encounters that perturb the Oort cloud / outer solar system (the "flyby" alternative
to a present-day Planet Nine, and a real catalog in its own right).

Pulls nearby stars with full 6D phase space (parallax + proper motion + radial
velocity), computes each one's closest linear approach to the Sun, and ranks by the
minimum distance. Sub-parsec encounters are the dynamically interesting ones (Gliese
710 grazes at ~0.06 pc in ~1.3 Myr; Scholz's star passed ~0.25 pc ~70 kyr ago).

  python scripts/gaia_encounters.py

Honest scope: linear (straight-line) approximation -- valid for recent encounters
(<~few Myr); beyond that the galactic tide bends trajectories and a full orbit
integration (Bailer-Jones et al.) is needed. This is a screen for the closest cases.
A flyby that could have detached the Sednoids was ~4 Gyr ago -- that star is long
gone and NOT findable here; this rules out / catalogs RECENT close encounters.
"""

from __future__ import annotations

import sys

import numpy as np

PC_PER_KMS_MYR = 1.022712  # km/s -> pc/Myr


def main():
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.gaia import Gaia

    Gaia.ROW_LIMIT = 60000
    print("=== Gaia DR3 close stellar-encounter search ===", flush=True)
    # nearby stars with full 6D phase space and good astrometry
    q = (
        "SELECT source_id, ra, dec, parallax, pmra, pmdec, radial_velocity "
        "FROM gaiadr3.gaia_source "
        "WHERE parallax > 8 AND radial_velocity IS NOT NULL "
        "AND parallax_over_error > 8 AND ruwe < 1.4"
    )
    print("  querying Gaia (nearby 6D sample)...", flush=True)
    t = Gaia.launch_job_async(q).get_results()
    print(f"  {len(t)} stars with 6D phase space within ~125 pc", flush=True)

    plx = np.array(t["parallax"], float)
    d_pc = 1000.0 / plx
    c = SkyCoord(
        ra=np.array(t["ra"]) * u.deg,
        dec=np.array(t["dec"]) * u.deg,
        distance=d_pc * u.pc,
        pm_ra_cosdec=np.array(t["pmra"]) * u.mas / u.yr,
        pm_dec=np.array(t["pmdec"]) * u.mas / u.yr,
        radial_velocity=np.array(t["radial_velocity"]) * u.km / u.s,
    )
    r = np.c_[
        c.cartesian.x.to(u.pc).value, c.cartesian.y.to(u.pc).value, c.cartesian.z.to(u.pc).value
    ]
    v = (
        np.c_[
            c.velocity.d_x.to(u.km / u.s).value,
            c.velocity.d_y.to(u.km / u.s).value,
            c.velocity.d_z.to(u.km / u.s).value,
        ]
        * PC_PER_KMS_MYR
    )
    rv = np.einsum("ij,ij->i", r, v)
    v2 = np.einsum("ij,ij->i", v, v)
    t_min = -rv / v2  # Myr
    d_min = np.linalg.norm(r - (rv / v2)[:, None] * v, axis=1)  # pc
    order = np.argsort(d_min)

    sub_pc = int(np.sum(d_min < 1.0))
    print(f"  encounters within 1 pc (Oort-cloud-grazing): {sub_pc}", flush=True)
    print(f"  {'d_min(pc)':>9} {'t_min(Myr)':>10} {'now(pc)':>8}  source_id", flush=True)
    for k in order[:20]:
        print(
            f"  {d_min[k]:9.3f} {t_min[k]:+10.2f} {d_pc[k]:8.1f}  {int(t['source_id'][k])}",
            flush=True,
        )
    closest = order[0]
    print(
        f"\n  closest approach: {d_min[closest]:.3f} pc at t={t_min[closest]:+.2f} Myr "
        f"(Gaia DR3 {int(t['source_id'][closest])})",
        flush=True,
    )
    print(
        "  (linear approx; sub-0.5 pc cases warrant full galactic-potential "
        "integration. The Sednoid-detaching flyby, if any, was ~4 Gyr ago and is "
        "unrecoverable -- this catalogs RECENT encounters only.)",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
