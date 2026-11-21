# CLAUDE.md

CZ-8A launch profiles: models and the 3D page. Git repo `exoplanet5/CZ-8A` (public). GitHub Pages serves `docs/` from `main` at https://exoplanet5.github.io/CZ-8A/.

**Layout.**
- Scripts and the shared model sit in the root.
- Each flight's NOTAM file, inputs (`outside-ref/`), `output/` and README live in `missions/<flight>/`, which is `cz8a.MISSION` and `cz8a.OUT`. `j2_plan.py` names the folder itself.
- The root `README.md` is the repository overview, and each mission README is that flight's record of results, assumptions and sources.
- `data/` holds the Space-Track caches shared by the flights.

**Missions.** One so far: **CZ-8A Y13** (SatNet / Guowang LEO group 30, i = 50°) from Hainan commercial LC1, launched 2026-11-20 19:59 BJT, in `missions/Y13_20261120/`.
- Lower stages fitted to the RPHI NOTAM zones and to four ground sightings (`sightings.py`).
- S2 long coast of 40.6 min onto the inserted orbit (`outside-ref/outside.tle`) with zero phase shift, then the parking orbit.
- Checked against the catalog phase of earlier flights.

**Revision 5 is the settled Y13 profile** (user, 2026-11-21: "settle down this version of Y13 profile is correct, we will only tweak it by more sighting data"). Do not refit, retune or re-argue the ascent, coast or vehicle, including the throttle open point. Change it only when new sightings come in: add them to `sightings.SIGHTINGS` and rerun the pipeline. **Y12 (2026-243) flew exactly the Y13 profile** (user: "Y12 is exactly following Y13, just remember no need did more"). Do no further Y12 modelling.

**More CZ-8A missions will be added.** Give each its own `missions/<flight>_<yyyymmdd>/`. Turn `cz8a.MISSION` and the Y13 constants in `build_traj.py`, `sightings.py` and `j2_plan.py` into a mission key, as the CZ-2D repository does with `build_traj.MISSIONS`. Add the flight to the page with `new_page.py --add`.

## Pipeline (order matters; run from the root; `output/` = `missions/Y13_20261120/output/`)

```
~/.venvs/astro313/bin/python fetch_data.py     # Space-Track, 3 queries, cached by the space-track skill
~/.venvs/astro313/bin/python j2_plan.py        # -> output/j2_plan.json (K_y13_deg is read by build_traj.py)
~/.venvs/astro313/bin/python build_traj.py     # --t0 ISO (default 11:59:59), --payload kg, --profile attached|direct|longcoast, --shift s, --u-apo deg, --full-thrust, --sens
~/.venvs/astro313/bin/python validate.py       # zone-predicted K on Y7/Y10/Y11/Y12 (full thrust, zone-only fit)
~/.venvs/astro313/bin/python validate.py phase # catalog phase lead at T0 per flight and profile (~3 min)
~/.venvs/astro313/bin/python compare_outside.py # model vs outside-ref/outside.tle, three profiles at build_traj.T0_NOMINAL
~/.venvs/astro313/bin/python coast_tradeoff.py  # phase shift and S2 left vs coast, T0 second, S2 load, full thrust (~1 min)
~/.venvs/astro313/bin/python make_map.py       # runs visibility.compute on the map extent, writes both PNGs
~/.venvs/astro313/bin/python web_export.py     # page data -> docs/data/cz8a-y13.json, then new_page.py --add to validate
```

## 3D page

- The page was made with the `launch-page` skill from `~/launch/launch-page-template` on 2026-11-21. Never copy `docs/` from an older page.
- Profile `satnet-50` ("SatNet · 50°", family `low-inclination`) holds mission `cz8a-y13` ("CZ-8A Y13 · 20 Nov 2026"), which is the default.
- Re-register after an export:
  ```
  new_page.py . --add --mission docs/data/cz8a-y13.json --id cz8a-y13 --label "CZ-8A Y13 · 20 Nov 2026" --profile satnet-50 --profile-label "SatNet · 50°" --family low-inclination
  ```
- Flight segments: `ascent` (to S2 cutoff 1), `coast`, `burn2` and `orbit` (2 revolutions after burn 2; `insertion_met` = cutoff 2). The fall tracks are `boosterfall`, `corefall` and `fairfall`.
- The CSVs round times to 1 ms, so `web_export.uniform` snaps boundary rows to the exact event times. The validator requires adjacent segments to share their endpoints.
- No `meta.jelly`: the powered ascent is in shadow, and the only sunlit burn is the apogee burn.
- `sky_mark` is the stage-1 separation, the timed sighting. `observer_home` is the sighting site, an extra city named "Lingshui sighting".
- Publish with the skill's `publish_page.py --repo exoplanet5/CZ-8A`. Do not publish `.shots/`.

`cz8a.py` is the model. It imports generics from `../launch_traj_NROL-97/nrol97.py` and the guided-burn kernel from `../launch_traj_CERES-2/ceres2.py`. `j2_plan.py` imports `../launch_traj_12/j2_tuning.py`. Never copy them.

`cz8a` puts those two folders at the end of `sys.path`, and moves NROL-97 back to the end after `ceres2` re-inserts it at the front. They hold their own `build_traj.py` and `make_map.py`, which would otherwise shadow this project's.

## Fit choices

- **Pad 19.5970°N 110.9307°E** (user, 2026-11-21: "initial ascent is from the offshore, should be the pad").
  - GCAT's `lp.tsv` row for HCSLS LC1 says 19.5685°N, which is 3.1 km south and 2.4 km out at sea.
  - The value used is GCAT's own site record (19.5968, 110.9308), which agrees with the surveyed pad in `../launch_traj_CZ-10B`.
- **T0 minute 11:59 UTC** = zone-1 B) time + 8 min, the lead on every CZ-8A flight except Y12. Xinhua confirmed 19:59 BJT.
  - The seconds are free since revision 5: zero phase shift holds from about 11:59:30 to 11:59:59, and 11:59:00 needs +26 s.
  - The nominal keeps 11:59:59, where the S2 margin is largest.
  - The J2 plan (Y12's K and D toward 2025-295) gives 11:58:18. The user suggested the 2025-295 target.
- **Y13 copies Y12 (2026-243), not Y10/Y11.** It parks at about 1049–1051 km mean with K = 309.70°. Y12/Y13 zones lie on the 140° corridor, Y10/Y11 on 142.8°, and the corridor predicts K (`validate.py`).
- **Zone 1 (B6120) = boosters; zone 2 (B6122) = core + fairing.** The user called the NOTAMs "first stage and fairing"; the model puts the boosters in zone 1.
- Lower stage: published sequence. Booster separation is 171.9 s (CZ-8 Y1); stage separation 198.5 s, S2 ignition 200.8 s and fairing 209.4 s come from the CZ-8A sequence. Propellant is consumption at the burn times.
- **S2 load 30.2 t** (user, 2026-11-21). Earlier loads: 26.0 t in revisions 1–2, 27.2 t in revisions 3–4. `cz8a.VEH` is the published full-thrust vehicle: lift-off 371.0 t.
- **Revision 5 (2026-11-21): lower stages fitted to the user's ground sightings.**
  - The observer is at N 18°26′14″ E 110°04′28″. The sightings are az/alt 45.5/13.5, 72.4/25.5, 89.4/28.6 (untimed) and the stage-1 separation at 94.6/28.9, timed at `t_csep`.
  - Steering law (`cz8a.fly_lower` with `t_q`):
    - a zero-AoA turn to t_q;
    - a linear-tangent pitch from the flight path at rate b_boost to booster separation;
    - the same law continued at rate b_core through the core phase, with no attitude step.
  - `fit_lower(..., sight=True)` fits kick, b_core, az, t_q and b_boost to 4 zone residuals (5 km) and 5 sighting residuals (0.3°). Start from `build_traj.X0_LOWER`: the fit has a second, bad minimum.
  - **Nominal vehicle `build_traj.VEH_NOM`: core engines at 77.5 % from T+10 s to booster separation.**
    - Lift-off is 341.3 t, against the published 371 t. Flag this; it is the open point of revision 5.
    - Only with that much less energy do the sightings (≤ 0.77°) and the zones (≤ 4.1 km) fit together. Full thrust fits the sightings to 2.3° at best; a constant altitude bias makes it worse.
    - The throttle stands in for whatever lowers the lower-stage energy; it is not claimed as fact.
    - Throttling from lift-off instead would contradict the published 480 t lift-off thrust.
  - Result: boosters off at 88 km (24°), core cutoff at 118 km, S2 ignition at 124 km, fairing at 135 km. Falls: boosters 12:11:27 and core 12:13:08 (both on their centroids), fairing 12:21:38, all inside their zones and windows.
  - `AIM_KM` is (0, 0). The 30 km booster aim of revision 4 existed only to lift the fairing over 100 km.
- **Profile: `attached` is nominal** (user's constraints from revision 4):
  - long coast 40–48 min after SECO-1;
  - the inserted orbit's plane and altitude met, its phase allowed to shift by tens of seconds;
  - zone fit kept;
  - T0 inside the minute.
  - Revision 5 meets all of them with zero shift: coast 40.6 min (apogee at u = 323°), burn 1 627 s (pitch 71° → −36°) to a 679 × 1061 km transfer, S2 1.01 t left.
  - The whole 40–48 min range closes at zero shift (0.18 t left at 48.6 min).
  - `build_traj.attach` returns shift 0 when the S2 can pay for it, else the smallest shift leaving `S2_RESERVE`.
- **Direct profile: no solution with 30.2 t.** The published 508.6 s burn 1 gives 3.80 km/s instead of 4.32 km/s with 27.2 t. `compare_outside.py` and `validate.py phase` record the failure instead of crashing.
  - Revision 2's inference that Y12 flew direct rested on the flat ascent. From the sighting-fitted ascent the long coast reaches Y12's orbit moved to Y13 with zero shift; the Y13 README states this as an inference.
- Earlier revisions, kept for the record in the Y13 README:
  - Revision 4 (flat ascent, 83 s ahead, 27.2 t).
  - Revision 3 (steep branch, booster separation at 136 km, rejected by the user).
  - Revision 2 (direct).
  - Revision 1 (Y10-type long coast, 405 s ahead).
  - `longcoast_check.py` is a revision 2 check (old pad, 26.0 t).
- The disposal is a propellant dump 300 s after deployment (attached). Deployment is at +120 s. For the direct and long-coast profiles the apex radius is iterated to Y12's first-TLE mean altitude (1049.15 km).
- **Maps are 1:1** (user, 2026-11-20): `make_map.equal_figure` sizes the figure from the extent (3960 px long side, 300 dpi), and `lock_aspect` sets equal aspect.
- **Visibility:** the powered ascent is entirely in shadow, so the ≥10° zone is the self-luminous flame.
  - It covers 105.1–144.9°E and 18.9°S–24.3°N, which includes Guangzhou and Hong Kong.
  - The regional extent is [103, 148, −22, 27] (`EXTENT` in both `make_map.py` and `visibility.py`), 3288 × 3960 px.
  - The launch-corridor inset (`make_map.INSET`, placed on the open ocean east of the Philippines because the upper left holds the pad) shows the observer, the sightlines and the staging.
  - The Sun −12° line crosses the south-west corner.
  - First sunlight is at T+30.0 min over NZ. The apogee burn is in daylight over the south-east Pacific.

## Numerical note

`cz8a.solve_attached` corrects its along-track bias at burn-2 ignition (fixed 2026-11-21: it used to leave a constant offset of half the burn-2 length, 12–20 km).

The ceres2 kernel switches thrust off for tau > 1. Burns ending on a RK4 stage could round past 1 and make burn-time derivatives noisy, which stalled least_squares. `cz8a.EPS_T` pads every thrust window by 1 µs; keep it.

## After launch

1. When the catalog has Y13, fetch its first TLEs (space-track skill).
2. Compare i, Ω, mean altitude and `cz8a.lead_at_t0` with the Y13 README. The model sits on the outside set (zero shift).
3. Solve T0 seconds from Ω (K 309.70).
4. Rerun `build_traj.py --t0 <actual>`, `compare_outside.py` and `make_map.py`.
5. Timed sightings (`sightings.SIGHTINGS`) would pin the climb rate and so test the core throttle.
