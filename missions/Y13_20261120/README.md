# CZ-8A Y13 — ascent, long coast, apogee burn and orbit

**3D page: https://exoplanet5.github.io/CZ-8A/?profile=satnet-50&mission=cz8a-y13** (profile SatNet · 50°). The repository overview is in the [root README](../../README.md).

**Status: revision 5 (2026-11-21), after launch; still a model, not telemetry.** The flight is CZ-8A Y13 (长征八号甲遥十三, 卫星互联网低轨30组). It lifted off from Hainan commercial pad LC1 at **2026-11-20 19:59 Beijing time** (Xinhua). It carries nine Guowang satellites into a 1050 km parking orbit, from which J2 drift brings them into the plane of **2025-295** (CZ-12 Y4, WHDW 16) as its second group, on the Y12 (2026-243) plan.

**Revisions.**

1. **Revision 1** flew Y10's long coast from a normal ascent. It ran 23.7° ahead of the inserted orbit (the outside element set that caught the payload near T+100 min).
2. **Revision 2** switched to the published direct-ascent sequence (cutoff T+709 s, relight T+1761.5 s), which matched Y12's catalog phase to 2°.
3. **Revision 3** met the inserted orbit through a steep lower-stage branch (booster separation at 136 km). Withdrawn: the staging altitude was not credible.
4. **Revision 4** kept a normal ascent and a 40 min coast and arrived 83 s ahead of the inserted orbit, with a 27.2 t S2 burned to 0.1 t.
5. **Revision 5, nominal (the user's data, 2026-11-21):**
   - the **pad moved to the real LC1** (19.5970°N 110.9307°E). GCAT's LC1 row put it 3.1 km south, 2.4 km out at sea;
   - **S2 load 30.2 t**;
   - the lower stages **fitted to four ground sightings** from Lingshui as well as both NOTAM zones. That needs a much slower, higher climb: boosters off at 88 km, S2 ignition at 124 km;
   - with that ascent, the **40.6 min long coast lands exactly on the inserted orbit (zero phase shift)** at T0 11:59:59, with 1.0 t of S2 propellant left;
   - the price: the sightings and the zones together fit only if the core runs at 77.5 % thrust from T+10 s to booster separation. Lift-off is then **341.3 t, not the published 371 t** ([Ground sightings](#ground-sightings)).

**Revision 5 is the settled Y13 profile** (2026-11-21): the ascent corrected with the sightings, the long coast, and Y12's inserted elements moved in plane to Y13. It will change only with more sighting data. Y12 (2026-243) flew the same profile, so it is not modelled separately.

## Result

**T0 = 2026-11-20 11:59:59 UTC**, inside the published minute (see [T0](#t0)).

| Event | MET | UTC | Altitude | Inertial speed | Where |
|---|---:|---|---:|---:|---|
| Booster separation | 171.9 s | 12:02:50.9 | 87.9 km | 3.16 km/s (2.87 Earth-relative), flight-path 24° | 111 km downrange |
| Core cutoff / separation | 196.5 / 198.5 s | 12:03:15.5 | 118.5 km | 4.07 km/s, 21° | 185 km |
| S2 ignition 1 | 200.8 s | 12:03:19.8 | 124.1 km | 4.05 km/s, 20° | 200 km |
| Fairing jettison | 209.4 s | 12:03:28.4 | 135.2 km | 4.06 km/s | 230 km |
| **S2 cutoff 1** | **828.1 s** | **12:13:47.1** | **717.4 km** | 7.57 km/s, flight-path +0.95° | 1.15°S 127.18°E (northern Moluccas), 2909 km |
| **Long coast** | **2438 s (40.6 min)** | | 717 → 1061 km | | over northern Australia, NZ, the South Pacific |
| **S2 ignition 2, apogee burn** | **3266.1 s** | **12:54:25.1** | **1054.5 km** | Δv 95 m/s in 5.3 s | 27.7°S 90.1°W (south-east Pacific), on the inserted orbit at u = 323° |
| Deployment (assumed +120 s) | 3391.3 s | 12:56:30 | 1053.8 km | | |
| Stage disposal (propellant dump, assumed +300 s) | 3691.3 s | 13:01:30 | 1052.8 km | perigee → 904.6 km | |

| Orbit | Value |
|---|---|
| Transfer (osculating at cutoff 1) | 678.7 × 1061.2 km, i 50.023°, Ω 189.191° |
| Parking (osculating after burn 2) | 1045.6 × 1057.7 km, i 50.008°, Ω 189.056° |
| Parking (SGP4 mean, fitted TLE) | **mean altitude 1049.48 km**, i 50.0032°, e 0.00123, ω 251.2° |
| Inserted orbit (outside set) | mean altitude 1049.44 km, i 50.0058°, e 0.00123, ω 252.3° |
| **Offset from the inserted orbit, T+70 to T+240 min** | **along-track within ±0.4 km (0.0 s), ΔΩ −0.004°, Δi −0.003°, Δr < 0.1 km** |
| Phase lead at T0 | −15.19° (inserted orbit −15.20°) |
| Stage after disposal | 904.6 × 1052.2 km |

**S2 budget:** burn 1 uses 28.91 t in 627.3 s, burn 2 0.24 t and the disposal 0.04 t. That leaves **1.01 t of the 30.2 t load**. Lift-off is **341.3 t** with the throttled core, against 371.0 t for the same vehicle at full thrust.

The TLE pair is [CZ-8A_Y13_prelaunch.tle](output/CZ-8A_Y13_prelaunch.tle). Catalog numbers 99913 and 99914 and the designator 26999 are placeholders.

## Ground sightings

The user supplied four sightings from **N 18°26′14″, E 110°04′28″** (near Lingshui, south-east Hainan), 157 km south-west of the pad: three positions on the climb and the stage-1 (core / S2) separation. The model times the separation at T+198.5 s; the climb sightings carry no time. [sightings.py](../../sightings.py) scores an untimed sighting by its angular distance from the model's sky path and the timed one by its azimuth and altitude error. The zone fit ([cz8a.fit_lower](../../cz8a.py), `sight=True`) weighs 0.3° of sighting against 5 km of zone miss.

| Sighting | Seen az / alt | Model az / alt | Off | Model time | Altitude | Range from observer |
|---|---|---|---:|---:|---:|---:|
| Climb 1 | 45.5° / 13.5° | 45.26° / 13.79° | 0.37° | T+122 s | 39.4 km | 158 km |
| Climb 2 | 72.4° / 25.5° | 72.19° / 26.24° | 0.77° | T+169 s | 84.4 km | 186 km |
| Climb 3 | 89.4° / 28.6° | 89.41° / 28.25° | 0.35° | T+191 s | 111.5 km | 229 km |
| Stage-1 separation | 94.6° / 28.9° | 94.47° / 28.19° | 0.72° | T+198.5 s | 121.1 km | 249 km |

For the untimed sightings, "model az / alt" is the closest point of the model's sky path. Revision 4's flat ascent misses the same four sightings by 2.6°, 6.6°, 8.8° and 13.2°: it was too low (93 km at separation) and too far downrange.

**What the sightings need.** Intersected with the zone corridor, the sightlines put the vehicle at about 41 km altitude 30 km downrange, 87 km at 108 km and 133 km at 202 km. The model now flies:

- a gravity turn after a 0.395° kick at T+10 s;
- from T+136 s, a linear-tangent pitch (rate −1.06 tan per 100 s), continued without an attitude step through the core phase (+0.51);
- azimuth 139.24°.

Boosters separate at 88 km and 24° flight-path angle, the core at 118 km and 21°. Both zones still hold to 4 km.

**Why the core is throttled.** At these flight-path angles the full-thrust, 371 t vehicle carries too much speed for the zones. The boosters would land 270 km beyond zone 1 and the core 530 km beyond zone 2. Holding the zones, the best full-thrust fit misses the sightings by up to 2.3°. The zone-plus-sightings fit wants about 0.3 km/s less at booster separation (2.87 km/s Earth-relative against 3.0–3.2). The model gets it by running the core at 77.5 %, the CZ-8 throttle level, from T+10 s to booster separation. Since the stage loads are consumption at the burn times, the core then carries 30 t less propellant:

| Core at 77.5 % | Lift-off | Largest sighting miss | Largest zone miss |
|---|---:|---:|---:|
| from T+0 s (lift-off thrust 4260 kN, against the published 4800 kN) | 339.5 t | 0.85° | 3.4 km |
| **from T+10 s (nominal)** | **341.3 t** | **0.77°** | **4.1 km** |
| from T+30 s | 345.0 t | 1.37° | 2.4 km |
| from T+60 s | 350.5 t | 1.87° | 6.2 km |
| from T+120 s | 361.5 t | 2.42° | 9.5 km |
| never (full thrust) | 371.0 t | 2.30° | 5.8 km |

The T+30 to T+120 s rows need the booster-phase pitch rate at its bound. A constant altitude error in the sightings does not rescue full thrust: lowering all four by 1–3° makes the fit worse.

**Open point.** Within this 3-DOF model, the published lift-off mass (371 t), full YF-100 thrust to the published separation times (171.9 s, 198.5 s), the zones and the sightings cannot all hold. The throttle is a stand-in for whatever lowers the lower-stage energy: a partial load, heavier stages, or a booster separation earlier than CZ-8 Y1's 171.9 s. `build_traj.py --full-thrust` flies the 371 t vehicle on the same law for comparison.

## Phase, coast and T0

The slower ascent is what makes the long coast arrive on time. It reaches S2 ignition 0.6 km/s slower and 30 km higher than revision 4. Burn 1 (pitch 71° → −36°, 627 s) cuts off at 717 km on a 679 × 1061 km transfer, and the 40.6 min coast ends on the inserted orbit. **The phase shift is zero, and the S2 keeps 1.0 t.**

The table shows what the margin buys ([coast_tradeoff.py](../../coast_tradeoff.py), [coast_tradeoff.json](output/coast_tradeoff.json)). A shift is solved only when the S2 cannot pay for zero; "S2 left" is after the disposal dump.

| Case | Phase shift | S2 left | Coast | Cutoff 1 |
|---|---:|---:|---:|---|
| **Nominal: T0 11:59:59, apogee u = 323°, 30.2 t** | **0 s** | **1.01 t** | 40.6 min | T+828 s, 717 km |
| Apogee u = 330° | 0 s | 0.80 t | 42.6 min | T+833 s, 732 km |
| Apogee u = 337° | 0 s | 0.58 t | 44.6 min | T+838 s, 746 km |
| Apogee u = 344° | 0 s | 0.38 t | 46.6 min | T+843 s, 759 km |
| Apogee u = 351° | 0 s | 0.18 t | 48.6 min | T+848 s, 773 km |
| T0 11:59:45 | 0 s | 0.61 t | 40.7 min | T+838 s, 749 km |
| T0 11:59:30 | 0 s | 0.18 t | 40.8 min | T+848 s, 782 km |
| T0 11:59:00 | +26 s | 0.09 t | 40.8 min | T+850 s, 788 km |
| S2 load 28.7 t (sightings then ≤ 1.2°) | 0 s | 0.36 t | 40.9 min | T+810 s, 740 km |
| S2 load 27.2 t (sightings then ≤ 1.5°) | +17 s | 0.09 t | 41.1 min | T+784 s, 724 km |
| Full thrust, 371 t (sightings ≤ 2.3°) | 0 s | 0.83 t | 40.6 min | T+832 s, 720 km |

So:

- the whole 40–48 min coast range closes with zero shift;
- **T0 seconds are no longer pinned by the phase**: zero shift holds from about 11:59:30 to 11:59:59, and 11:59:00 needs only +26 s. The nominal stays at 11:59:59, the largest margin;
- the plane cannot pin them either. The zone corridor reaches S2 ignition at i = 47.7°, so the S2 yaws 2.4° into the 50° plane whatever the second; across the minute that changes by 0.03°.

**The direct profile no longer closes.** With 30.2 t loaded, the published 508.6 s first burn gives 3.80 km/s instead of 4.32 km/s at 27.2 t, and from the slower staging no steering reaches a 1050 km apex in time for the T+1761.5 s relight. The best try peaks at 966 km on a sub-orbital arc. `--profile direct` and `validate.py phase` now report this instead of a solution.

## Phase check against the outside orbit and the catalog

**The outside element set is Y12's orbit moved to Y13.** [outside-ref/outside.tle](outside-ref/outside.tle) (catalog number A1189 = 101189) has the first catalog TLE of 2026-243E (Y12 satellite 05, epoch 2026-10-25 11:59:22.887) with i, e, ω, n and B* unchanged. Only three things were changed:

- the epoch, shifted +26 d 09:57:00;
- Ω, shifted by exactly the GMST change over that interval;
- M, advanced by 1.5523° (27.5 s).

It reported a catch near T+100 min, so it stands for "Y13 flew like Y12".

**Comparison with the outside orbit** ([compare_outside.py](../../compare_outside.py), [outside_comparison.json](output/outside_comparison.json)). Offsets are model minus outside at T+100 min. A pure T0 error dt would move the node by +0.2507°/min × dt and the phase by −n·dt, so both must give the same dt.

| Model | Along-track | Δu | ΔΩ | Δi | T0 shift implied by phase / by node |
|---|---:|---:|---:|---:|---:|
| **Rev. 5, attached long coast, T0 11:59:59 (nominal)** | +0.05 km (0.0 s) | +0.003° | −0.004° | −0.003° | 0 s / +1 s |
| Rev. 5 lower stages, Y10-type long coast (200 km perigee) | +2126 km (+290 s) | +16.65° | −0.017° | −0.008° | +295 s / +4 s |
| Rev. 5 lower stages, direct ascent | no solution with the 30.2 t S2 | | | | |
| Rev. 4, attached long coast, T0 11:59:59 | +609 km (+83 s) | +4.71° | −0.010° | −0.001° | +83 s / +3 s |
| Rev. 4 lower stages, direct ascent, T0 11:59:59 | +60 km (+8 s) | +0.43° | +0.053° | −0.008° | +8 s / −13 s |
| Rev. 1, long coast, T0 11:59:00 | +2967 km (+405 s) | +23.73° | −0.271° | −0.010° | +420 s / +65 s |
| Rev. 2, direct ascent, T0 11:59:00 | +480 km (+65.5 s) | +3.83° | −0.198° | −0.007° | +68 s / +47 s: consistent |
| Rev. 3, steep ascent and 49.7 min coast, T0 11:59:00 (withdrawn) | +12 km (+1.6 s) | +0.10° | −0.005° | −0.001° | +1.7 s / +1.2 s |

Revision 4 and earlier rows used the old pad (19.5685°N) and a 26.0–27.2 t S2. The offsets are steady from T+70 to T+240 min, so the mean motions agree. The minute between revision 2's T0 and its best fit is partly the outside set's own construction. Its node assumes Y13 − Y12 = 26 d 09:57:00, i.e. 12:00:00 against Y12's 02:03:00, while its phase assumes 09:56:32.5. The rest is the model's 2° phase residual on Y12 (below). The published minute allows 11:59:00–11:59:59, so T0 seconds cannot be separated from that residual.

**Catalog phase of each flight** ([validate.py phase](../../validate.py), [phase_validation.json](output/phase_validation.json)). The quantity is the phase lead at T0: the argument of latitude of the parking orbit, propagated back to the published T0, minus the pad's, in the same plane. Each run uses that flight's own zones, K, published T0 minute and first-TLE mean altitude. These flights have no sightings, so they fly the zone-only flat ascent at full thrust:

| Flight | Parking | Catalog lead | Long-coast model, rev. 5 (old) | Direct-ascent model, rev. 5 (old: 27.2 t, old pad) |
|---|---:|---:|---:|---:|
| Y10 (2026-176) | 1100.4 km | +4.57° | **+5.14° (+0.57°, 10 s)** (+6.01°) | no solution (−14.38°) |
| Y12 (2026-243) | 1049.1 km | −16.81° | +3.97° (+20.8°) (+4.89°) | no solution (**−14.69°, +2.1°**) |
| Y7 (2026-007, older 920 km mission) | 919.4 km | −23.72° | | no solution (−15.19°) |
| Y13 (outside set), T0 11:59:59 | 1049.4 km | −15.20° | +1.45° (sighting-fitted ascent, perigee cutoff); **−15.19° attached, 40.6 min coast** | no solution |

Y10 flew the hour-long coast, as Beijing News reported, relighting near T+60 min, and the long coast still matches its catalog phase to 0.6°. **Revision 2's conclusion that Y12 flew the direct sequence rested on the flat ascent.** Y13's outside set is Y12's first catalog TLE moved to Y13, and from the sighting-fitted ascent the 40.6 min long coast reaches it with zero shift. So if Y12 climbed like Y13, its catalog phase fits the long coast too. That is an inference: Y12 has no sightings. Y7 was 8.6° off in revision 2: the older 920 km mission evidently used other timing.

**Revision 2's evidence for the direct profile** (26.0 t load, flat ascent, old pad):

- The S2 now burns 23.4 t + 1.2 t and keeps **1.27 t**; the long coast kept 5.1 t.
- The fitted parking ω (252.6°) matches the outside/Y12 ω (252.3°).
- Disposal at Y12's stage apogee gives Y12's stage orbit, about 905 × 1051–1061 km with ω about 152°.
- The relight time is the published one; a later relight (1950 s) gives −13.4°, further from the catalog.

### Can a 45–50 min coast reach the inserted orbit from the flat ascent? (revisions 2–4)

**Not exactly, from the flat ascent with T0 inside 19:59 BJT**; revision 4 got to within 83 s. The sighting-fitted ascent of revision 5 removes the problem (above). The check below was run for revision 2 (T0 11:59:00, 26.0 t load, old pad, perigee cutoffs). [longcoast_check.py](../../longcoast_check.py) ([longcoast_check.json](output/longcoast_check.json)) flies S2 burn 1 from the zone-fitted ascent into a transfer orbit with apogee at the parking radius. Cutoff varies from 150 to 600 km and the flight-path angle from −10° to +20°. Each case then coasts to apogee and circularizes. The inserted orbit's lead at T0 is −18.39°.

| Cutoff | Coast | Lead | Miss | T0 shift that would absorb it | S2 left |
|---|---:|---:|---:|---:|---:|
| 150 km perigee | 48.1 min | +5.95° | +24.3° | +7.2 min | 5.3 t |
| 200 km perigee | 48.4 min | +4.78° | +23.2° | +6.8 min | 5.2 t |
| 300 km perigee | 48.9 min | +2.27° | +20.7° | +6.1 min | 4.5 t |
| 450 km perigee | 49.8 min | −1.82° | +16.6° | +4.9 min | 2.4 t |
| 600 km perigee | 50.6 min | −6.08° | +12.3° | +3.6 min | 0.1 t |
| 300 km, flight-path +10° (sub-orbital) | 17.6 min | −12.13° | +6.3° | +1.9 min | 3.0 t |

The reason is geometric. The lead splits into a cutoff term (u_cut − n·t_cut − u_pad), which the ascent pins at −10° to −14°, and a coast term (Δu − n·Δt). A coast of 45 min or more from a low perigee to a 1050 km apogee covers about 180° faster than the circular orbit moves, so its coast term is +7° to +16° and never negative. Matching −18.4° would need a coast slower than the 1050 km circular orbit, which means spending the coast above 1050 km. Even using all the S2 propellant (600 km perigee) leaves 12°. From that ascent, reconciling a long coast with the inserted orbit would need T0 3.6–7 min later (12:02:40–12:06), against Xinhua's 19:59. Revision 3 slowed the lower stages instead (withdrawn: staging at 136 km). Revision 4 slowed the S2 burn (cutoff at 643 km after 581 s) as far as the propellant allowed. Revision 5's lower stages, fitted to the sightings, are slow enough on their own.

The official "hour-level" coast was announced for Y10, and Y10's catalog phase does match the long-coast model, to 1.8°. In the direct profile the S2 coasts from cutoff 1 (T+11.8 min) to its last burn at Y12's stage apsis (T+56.6 min). That is 44.8 min, interrupted only by the 27 s apex burn at T+29.4 min.

## Stage falls and NOTAM attribution

The NOTAMs give no stage names; the attributions below are the model's.

| Object | Lands | UTC | Range | Zone and window |
|---|---|---|---:|---|
| 2 × K2-1 boosters (β 300 kg/m²) | T+689 s, 12.69°N 116.73°E | 12:11:27 | 986 km | **B6120/26** (11:51–12:17), centroid 979 km |
| K3 core stage (β 400) | T+789 s, 7.58°N 120.71°E | 12:13:08 | 1702 km | **B6122/26** (11:52–12:28), centroid 1699 km |
| Fairing halves (β 30) | T+1299 s, 7.55°N 120.74°E | 12:21:38 | 1707 km | **B6122/26**, 6 min before the window closes |

All three land inside their zones and windows; both stages are aimed at the zone centroids. The fairing jettisons at 135 km. Revision 4 had aimed the boosters 30 km long to lift it to 102 km, and the sightings make that unnecessary.

**Zone 1 is the booster zone, not the core stage.** CZ-8A flies CZ-7's K2-1 boosters, which separate at T+172 s at about 3 km/s and reach about 980 km. CZ-7A's booster zone sits at 884 km (`~/launch/launch_traj_7A_GTO`). The core stage and the fairing both fall in zone 2. The fairing jettisons 11 s after the core stage and comes down 8–9 min after it, slowed by its low ballistic coefficient. That fits zone 2 staying open 11 min longer than zone 1 (36 min against 25–26 min on every CZ-8A flight).

## T0

**NOTAM lead.** On Y1–Y11, the zone-1 B) time opened exactly 8 min before the published lift-off. Y12 is the exception: its windows were 50–60 min and it launched 21 min after B). Y13's windows are the usual 25/36 min, with B) at 11:51, so T0 = 11:59.

**T0 seconds.** Xinhua gives the minute only. In revision 4 the phase pinned them to the end of the minute. In revision 5 the long coast reaches the inserted orbit with zero shift from about 11:59:30 to 11:59:59 ([Phase, coast and T0](#phase-coast-and-t0)), and the S2 yaws into the plane whatever the second. The nominal keeps **11:59:59 UTC**, where the S2 margin is largest (1.0 t); any second of the minute fits within 26 s of phase.

**J2 plan.** The script [j2_plan.py](../../j2_plan.py) imports `../launch_traj_12/j2_tuning.py` and measures, for each recent group, the site constant K = Ω_inj(T0) − GMST(T0) and the offset D = Ω_inj − Ω_target at T0:

| Flight | Group | Published T0 (UTC) | Park | K | Target | D |
|---|---|---|---:|---:|---|---:|
| Y7 | 2026-007 | 01-13 15:25 | 920 km | 310.11° | — | — |
| Y10 | 2026-176 | 08-04 08:52 | 1100.8 km | 308.47° | 2025-187 | +1.31° |
| Y11 | 2026-221 | 09-23 13:32 | 1100.7 km | 308.70° | 2025-285 | +0.95° |
| Y12 | 2026-243 | 10-25 02:03 | 1050.7 km | 309.70° | 2025-311 | +3.98° |

Joining 2025-295 with the **Y12 plan gives T0 = 11:58:18 UTC**, 42 s from the NOTAM T0. The Y10/Y11 plan would put T0 at 11:51:25, the moment zone 1 opens, which no CZ-8A flight has done. Y13 copies Y12 on the evidence: the NOTAM zones, the parking altitude and the plane.

**Why K differs between Y10/Y11 and Y12.** The zones fix the lower-stage corridor, and the corridor fixes K. The script [validate.py](../../validate.py) fits the same full-thrust vehicle (30.2 t S2, corrected pad, zone-only flat ascent) to each flight's own two zones. It then takes the inertial plane at S2 ignition. This prediction uses no TLE:

| Flight | Zone bearing | Fitted azimuth | Kick | K from zones | K measured | Difference |
|---|---:|---:|---:|---:|---:|---:|
| Y7 | 140.7° | 139.79° | 0.406° | 309.62° | 310.11° | −0.49° |
| Y10 | 142.9° | 142.16° | 0.411° | 308.32° | 308.47° | −0.15° |
| Y11 | 143.0° | 142.27° | 0.413° | 308.27° | 308.70° | −0.43° |
| Y12 | 140.1° | 139.41° | 0.408° | 309.84° | 309.70° | +0.14° |
| **Y13** | **140.2°** | **139.43°** | **0.398°** | **309.82°** | (Y12 adopted) | |

Y10 and Y11 flew a corridor 2.8° more southerly, with zone 2 at 1530 km instead of 1700 km. Y12 and Y13 fly the old Y1–Y8 corridor. One pitch kick (0.40–0.41° with the 30.2 t S2; 0.51–0.53° with 27.2 t) reproduces all five booster zones. The sighting-fitted Y13 ascent climbs more steeply, and its plane at S2 ignition gives K = 310.32°.

## Model

[cz8a.py](../../cz8a.py) holds the 3-DOF point-mass model. Frames, J2–J4 gravity, the atmosphere, coasts and falls, elements and the SGP4 fit are imported from `../launch_traj_NROL-97/nrol97.py`. The guided-burn kernel, pad geometry and NOTAM parser come from `../launch_traj_CERES-2/ceres2.py`.

| Item | Value | Source |
|---|---|---|
| Pad | 19.5970°N 110.9307°E | GCAT HCSLS site record (19.5968, 110.9308) and the surveyed pad of `~/launch/launch_traj_CZ-10B`; GCAT's LC1 row (19.5685°N) lies 3.1 km south, at sea |
| YF-100 | 1200 / 1340 kN, Isp 300 / 335 s, 407.9 kg/s | Wikipedia, GCAT |
| K2-1 booster ×2 | dry 6.0 t (GCAT); 70.1 t burned to 171.9 s; 1.0 t residual (assumed) | CZ-8 Y1 booster separation 171.88 s |
| K3 core | dry 12.0 t (GCAT); 1.5 t residual (assumed). Full thrust: 160.3 t burned to 196.5 s. **Nominal: both engines at 77.5 % from T+10 s to booster separation, 130.6 t burned** | CZ-8A separation 198.5 s; CZ-8 throttle level; the sightings (above) |
| S2 (2 × YF-75H/DB) | 200 kN, Isp 442.6 s; dry 3.6 t incl. multi-function cabin (assumed); **30.2 t load** (user, revision 5; 26.0 t in revisions 1–2, 27.2 t in 3–4); three ignitions since Y9 | Wikipedia, GCAT |
| Fairing (5.2 m) | 2.9 t (assumed) | |
| Payload | 9 × 700 kg = 6.3 t | GCAT estimate |
| Sequence | S2 ignition 200.8 s, fairing 209.4 s, S2 cutoff 709.4 s, S2 relight 1761.5 s | CZ-8A typical sequence (huijiwiki) |
| Lift-off mass | **341.3 t** nominal (core throttled), T/W 1.43; 371.0 t at full thrust, the published 371 t | |
| S2 stage size | 12.375 m × 3.35 m | Wikipedia; no published propellant load found |

The flight law:

1. Vertical rise to T+10 s, then a pitch kick.
2. A zero-incidence gravity turn on an Earth-fixed great circle to t_q.
3. From t_q a linear-tangent pitch, starting at the flight path, to booster separation, continued without an attitude step through the core-only phase at a second rate. Without t_q (the zone-only fits of `validate.py`), the gravity turn runs to booster separation and the core pitch starts from the flight path.
4. S2 burn 1 with linear-tangent pitch and yaw in the target-plane frame, blended from the flight path over 10 s.
   - **Nominal (attached long coast, `--profile attached`):** Burn 1 is solved (`cz8a.solve_attached`: four steering tangents and the burn time) so that the transfer orbit's apogee lies where the inserted orbit passes u = 323°, at that moment and in its plane. A two-body apogee prediction is corrected by the integrated (J2) apogee time and radius and by the along-track miss at burn-2 ignition. Burn 2 is centred on apogee and ends on the inserted orbit's radial, horizontal and normal velocity, so the stack does not simply circularize. `build_traj.attach` keeps the vehicle on the element set (zero shift) when the S2 can pay for it, and otherwise solves the smallest shift ahead.
   - **Direct ascent (`--profile direct`; no solution with the 30.2 t S2):** burn 1 lasts the published 508.6 s, to T+709.4 s. The four steering tangents put the apex at the parking radius at the time that makes burn 2 light at T+1761.5 s. Of the two branches that do this, the lofted one is kept (cutoff 582 km). The flat one skims at 100–160 km for 500 s.
   - **Long coast (`--profile longcoast`, Y10 type):** burn 1 ends at a 200 km perigee with flight-path angle 0, followed by half a transfer orbit of coast.
5. A coast under J2–J4.
6. Burn 2 centred on the apex or apogee. For the direct and long-coast profiles its radius is iterated until the fitted SGP4 mean altitude equals Y12's first TLEs (1049.15 km).

**Nominal fit (revision 5).** Five parameters fit nine residuals: the zones' four along- and cross-track misses and the sightings' five angles.

- The parameters are kick, azimuth, t_q, and the booster- and core-phase pitch rates.
- The result is kick 0.395°, azimuth 139.24°, t_q 135.8 s, rates −1.06 and +0.51.
- Zone misses are 4.0 km or less; sightings miss by 0.35–0.77°.
- Five starting points find two minima. The other one, kick 0.49° with both rates at their bounds, misses by 3.1° and 54 km.

**Zone-only fits** (`validate.py`, revisions 1–4): three parameters (kick, core pitch rate, azimuth) for the four zone misses, at full thrust. Booster range against kick has a maximum, so the zones have two solutions. The flat branch put booster separation at 71 km in revision 4. The steep branch put it at 136 km in revision 3 (withdrawn). Between them, booster range rises to 1350 km. The sighting-fitted ascent sits between the two branches; it lands the boosters in zone 1 only because the throttled core takes about 0.3 km/s off.

S2 yaw takes the corridor's plane into the 50.0° target plane while keeping the node.

Other checks:

- Max-q is 26.8 kPa at T+69 s (Mach 1.30), with the core at 77.5 % (33.2 kPa in revision 4).
- **Mass closure (attached):** S2 uses 28.91 + 0.24 + 0.04 t, leaving **1.01 t of 30.2 t** for the disposal dump.
- **Sensitivity**, revision 2 direct profile at 26.0 t ([sensitivity.json](output/sensitivity.json); not rerun since):
  - With +3 t of payload the zones refit at kick 0.43° and 0.66 t remains.
  - S2 at 176 kN with 23.5 t, core β 200 and K 308.58 all keep the apex burn within 0.5° of 40.8°S 168.4°E and the lead within −13.9° to −14.7°.
  - Only a different relight time moves it: at 1950 s the burn is at 46.3°S 177.7°W.

## Visibility

Lift-off is at 19:59 Beijing time. The Sun is 30° below the horizon at booster separation and 40° below under S2 cutoff. **The whole powered ascent is in Earth's shadow: no sunlit plume during burn 1.**

The **≥10° zone** on the [regional map](output/CZ-8A_Y13_visibility.png) is therefore the **self-luminous engine flame**. It covers points that see the vehicle at ≥ 10° elevation between T+0 and S2 cutoff (T+828 s at 717 km), while their own Sun-centre altitude is ≤ 0° at that moment. It is clipped to the night side at the S2-cutoff epoch (12:13:47 UTC); the Sun −12° line at that epoch crosses the map's south-west corner.

Because the climb is steeper, the zone covers **15.7 million km²**, from 105.1–144.9°E and 18.9°S–24.3°N (13.5 million km² in revision 4). It takes in:

- Hainan, Guangzhou, Hong Kong and Kaohsiung;
- central and southern Vietnam;
- the Philippines, Borneo, Sulawesi, the Moluccas, Timor, Palau and the Top End of Australia.

Hanoi, Taipei, Singapore and Jakarta stay outside it.

The regional map is 1:1, extent 103–148°E, 22°S–27°N (3288 × 3960 px). A launch-corridor inset shows the observer, the four sightlines (observed azimuth out to the model's sub-vehicle point) and the staging. The zone uses the geometric Sun centre against the WGS84 ellipsoid with JPL DE440s, and leaves out brightness, extinction, weather and plume size. Data: [visibility.json](output/visibility.json), [illumination.csv](output/illumination.csv), [visibility_boundary.csv](output/visibility_boundary.csv).

**Attached long coast (nominal):** first sunlight comes during the coast at **T+30.0 min (12:29:59 UTC), 41.6°S 169.8°E, 891 km**, over New Zealand's South Island, with no engine running. The apogee burn at 12:54:25 UTC over the south-east Pacific (27.7°S 90.1°W) is in morning daylight. No sunlit plume is expected.

**Direct ascent, for reference (revision 2 values, T0 11:59:00):** its apex burn is sunlit over a dark New Zealand. First sunlight there comes at T+28.5 min, 34 s before the second ignition. The 27 s burn (12:28:21–12:28:49 UTC, 01:28 NZDT) is entirely sunlit while the ground below is near local midnight:

| Site | Sun altitude | Vehicle elevation | Azimuth |
|---|---:|---:|---:|
| Christchurch | −27° | 68° | 312° |
| Hokitika | −28° | 78° | 320° |
| Wellington | −29° | 63° | 268° |
| Auckland | −33° | 52° | 221° |
| Sydney | −34° | 20° | 123° |
| Hobart | −24° | 20° | 92° |

A sunlit LH2/LOX plume ("jelly") seen from New Zealand around 01:28 NZDT would therefore point to the direct profile; its absence fits the long coast. Brightness is not modeled.

## Files

The scripts live at the repository root and run from there. Their `output/` is this folder's [output/](output/).

```
~/.venvs/astro313/bin/python fetch_data.py     # Space-Track: Y12 history, 2025-295 current gp (cached)
~/.venvs/astro313/bin/python j2_plan.py        # K, D, T0 solutions -> output/j2_plan.json
~/.venvs/astro313/bin/python build_traj.py     # [--t0 ISO] [--payload kg] [--profile attached|direct|longcoast] [--shift s] [--u-apo deg] [--full-thrust] [--sens]
~/.venvs/astro313/bin/python validate.py       # zone-predicted K on Y7/Y10/Y11/Y12 -> output/validation.json
~/.venvs/astro313/bin/python validate.py phase # catalog phase lead per flight and profile -> output/phase_validation.json
~/.venvs/astro313/bin/python compare_outside.py # model vs outside-ref/outside.tle -> output/outside_comparison.json
~/.venvs/astro313/bin/python coast_tradeoff.py # phase shift and S2 left vs coast, T0 second, S2 load, full thrust -> output/coast_tradeoff.json
~/.venvs/astro313/bin/python longcoast_check.py # revision 2 check: perigee-cutoff long coasts from T0 11:59:00 -> output/longcoast_check.json
~/.venvs/astro313/bin/python make_map.py       # visibility zone + both PNGs (3960 px long side, 300 dpi, 1:1 aspect)
~/.venvs/astro313/bin/python web_export.py     # mission data for the 3D page -> docs/data/cz8a-y13.json
```

[sightings.py](../../sightings.py) holds the observer, the four sightings and their residuals; `cz8a.fit_lower(..., sight=True)` uses it.

- [CZ-8A_Y13_traj.csv](output/CZ-8A_Y13_traj.csv): lift-off to deployment. Phases are boost, core, staging, S2 burn 1, coast, S2 burn 2 and orbit. Columns: ECEF x/y/z (km), WGS84 lat/lon/alt, inertial and Earth-relative speed, flight-path angle, downrange.
- [CZ-8A_Y13_orbit.csv](output/CZ-8A_Y13_orbit.csv): the payload for two revolutions after deployment.
- [CZ-8A_Y13_stage.csv](output/CZ-8A_Y13_stage.csv): the stage coast, disposal burn and stage orbit.
- [CZ-8A_Y13_falls.csv](output/CZ-8A_Y13_falls.csv): fall trajectories of the boosters, core and fairing.
- [CZ-8A_Y13_traj.kml](output/CZ-8A_Y13_traj.kml): the 3D track, falls, impact points and NOTAM polygons.
- [CZ-8A_Y13_solution.json](output/CZ-8A_Y13_solution.json): all events, elements, steering, the sighting residuals, the budget and the assumptions.

## After launch

1. When the catalog has Y13, fetch its first TLEs with the space-track skill.
2. Compare i, Ω, mean altitude and the phase lead at T0 (`cz8a.lead_at_t0`) with the tables above.
3. Solve T0 seconds from Ω with K = 309.70°. The sighting-fitted corridor gives K 310.32° at S2 ignition, but the S2 yaws 2.4° after that, so the node pins the seconds only through Y12's K.
4. Compare the catalog's along-track position with the outside set. The model is on it (zero shift).
5. Rerun `build_traj.py --t0 <actual>`, `compare_outside.py` and `make_map.py`. In the attached profile T0 also changes the S2 burn, because the arrival time on the inserted orbit is fixed.
6. If timed sightings appear (video with a clock), add the times to `sightings.SIGHTINGS`. They would pin the climb rate and so test the throttle.

## Sources

- Element data: Space-Track.org.
- Launch and object catalog: GCAT (J. McDowell, planet4589.org), release 1.8.9.
- Y13 T0 19:59 BJT: [Xinhua](http://www.xinhuanet.com/20261120/da600aa6b4e94c90a93b2dd59f9081f3/c.html).
- Y12 T0 10:03 BJT: [CCTV](https://news.cctv.com/2026/10/25/ARTIM9NOyzQLvtahDjw0iMWc261025.shtml) and [China Daily / Xinhua](https://cn.chinadaily.com.cn/a/202610/25/WS6add9bb6e4b09a165c791524.html).
- Y10 hour-level LH2/LOX coast and direct insertion at 1100 km: [Beijing News](https://m.bjnews.com.cn/detail/1785983705129729.html).
- Flight sequence: [huijiwiki 长征八号甲](https://sat.huijiwiki.com/wiki/%E9%95%BF%E5%BE%81%E5%85%AB%E5%8F%B7%E7%94%B2) and [长征八号](https://sat.huijiwiki.com/wiki/%E9%95%BF%E5%BE%81%E5%85%AB%E5%8F%B7) (Y1 booster separation 171.88 s).
- Engines and stages: [Wikipedia, Long March 8](https://en.wikipedia.org/wiki/Long_March_8).
- CZ-8 throttle: [Baidu Baike 长征八号](https://baike.baidu.com/item/%E9%95%BF%E5%BE%81%E5%85%AB%E5%8F%B7/55468678).
- NOTAMs: RPHI B6120/26 and B6122/26; earlier flights from `../notam/notam_data.json`.
- Outside element set: `outside-ref/outside.tle`, supplied by the user (derived from Space-Track's 2026-243E).
- Ground sightings from N 18°26′14″, E 110°04′28″: supplied by the user (2026-11-21).
- Pad: GCAT `sites.tsv` (HCSLS) and the CZ-10B project's pad; GCAT `lp.tsv` lists LC1 at 19.5685°N, which is offshore.
