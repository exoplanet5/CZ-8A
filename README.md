# CZ-8A launch profiles: models and 3D page

Ascent-to-orbit models of Long March 8A (长征八号甲) launches from the Hainan commercial launch site, with an interactive 3D globe page. They are built from the NOTAM drop zones, ground sightings and catalog orbits. They are models, not telemetry and not official ephemerides.

**3D page: https://exoplanet5.github.io/CZ-8A/**

| Mission | Lift-off (UTC) | Page profile | Orbit modelled | Notes | 3D page |
|---|---|---|---|---|---|
| [CZ-8A Y13, SatNet LEO group 30](missions/Y13_20261120/README.md) | 2026-11-20 11:59 (Xinhua: 19:59 BJT; model 11:59:59) | SatNet · 50° | 1049 km mean, i 50.00° | Flown. Ascent fitted to four ground sightings; 40.6 min long coast to Y12's inserted orbit | [open](https://exoplanet5.github.io/CZ-8A/?profile=satnet-50&mission=cz8a-y13) |

Y12 (2026-243) flew the same profile as Y13. Its first element set, moved in plane to Y13, is Y13's inserted orbit, so it has no separate model.

Each mission folder holds its NOTAM text, inputs, results, maps and a README with the assumptions, uncertainty and sources.

## Layout

```
cz8a.py              CZ-8A vehicle and 3-DOF ascent model; the mission folder paths (MISSION, OUT, NOTAM)
sightings.py         Y13 ground sightings and their residuals
build_traj.py        trajectory, transfer and parking orbits, falls, TLEs and KML
j2_plan.py           site constant K and the J2 drift plan toward the target plane
validate.py          zone-predicted K and catalog phase on earlier CZ-8A flights
compare_outside.py   model against the inserted element set
coast_tradeoff.py    phase shift and S2 margin against coast length, T0 second and S2 load
visibility.py        10° visibility zone of the powered ascent
make_map.py          visibility and orbit maps
web_export.py        mission data for the 3D page
fetch_data.py        Space-Track queries -> data/
longcoast_check.py   revision 2 check, kept for the record
data/                Space-Track element sets: Y12 (2026-243) history, target plane 2025-295
missions/<flight>/   NOTAM, inputs, output/, maps and README of one flight
docs/                the 3D page (GitHub Pages serves main -> /docs)
tools/serve.py       local server for the page
```

## Model in brief

- **Vehicle:** core stage plus two K2 liquid boosters, and a 3.35 m hydrolox second stage with a 5.2 m fairing. Burn times follow the published CZ-8A sequence. Stage loads are consumption at those burn times. The S2 load (30.2 t) is assumed.
- **Lower stages:**
  - The model flies the published sequence.
  - It is fitted to the centroids of the two RPHI NOTAM zones: zone 1 for the boosters, zone 2 for the core and fairing.
  - For Y13 it is also fitted to four sightings from near Lingshui, which it matches to 0.77°.
  - The zones and the sightings fit together only with a core throttled to 77.5 % from T+10 s to booster separation. Lift-off is then 341 t, not the published 371 t. This is the open point of the model.
- **Second stage:**
  - Burn 1 ends on a 679 × 1061 km transfer orbit.
  - A 40.6 min long coast follows.
  - A 5 s apogee burn then lands on the inserted orbit with zero phase shift.
  - A propellant dump disposes of the stage.
- **Checks:** the zones alone predict the site constant K (the orbit's node relative to sidereal time at lift-off) of Y7, Y10, Y11 and Y12 to within 0.5° of the catalog. No TLE is used for that prediction.

## Running

The outputs are committed, so reading the results or the page needs nothing installed. To serve the page locally:

```sh
python3 tools/serve.py        # then open http://127.0.0.1:8766/
```

Rerunning the model needs Python 3.13 with numpy, scipy, matplotlib, skyfield, sgp4 and astropy. It also needs shared modules that are not in this repository:

- `~/launch/launch_traj_NROL-97/nrol97.py` (frames, gravity, atmosphere, coasts, elements);
- `~/launch/launch_traj_CERES-2/ceres2.py` (guided-burn kernel, pad geometry, NOTAM parser);
- `~/launch/launch_traj_12/j2_tuning.py` (J2 plan);
- `~/sda/starship-ITF-14/work/ift14_style.py` (map style);
- the JPL DE440s kernel at `~/JPLDE/de440s.bsp`.

The pipeline order is in [CLAUDE.md](CLAUDE.md) and in each mission README.

## Adding a mission

1. Make `missions/<flight>_<yyyymmdd>/` with the NOTAM file and inputs.
2. The scripts serve one flight through `cz8a.MISSION`. With a second flight, turn that into a mission key, as `build_traj.MISSIONS` does in the CZ-2D repository.
3. Export the page data with `web_export.py`. Register it with the launch-page skill's `new_page.py --add`, under an existing profile or a new one such as a different inclination.

## Sources and licence

- NOTAMs: RPHI (Manila FIR).
- Launch facts: Xinhua and CCTV.
- Catalog orbits: Space-Track.org.
- Launch site and vehicle: GCAT (J. McDowell, planet4589.org/space/gcat), Wikipedia and huijiwiki.

Per-mission sources are listed in each mission README. The code is under the MIT licence ([LICENSE](LICENSE)). The page's viewer, map imagery, borders and star data carry their own licences in `docs/LICENSE-CZ-6A-viewer.txt`, `docs/tex/TILES-LICENSE.txt` and `docs/data/*-LICENSE.txt`.
