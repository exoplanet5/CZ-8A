"""Maps for CZ-8A Y13, SatObserver-MX style (../../sda/starship-ITF-14/work/ift14_style.py).

    ~/.venvs/astro313/bin/python make_map.py   -> output/CZ-8A_Y13_visibility.png, output/CZ-8A_Y13_orbit.png

Both 300 dpi with a 3960 px long side, map axes 1:1 (equal degrees in longitude and latitude).  The regional map carries the ascent, the NOTAM zones and stage falls, the
10 deg zone of visibility.py, night shading and the Sun 0 / -12 deg lines at the S2 cutoff epoch.
"""
import csv
import importlib.util
import json
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import patheffects as pe  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Polygon, Rectangle  # noqa: E402

import cz8a as C  # noqa: E402
import visibility as V  # noqa: E402

spec = importlib.util.spec_from_file_location("ift14_style", Path("~/sda/starship-ITF-14/work/ift14_style.py").expanduser())
s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)
OUT = C.OUT
BLUE, ORANGE, GREEN, FALL, SIGHT = "#3987e5", "#d95926", "#199e70", "#c9a27a", "#c39bd3"
INK, DIM, BG = s.FG, s.FGDIM, s.BG
TRACK_HALO = [pe.Stroke(linewidth=3.2, foreground="#080b10"), pe.Normal()]
EXTENT = [103.0, 148.0, -22.0, 27.0]
INSET = [109.3, 113.3, 17.5, 20.3]           # launch corridor: pad, observer, sightlines, staging (1:1)
INSET_AX = [0.556, 0.469, 0.36]              # x, y, width in map-axes fraction (height follows from 1:1)
GLOBAL = [-180.0, 180.0, -70.0, 70.0]
CITIES = [("Haikou", 110.33, 20.03, (.3, .25, "left")), ("Sanya", 109.51, 18.25, (-.3, -.55, "right")),
          ("Hong Kong", 114.17, 22.32, (.3, .2, "left")), ("Nanning", 108.37, 22.82, (-.3, -.75, "right")),
          ("Hanoi", 105.85, 21.03, (-.3, .25, "right")), ("Da Nang", 108.2, 16.05, (.3, .25, "left")),
          ("Ho Chi Minh City", 106.7, 10.78, (.3, .25, "left")), ("Phnom Penh", 104.92, 11.56, (-.3, .25, "right")),
          ("Kaohsiung", 120.3, 22.6, (.3, .15, "left")), ("Manila", 120.98, 14.6, (.3, .25, "left")),
          ("Cebu", 123.9, 10.3, (.3, .25, "left")), ("Davao", 125.6, 7.07, (.3, .25, "left")),
          ("Zamboanga", 122.08, 6.92, (.3, -.65, "left")), ("Puerto Princesa", 118.74, 9.74, (-.3, .3, "right")),
          ("Kota Kinabalu", 116.07, 5.98, (-.3, .25, "right")), ("Kuching", 110.35, 1.55, (.3, .25, "left")),
          ("Singapore", 103.82, 1.35, (.3, .25, "left")), ("Manado", 124.84, 1.47, (.3, .25, "left")),
          ("Makassar", 119.42, -5.14, (.3, .25, "left")), ("Balikpapan", 116.83, -1.24, (-.3, .25, "right")),
          ("Ambon", 128.18, -3.7, (.3, .25, "left")), ("Jakarta", 106.85, -6.21, (.3, .25, "left")),
          ("Koror", 134.48, 7.34, (.3, .25, "left")), ("Darwin", 130.84, -12.46, (.3, .25, "left")),
          ("Dili", 125.57, -8.56, (.3, .25, "left")), ("Guangzhou", 113.26, 23.13, (-.3, .25, "right")),
          ("Taipei", 121.56, 25.04, (.3, .2, "left")), ("Xiamen", 118.09, 24.48, (-.3, .25, "right")),
          ("Kunming", 102.71, 25.04, (.3, .25, "left")), ("Changsha", 112.94, 28.23, (.3, .25, "left")),
          ("Wuhan", 114.31, 30.59, (.3, .25, "left")), ("Shanghai", 121.47, 31.23, (.3, .25, "left")),
          ("Bangkok", 100.50, 13.75, (.3, .25, "left")), ("Naha", 127.68, 26.21, (.3, .25, "left")),
          ("Yangon", 96.20, 16.84, (.3, .25, "left"))]


def read(path):
    return V.read(path)


def dest(p, brg, d_km):
    """Great-circle destination from p=(lat, lon) on bearing brg (deg) after d_km (spherical)."""
    la, lo, b, a = np.radians(p[0]), np.radians(p[1]), np.radians(brg), np.asarray(d_km) / C.RE
    la2 = np.arcsin(np.sin(la) * np.cos(a) + np.cos(la) * np.sin(a) * np.cos(b))
    lo2 = lo + np.arctan2(np.sin(b) * np.sin(a) * np.cos(la), np.cos(a) - np.sin(la) * np.sin(la2))
    return np.degrees(la2), np.degrees(lo2)


def track(ax, lon, lat, color, ls="-", lw=1.6, alpha=1.0):
    s.plot_track(ax, lon, lat, color=color, lw=lw, ls=ls, alpha=alpha, zorder=5, path_effects=TRACK_HALO)


def dotlabel(ax, x, y, text, dx=0, dy=0, marker="o", color=INK, size=9, ha="left", ms=None):
    ax.plot(x, y, marker, ms=ms or (5 if marker == "o" else 7.5), mfc=color, mec=BG, mew=1.2, zorder=7)
    ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points", fontsize=size, color=INK, ha=ha,
                va="center", zorder=8, path_effects=[s.HALO],
                arrowprops=dict(arrowstyle="-", color=DIM, linewidth=.7) if abs(dx) + abs(dy) > 18 else None)


def equal_figure(extent, title, subtitle, top=0.65, bottom=0.9, width=13.2):
    """Figure sized so the map axes are 1:1 (1 deg longitude = 1 deg latitude on the page, plate carree), with the
    long side of the figure `width` inches (3960 px at 300 dpi) in either orientation."""
    ratio = (extent[3] - extent[2]) / (extent[1] - extent[0])
    aw = 0.976 * width
    ah = aw * ratio
    h = ah + top + bottom
    if h > width:                      # portrait: the height is the long side
        h = width
        ah = h - top - bottom
        width = ah / ratio / 0.976
    fig, ax = s.new_figure(width, h, title, subtitle, "")
    fig.texts[0].set_y(1 - 0.32 / h)
    fig.texts[1].set_y(1 - 0.32 / h)
    ax.set_position([.012, bottom / h, .976, ah / h])
    return fig, ax


def lock_aspect(ax, extent):
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect("equal", adjustable="box")


def corridor_inset(ax, sol, tr, t0):
    """Upper-stage-free close-up of the launch corridor: pad, observer, the four sightlines (observed azimuth out to the
    model's sub-vehicle point at the matched time), booster separation and S2 ignition."""
    x, y, w = INSET_AX
    h = w * (EXTENT[1] - EXTENT[0]) / (EXTENT[3] - EXTENT[2]) * (INSET[3] - INSET[2]) / (INSET[1] - INSET[0])
    ins = ax.inset_axes([x, y, w, h], zorder=9)
    for sp in ins.spines.values():
        sp.set_edgecolor("#9aa4ae")
    s.draw_base(ins, INSET, t0.isoformat(), grat=1, label_step=1, night=False, sun=False)
    ins.add_patch(Rectangle((0, 0), 1, 1, transform=ins.transAxes, color=BG, alpha=.22, zorder=2))
    ax.add_patch(Rectangle((INSET[0], INSET[2]), INSET[1] - INSET[0], INSET[3] - INSET[2], fill=False,
                           edgecolor="#9aa4ae", lw=0.8, zorder=6))
    ins.text(.025, .955, "Launch corridor", transform=ins.transAxes, color=INK, fontsize=8.5, va="top",
             path_effects=[s.HALO], zorder=9)
    ph, t = tr["phase"], tr["t"]
    low = np.isin(ph, ["boost", "core", "staging"])
    b1 = (ph == "S2 burn 1") & (t < 260)
    track(ins, tr["lon"][low], tr["lat"][low], ORANGE, lw=2.0)
    j = np.flatnonzero(b1)
    track(ins, tr["lon"][np.r_[j[0] - 1, j]], tr["lat"][np.r_[j[0] - 1, j]], GREEN, lw=2.0)
    for name, lo, la, (dx, dy, ha) in CITIES:
        if INSET[0] < lo < INSET[1] and INSET[2] < la < INSET[3]:
            ins.plot(lo, la, "o", ms=3.5, mfc="#c9ccd1", mec=BG, mew=.8, zorder=7)
            s.label(ins, lo + .07, la - .02, name, color=INK, size=8, ha="left", va="center", clip_on=True)
    sg = sol["sightings"]
    obs = (sg["observer"]["lat"], sg["observer"]["lon"])
    for d in sg["fit"]:
        sla, slo = np.interp(d["t"], t, tr["lat"]), np.interp(d["t"], t, tr["lon"])
        la, lo = dest(obs, d["az_obs"], np.linspace(0.0, float(C.gc_dist(obs, (sla, slo))), 60))
        ins.plot(lo, la, color=SIGHT, lw=1.2, ls=(0, (2, 2)), zorder=5.5, path_effects=TRACK_HALO)
        ins.plot(slo, sla, "o", ms=6, mfc="none", mec=SIGHT, mew=1.4, zorder=7.5)
    ev = sol["events"]
    dotlabel(ins, C.PAD[1], C.PAD[0], "LC1  T+0", dx=10, dy=4, marker="^", color=ORANGE, size=8.5)
    dotlabel(ins, obs[1], obs[0], "Observer", dx=-8, dy=8, marker="D", color=SIGHT, size=8.5, ha="right", ms=5.5)
    e = ev["booster_sep"]
    dotlabel(ins, e["lon"], e["lat"], f"Boosters off  {e['alt_km']:.0f} km", dx=10, dy=6, marker="s", color=ORANGE,
             size=8.5)
    e = ev["s2_ignition"]
    dotlabel(ins, e["lon"], e["lat"], f"S2 ignition  {e['alt_km']:.0f} km", dx=-6, dy=-14, marker="s", color=GREEN,
             size=8.5, ha="right")
    ins.set_xlim(INSET[0], INSET[1])
    ins.set_ylim(INSET[2], INSET[3])
    ins.set_aspect("equal", adjustable="box")


def regional(sol, tr, fl, vis):
    t0 = datetime.fromisoformat(sol["t0"])
    fig, ax = equal_figure(EXTENT, "CZ-8A Launch visibility", t0.strftime("%Y-%m-%d %H:%M UTC"))
    fig.texts[1].set_x(.245)
    s.draw_base(ax, EXTENT, t0.isoformat(), grat=10, label_step=10, sun=False, night=False)
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, color=BG, alpha=.16, zorder=2))
    ref = vis["ref_sun"]
    shade = np.zeros(ref.shape + (4,))
    shade[..., 2] = 10 / 255
    shade[..., 3] = np.where(ref <= 0, .35, 0.0)
    ax.imshow(shade, extent=EXTENT, origin="lower", aspect="auto", interpolation="bilinear", zorder=2.1)
    ax.contourf(vis["lon"], vis["lat"], vis["field"], levels=[0.0, 100.0], colors=[INK], alpha=.10, zorder=3.2)
    bnd = ax.contour(vis["lon"], vis["lat"], vis["field"], levels=[0.0], colors=[INK], linewidths=1.1,
                     linestyles="--", zorder=4)
    lines = []
    for lev, col in ((0.0, s.SUN), (-12.0, s.ACCENT)):
        if ref.min() < lev < ref.max():
            c = ax.contour(vis["lon"], vis["lat"], ref, levels=[lev], colors=[col], linewidths=1.8, linestyles="-",
                           zorder=4)
            c.set_path_effects([pe.Stroke(linewidth=3.0, foreground=BG), pe.Normal()])
            lines.append((lev, c))
    for lev, c in lines:
        seg = max(c.allsegs[0], key=len)
        p = seg[np.argmin(abs(seg[:, 1] - 12.0))]
        ax.annotate("Sun −12°" if lev < 0 else "Sun 0°", p, xytext=(6, 0), textcoords="offset points", ha="left",
                    va="center", color=INK, fontsize=8.5, path_effects=[s.HALO], zorder=8)
    with open(OUT / "visibility_boundary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["segment", "lon_deg", "lat_deg"])
        for k, pts in enumerate(bnd.allsegs[0]):
            w.writerows([[k, f"{x:.3f}", f"{y:.3f}"] for x, y in pts])
    # NOTAM zones
    zlab = {0: "Boosters", 1: "Core + fairing"}
    for k, z in enumerate(sol["zones"]):
        ring = np.array([[lo, la] for la, lo in C.load_zones()[k]["ring"]])
        ax.add_patch(Polygon(ring, facecolor=(1, 1, 1, .10), edgecolor=INK, lw=1, zorder=6))
        la, lo = z["centroid"]
        s.label(ax, lo + (1.0 if k == 0 else -1.0), la + (0.55 if k == 0 else -0.35),
                f"{z['id']}  {zlab[k]}\n{z['b'][-4:]}–{z['c'][-4:]} UTC", size=8.5, va="center",
                ha="left" if k == 0 else "right")
    # stage falls
    for key in ("boosters", "core", "fairing"):
        m = fl["phase"] == key
        track(ax, fl["lon"][m], fl["lat"][m], FALL, lw=1.0, alpha=.9)
        f = sol["falls"][key]
        ax.plot(f["lon"], f["lat"], "x", ms=7, mew=1.8, color=FALL, zorder=7, path_effects=TRACK_HALO)
    # ascent and coast
    ph, t = tr["phase"], tr["t"]
    lower = np.isin(ph, ["boost", "core", "staging"])
    s2 = ph == "S2 burn 1"
    coast = (ph == "coast") & (t < 1500)
    track(ax, tr["lon"][lower], tr["lat"][lower], ORANGE)
    track(ax, tr["lon"][s2 | (np.r_[False, s2[:-1]])], tr["lat"][s2 | (np.r_[False, s2[:-1]])], GREEN)
    j0 = np.flatnonzero(s2)[-1]
    cm = coast.copy()
    cm[j0] = True
    track(ax, tr["lon"][cm], tr["lat"][cm], BLUE, ls=(0, (5, 3)))
    ev = sol["events"]
    dotlabel(ax, C.PAD[1], C.PAD[0], "Hainan LC1  T+0", dx=12, dy=3, marker="^", color=ORANGE)
    for key, col in (("booster_sep", ORANGE), ("s2_ignition", GREEN)):
        ax.plot(ev[key]["lon"], ev[key]["lat"], "s", ms=6, mfc=col, mec=BG, mew=1.0, zorder=7)
    cut = ev["s2_cutoff_1"]
    dotlabel(ax, cut["lon"], cut["lat"], f"S2 cutoff  T+{cut['t']:.0f} s\n{cut['alt_km']:.0f} km · "
             f"{(t0 + timedelta(seconds=cut['t'])):%H:%M:%S} UTC", dx=-16, dy=-30, marker="s", color=GREEN,
             ha="right", size=8.5)
    for tm in (300,):
        j = np.argmin(abs(t - tm))
        dotlabel(ax, tr["lon"][j], tr["lat"][j], f"T+{tm // 60} min", 9, 8, size=8.5)
    j = np.argmin(abs(t - (cut["t"] + 90)))
    s.label(ax, tr["lon"][j] + .5, tr["lat"][j] + .3, "Coast to apogee", size=8.5)
    s.label(ax, 126.0, 17.0, f"≥{V.ELEVATION_DEG:.0f}° ZONE", size=9, ha="center")
    # cities
    for name, lo, la, (dx, dy, ha) in CITIES:
        if not (EXTENT[0] < lo < EXTENT[1] and EXTENT[2] < la < EXTENT[3]):
            continue
        ax.plot(lo, la, "o", ms=3.5, mfc="#c9ccd1", mec=BG, mew=.8, zorder=7)
        s.label(ax, lo + dx, la + dy, name, color=INK, size=8, ha=ha, va="bottom", clip_on=True)
    corridor_inset(ax, sol, tr, t0)
    tro, pk = sol["transfer_orbit"], sol["parking_orbit_osculating"]
    s.hud(ax, [("CZ-8A Y13 · GUOWANG", INK, True),
               ((f"Transfer {tro['hp']:.0f} × {tro['ha']:.0f} km" if tro['hp'] > 0 else
                 f"Cutoff {cut['alt_km']:.0f} km · apex {tro['ha']:.0f} km"), INK, False),
               (f"Apogee burn T+{sol['apogee_burn']['t_ign'] / 60:.0f} min", INK, False),
               (f"Orbit {pk['ha']:.0f} km · i 50.0° · Ω {sol['plane']['raan_deg']:.1f}°", DIM, False)],
          x0=.725, y0=.875, lh=.026, width=.265, fontsize=9.5)
    handles = [Line2D([], [], color=ORANGE, lw=2, label="Boosters + core"),
               Line2D([], [], color=GREEN, lw=2, label="S2 burn 1"),
               Line2D([], [], color=BLUE, lw=2, ls=(0, (5, 3)), label="Coast"),
               Line2D([], [], color=FALL, lw=1.2, marker="x", ms=6, mew=1.6, label="Stage fall"),
               Patch(facecolor=(1, 1, 1, .10), edgecolor=INK, label="NOTAM zone"),
               Line2D([], [], color=INK, lw=1.2, ls="--", label=f"≥{V.ELEVATION_DEG:.0f}° zone"),
               Line2D([], [], color=SIGHT, lw=1.2, ls=(0, (2, 2)), marker="o", ms=5, mfc="none", mew=1.2,
                      label="Sightline")]
    for lev, _ in lines:
        handles.append(Line2D([], [], color=s.ACCENT if lev < 0 else s.SUN, lw=1.8,
                              label="Sun −12°" if lev < 0 else "Sun 0°"))
    leg = ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(.005, .012), ncols=int(np.ceil(len(handles) / 2)),
                    fontsize=8.5, facecolor=BG, edgecolor="#37404a", framealpha=.92)
    for txt in leg.get_texts():
        txt.set_color(INK)
    fig.text(.988, .022, "Illustration: Mickey  ·  base map NASA Blue Marble", color=INK, fontsize=7.5, ha="right",
             va="bottom")
    lock_aspect(ax, EXTENT)
    fig.savefig(OUT / "CZ-8A_Y13_visibility.png", dpi=300, facecolor=BG)
    plt.close(fig)


def global_track(sol, tr, orb):
    t0 = datetime.fromisoformat(sol["t0"])
    fig, ax = equal_figure(GLOBAL, "CZ-8A Y13 Ascent and orbit", t0.strftime("%Y-%m-%d %H:%M UTC"), bottom=0.5)
    fig.texts[1].set_x(.31)
    s.draw_base(ax, GLOBAL, (t0 + timedelta(seconds=sol["apogee_burn"]["t_ign"])).isoformat(),
                grat=30, label_step=30, sun=False, night=True)
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, color=BG, alpha=.12, zorder=2))
    ph, t = tr["phase"], tr["t"]
    for mask, col, ls in ((np.isin(ph, ["boost", "core", "staging"]), ORANGE, "-"), (ph == "S2 burn 1", GREEN, "-"),
                          (ph == "coast", BLUE, (0, (5, 3))), (ph == "S2 burn 2", GREEN, "-")):
        track(ax, tr["lon"][mask], tr["lat"][mask], col, ls=ls, lw=1.8 if col == GREEN else 1.5)
    rev = orb["t"] <= orb["t"][0] + sol["parking_orbit_osculating"]["period"]
    track(ax, orb["lon"][rev], orb["lat"][rev], "#9aa4ae", lw=1.1, alpha=.9)
    ev = sol["events"]
    dotlabel(ax, C.PAD[1], C.PAD[0], "Hainan", dx=-8, dy=10, marker="^", color=ORANGE, ha="right", size=8.5)
    cut = ev["s2_cutoff_1"]
    dotlabel(ax, cut["lon"], cut["lat"], f"S2 cutoff  T+{cut['t'] / 60:.1f} min", dx=10, dy=-12, marker="s",
             color=GREEN, size=8.5)
    b2 = ev["s2_ignition_2"]
    dotlabel(ax, b2["lon"], b2["lat"], f"Apogee burn  T+{b2['t'] / 60:.1f} min\n"
             f"{(t0 + timedelta(seconds=b2['t'])):%H:%M} UTC · {b2['alt_km']:.0f} km", dx=16, dy=-22, marker="s",
             color=GREEN, size=8.5, ha="left")
    p = json.loads((OUT / "visibility.json").read_text())["sunlight_pickup"]
    ax.plot(p["lon"], p["lat"], "*", ms=11, mfc=INK, mec=BG, mew=1, zorder=7)
    ax.annotate(f"Sunlight  T+{p['t_s'] / 60:.1f} min", (p["lon"], p["lat"]), xytext=(-12, 10), textcoords="offset points", ha="right",
                fontsize=8.5, color=INK, path_effects=[s.HALO], zorder=8)
    for tm in (1200, 6000):
        src = tr if tm <= t.max() else orb
        j = np.argmin(abs(src["t"] - tm))
        dotlabel(ax, src["lon"][j], src["lat"][j], f"T+{tm // 60} min", 7, 8, size=8)
    handles = [Line2D([], [], color=ORANGE, lw=2, label="Boosters + core"),
               Line2D([], [], color=GREEN, lw=2, label="S2 burns"),
               Line2D([], [], color=BLUE, lw=2, ls=(0, (5, 3)), label="Coast"),
               Line2D([], [], color="#9aa4ae", lw=1.4, label="First orbit"),
               Line2D([], [], color=INK, marker="*", ms=10, lw=0, label="Sunlight"),
               Line2D([], [], color=(1, 205 / 255, 110 / 255, .6), lw=1.2, label="Terminator")]
    leg = ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(.36, .015), ncols=3, fontsize=8.5,
                    facecolor=BG, edgecolor="#37404a", framealpha=.92)
    for txt in leg.get_texts():
        txt.set_color(INK)
    fig.text(.988, .025, "Illustration: Mickey  ·  base map NASA Blue Marble", color=INK, fontsize=7.5, ha="right",
             va="bottom")
    lock_aspect(ax, GLOBAL)
    fig.savefig(OUT / "CZ-8A_Y13_orbit.png", dpi=300, facecolor=BG)
    plt.close(fig)


def main():
    sol = json.loads((OUT / "CZ-8A_Y13_solution.json").read_text())
    tr = read(OUT / "CZ-8A_Y13_traj.csv")
    fl = read(OUT / "CZ-8A_Y13_falls.csv")
    orb = read(OUT / "CZ-8A_Y13_orbit.csv")
    vis = V.compute(extent=EXTENT)
    regional(sol, tr, fl, vis)
    global_track(sol, tr, orb)
    print("wrote", OUT / "CZ-8A_Y13_visibility.png", OUT / "CZ-8A_Y13_orbit.png")


if __name__ == "__main__":
    main()
