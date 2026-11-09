"""Rebuild the English figures used in the README from the processed CSVs.

Run from the repository root:  python scripts/make_readme_figures.py
Inputs : data/processed/*.csv, data/boundaries/nouakchott_moughataa.geojson,
         figures/gt_04c_sar_2016_vs_2025.png (map panels exported by the notebook)
Outputs: figures/fig1_*.png ... figures/fig4_*.png
"""
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA, FIG = ROOT / "data", ROOT / "figures"

# Palette (single accent + neutrals)
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SURFACE = "#fcfcfb"
SEQ = LinearSegmentedColormap.from_list(
    "blue_seq", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans", "font.size": 10.5,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False,
})


def area_weights(gdf):
    a = gdf.to_crs(32628).set_index("ADM2_EN").area
    return a / a.sum()


def sen_slope(t, y):
    s = [(y[j] - y[i]) / (t[j] - t[i]) for i in range(len(t)) for j in range(i + 1, len(t))]
    return float(np.median(s))


gdf = gpd.read_file(DATA / "boundaries" / "nouakchott_moughataa.geojson")
w = area_weights(gdf)
q = pd.read_csv(DATA / "processed" / "quarterly_uf_stl.csv", parse_dates=["date"])
annual = pd.read_csv(DATA / "processed" / "sar_urban_fraction_timeseries.csv").set_index("moughataa")
eco = pd.read_csv(DATA / "processed" / "eco_quarterly_tests.csv", parse_dates=["date"], index_col="date")

# ── Fig 1 — before / after radar panels ────────────────────────────────────
src = Image.open(FIG / "gt_04c_sar_2016_vs_2025.png").convert("RGB")
W, H = src.size
pw = W // 4
top = int(H * 0.11)
def trim_to_map(im, thr=235, frac=0.35):
    """Keep only rows/columns that are mostly map (drops frame, margins, in-figure legend)."""
    a = np.asarray(im).min(axis=2) < thr
    rows = np.where(a.mean(axis=1) > frac)[0]
    cols = np.where(a.mean(axis=0) > frac)[0]
    return im.crop((cols[0], rows[0], cols[-1] + 1, rows[-1] + 1))


def inner(im, k=6):
    w_, h_ = im.size
    return im.crop((k, k, w_ - k, h_ - k))


panels = [trim_to_map(inner(trim_to_map(src.crop((i * pw, top, (i + 1) * pw, H)))), frac=0.5)
          for i in (0, 1, 3)]
titles = ["Radar backscatter, 2016", "Radar backscatter, 2025",
          "Change 2016 → 2025 (red = new built-up)"]
fig, axes = plt.subplots(1, 3, figsize=(15, 7.4))
for ax, im, t in zip(axes, panels, titles):
    ax.imshow(im); ax.set_axis_off(); ax.set_title(t, loc="center", fontsize=12)
fig.suptitle("Nouakchott seen by Sentinel-1: built-up areas reflect strongly, sand stays dark",
             x=0.01, ha="left", fontsize=14, fontweight="bold", color=INK)
fig.text(0.01, 0.015, "Right panel: red = stronger radar return (new buildings), blue = weaker. Yellow lines: district (moughataa) boundaries. Dots: ground-truth calibration points.",
         color=INK2, fontsize=10)
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
fig.savefig(FIG / "fig1_radar_2016_vs_2025.png", dpi=110)
plt.close(fig)

# ── Fig 2 — growth map by district (Theil–Sen slope, pp/year) ──────────────
rows = []
for z, d in q.groupby("zone"):
    t = ((d.date - d.date.min()).dt.days / 365.25).to_numpy()
    y = d.uf_trend.to_numpy()
    tau, p = stats.kendalltau(t, y)
    rows.append({"ADM2_EN": z, "slope": sen_slope(t, y), "tau": tau, "p": p,
                 "start": y[0], "end": y[-1]})
growth = pd.DataFrame(rows)
growth.round(4).to_csv(DATA / "processed" / "district_growth_summary.csv", index=False)
m = gdf.merge(growth, on="ADM2_EN").to_crs(32628)

fig, ax = plt.subplots(figsize=(8.5, 8.5))
m.plot(column="slope", cmap=SEQ, vmin=0, vmax=2.3, edgecolor=SURFACE, linewidth=2, ax=ax)
for _, r in m.iterrows():
    c = r.geometry.representative_point()
    dark = r.slope > 1.2
    ax.annotate(f"{r.ADM2_EN}\n+{r.slope:.2f} pp/yr", (c.x, c.y), ha="center", va="center",
                fontsize=9.5, fontweight="bold", color="white" if dark else INK)
ax.set_axis_off()
ax.set_title("Where Nouakchott is growing fastest\nBuilt-up share growth per year, 2016–2026 (Theil–Sen slope)",
             loc="left", fontsize=13)
sm = plt.cm.ScalarMappable(cmap=SEQ, norm=plt.Normalize(0, 2.3))
cb = fig.colorbar(sm, ax=ax, shrink=0.55, pad=0.02)
cb.set_label("percentage points per year", color=INK2); cb.outline.set_visible(False)
ax.text(0, -0.02, "Mann–Kendall trend significant in all 9 districts (p < 0.01).",
        transform=ax.transAxes, color=INK2, fontsize=9.5)
fig.tight_layout()
fig.savefig(FIG / "fig2_growth_map_by_district.png", dpi=120)
plt.close(fig)

# ── Fig 3 — quarterly index: city + districts ──────────────────────────────
piv = q.pivot_table(index="date", columns="zone", values="uf_trend")
city = piv[w.index].mul(w, axis=1).sum(axis=1)
fig, ax = plt.subplots(figsize=(11, 5.6))
for z in piv.columns:
    ax.plot(piv.index, piv[z], color=AXIS, lw=1.3)
    if z in ("Arafat", "Dar Naim", "Ksar", "Riyadh"):
        ax.annotate(f"{z} {piv[z].iloc[-1]:.0f}%", (piv.index[-1], piv[z].iloc[-1]),
                    xytext=(6, 0), textcoords="offset points", va="center", color=INK2, fontsize=9.5)
ax.plot(city.index, city, color=BLUE, lw=3)
ax.annotate(f"Whole city {city.iloc[-1]:.0f}%", (city.index[-1], city.iloc[-1]), xytext=(6, 0),
            textcoords="offset points", va="center", color=INK, fontweight="bold", fontsize=10)
ax.set_xlim(piv.index[0], piv.index[-1] + pd.Timedelta(days=500))
ax.set_ylim(25, 100)
ax.set_ylabel("Built-up share of district (%)")
ax.set_title("Quarterly, seasonally adjusted urbanisation index, 2016 Q1 – 2026 Q2\n"
             "Blue: area-weighted city index · grey: the 9 districts")
fig.tight_layout()
fig.savefig(FIG / "fig3_quarterly_index.png", dpi=120)
plt.close(fig)

# ── Fig 4 — spurious in levels vs. real in changes ─────────────────────────
d = eco[["sar_sa", "icc_global"]].dropna()
z = (d - d.mean()) / d.std()
r_lvl = stats.pearsonr(d.sar_sa, d.icc_global)[0]
dd = d.diff().dropna()
r_dif = stats.pearsonr(dd.sar_sa, dd.icc_global)[0]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.8))
a1.plot(z.index, z.sar_sa, color=BLUE, lw=2, label="Satellite urban index")
a1.plot(z.index, z.icc_global, color=ORANGE, lw=2, label="Construction Cost Index (ANSADE)")
a1.set_ylabel("standardised (z-score)")
a1.set_title(f"Levels: both trend upward, r = {r_lvl:.2f}\n(a misleading correlation)")
a1.legend(loc="upper left", fontsize=9.5)
a2.scatter(dd.icc_global, dd.sar_sa, s=42, color=BLUE, edgecolor=SURFACE, linewidth=1.5)
a2.axhline(0, color=AXIS, lw=1); a2.axvline(0, color=AXIS, lw=1)
a2.set_xlabel("Quarterly change in construction cost index (points)")
a2.set_ylabel("Quarterly change in built-up share (pp)")
a2.set_title(f"Quarter-on-quarter changes: r = {r_dif:.2f}\n(the honest starting point for the VECM)")
fig.tight_layout()
fig.savefig(FIG / "fig4_spurious_vs_real_correlation.png", dpi=120)
plt.close(fig)

print(growth.sort_values("slope", ascending=False).round(2).to_string(index=False))
print(f"city: {city.iloc[0]:.1f}% -> {city.iloc[-1]:.1f}%  | r levels {r_lvl:.2f}  r diffs {r_dif:.2f}")
