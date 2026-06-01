import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx

from shapely.geometry import Polygon
from scipy.spatial import Voronoi
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# ======================================================
# PATHS
# ======================================================

BASE = r"C:\Users\Jona\Desktop\Urban_RF_EMF_Tirana"

INPUT_CSV = os.path.join(
    BASE, "data", "processed", "datasets", "urban_rf_emf_dataset.csv"
)

OUT_DIR = os.path.join(
    BASE, "results", "figures", "paper_voronoi_final_renamed"
)

os.makedirs(OUT_DIR, exist_ok=True)


# ======================================================
# SETTINGS
# ======================================================

ZONE_COL = "Zone"
LAT_COL = "Lat"
LON_COL = "Lon"

ZONE_RENAME = {
    "Zogu_i_Zi": "Sheshi Karl Topia",
    "Zogu i Zi": "Sheshi Karl Topia",
    "21_Dhjetori": "Sheshi Mustafa Qemal Ataturk",
    "21 Dhjetori": "Sheshi Mustafa Qemal Ataturk"
}

BANDS = [
    "E_total",
    "LTE800",
    "GSM900",
    "LTE1800",
    "UMTS2100",
    "LTE2600",
    "NR3500"
]

CRS_WGS84 = "EPSG:4326"
CRS_WEB = "EPSG:3857"

N_CLASSES = 5
INFLUENCE_RADIUS_M = 120
MAP_PAD_M = 120

COLORS = [
    "#eff3ff",
    "#bdd7e7",
    "#6baed6",
    "#3182bd",
    "#08519c"
]

CELL_ALPHA = 0.86
CELL_EDGE_COLOR = "white"
CELL_EDGE_WIDTH = 0.55
POINT_SIZE = 13
POINT_COLOR = "black"
BASEMAP_ALPHA = 0.68


# ======================================================
# HELPERS
# ======================================================

def clean_filename(text):
    return (
        text.replace(" ", "_")
        .replace("/", "_")
        .replace("ë", "e")
        .replace("ç", "c")
    )


def voronoi_finite_polygons_2d(vor, radius=None):
    new_regions = []
    new_vertices = vor.vertices.tolist()
    center = vor.points.mean(axis=0)

    if radius is None:
        radius = np.ptp(vor.points, axis=0).max() * 3

    all_ridges = {}

    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    for p1, region_index in enumerate(vor.point_region):
        vertices = vor.regions[region_index]

        if all(v >= 0 for v in vertices):
            new_regions.append(vertices)
            continue

        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1

            if v1 >= 0:
                continue

            tangent = vor.points[p2] - vor.points[p1]
            tangent = tangent / np.linalg.norm(tangent)

            normal = np.array([-tangent[1], tangent[0]])
            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, normal)) * normal

            far_point = vor.vertices[v2] + direction * radius

            new_vertices.append(far_point.tolist())
            new_region.append(len(new_vertices) - 1)

        vertices_array = np.asarray([new_vertices[v] for v in new_region])
        centroid = vertices_array.mean(axis=0)

        angles = np.arctan2(
            vertices_array[:, 1] - centroid[1],
            vertices_array[:, 0] - centroid[0]
        )

        new_region = np.array(new_region)[np.argsort(angles)]
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


def create_boundary(zone_gdf):
    point_buffers = zone_gdf.geometry.buffer(INFLUENCE_RADIUS_M)
    buffer_union = point_buffers.union_all()

    hull = zone_gdf.geometry.union_all().convex_hull.buffer(60)
    boundary = buffer_union.intersection(hull)

    return boundary


def create_voronoi_cells(zone_gdf):
    coords = np.array([[geom.x, geom.y] for geom in zone_gdf.geometry])

    if len(coords) < 4:
        raise ValueError("Voronoi kërkon të paktën 4 pika për zonë.")

    vor = Voronoi(coords)
    regions, vertices = voronoi_finite_polygons_2d(vor)

    boundary = create_boundary(zone_gdf)

    polygons = []

    for i, region in enumerate(regions):
        polygon = Polygon(vertices[region])
        clipped = polygon.intersection(boundary)

        if not clipped.is_empty:
            polygons.append({
                "point_id": i,
                "geometry": clipped
            })

    vor_gdf = gpd.GeoDataFrame(polygons, crs=CRS_WEB)

    point_data = zone_gdf.reset_index(drop=True).copy()
    point_data["point_id"] = point_data.index

    vor_gdf = vor_gdf.merge(
        point_data.drop(columns="geometry"),
        on="point_id"
    )

    return vor_gdf, boundary


def add_basemap(ax):
    ctx.add_basemap(
        ax,
        source=ctx.providers.Esri.WorldImagery,
        crs=CRS_WEB,
        attribution=False,
        zoom=17,
        alpha=BASEMAP_ALPHA
    )


def set_extent(ax, boundary):
    xmin, ymin, xmax, ymax = boundary.bounds
    ax.set_xlim(xmin - MAP_PAD_M, xmax + MAP_PAD_M)
    ax.set_ylim(ymin - MAP_PAD_M, ymax + MAP_PAD_M)


def add_scale_bar(ax, length=200):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    x0 = xlim[0] + (xlim[1] - xlim[0]) * 0.07
    y0 = ylim[0] + (ylim[1] - ylim[0]) * 0.07

    ax.plot(
        [x0, x0 + length],
        [y0, y0],
        color="black",
        linewidth=4,
        solid_capstyle="butt"
    )

    ax.text(
        x0 + length / 2,
        y0 + (ylim[1] - ylim[0]) * 0.022,
        f"{length} m",
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2)
    )


def add_north_arrow(ax):
    ax.annotate(
        "N",
        xy=(0.94, 0.87),
        xytext=(0.94, 0.77),
        xycoords="axes fraction",
        arrowprops=dict(
            facecolor="black",
            edgecolor="black",
            width=3,
            headwidth=12,
            headlength=12
        ),
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold"
    )


def make_bins(values):
    bins = np.quantile(values, np.linspace(0, 1, N_CLASSES + 1))
    bins = np.unique(bins)

    if len(bins) < N_CLASSES + 1:
        bins = np.linspace(values.min(), values.max(), N_CLASSES + 1)

    return bins


def make_labels(bins):
    return [f"{bins[i]:.4f} – {bins[i + 1]:.4f}" for i in range(len(bins) - 1)]


def plot_legend(ax, labels, band):
    handles = []

    for i, label in enumerate(labels):
        handles.append(
            Patch(
                facecolor=COLORS[i],
                edgecolor="gray",
                label=label,
                alpha=CELL_ALPHA
            )
        )

    handles.append(
        Line2D(
            [0],
            [0],
            marker="o",
            color="black",
            markerfacecolor="black",
            linestyle="None",
            markersize=6,
            label="Measurement point"
        )
    )

    ax.legend(
        handles=handles,
        title=f"{band} RF-EMF exposure (V/m)",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=3,
        fontsize=9,
        title_fontsize=10,
        frameon=True
    )


# ======================================================
# LOAD DATA
# ======================================================

df = pd.read_csv(INPUT_CSV)

df = df.dropna(subset=[LAT_COL, LON_COL, ZONE_COL])

df[ZONE_COL] = df[ZONE_COL].replace(ZONE_RENAME)

gdf = gpd.GeoDataFrame(
    df,
    geometry=gpd.points_from_xy(df[LON_COL], df[LAT_COL]),
    crs=CRS_WGS84
).to_crs(CRS_WEB)

zones = sorted(gdf[ZONE_COL].unique())

print("Detected zones:")
print(zones)


# ======================================================
# MAIN LOOP
# ======================================================

for band in BANDS:

    if band not in gdf.columns:
        print(f"Skipped {band}: column not found.")
        continue

    print(f"\nProcessing {band}...")

    gdf_band = gdf.dropna(subset=[band]).copy()
    gdf_band[band] = gdf_band[band].astype(float)

    bins = make_bins(gdf_band[band])
    labels = make_labels(bins)

    for zone in zones:

        zone_gdf = gdf_band[gdf_band[ZONE_COL] == zone].copy()

        if len(zone_gdf) < 4:
            print(f"Skipped {zone}: not enough points.")
            continue

        print(f"  Zone: {zone}")

        vor_gdf, boundary = create_voronoi_cells(zone_gdf)

        vor_gdf["class"] = pd.cut(
            vor_gdf[band],
            bins=bins,
            labels=labels,
            include_lowest=True
        )

        fig, ax = plt.subplots(figsize=(8, 8), dpi=300)

        set_extent(ax, boundary)
        add_basemap(ax)

        gpd.GeoSeries([boundary], crs=CRS_WEB).boundary.plot(
            ax=ax,
            color="black",
            linewidth=0.8,
            alpha=0.65,
            zorder=3
        )

        for i, label in enumerate(labels):
            subset = vor_gdf[vor_gdf["class"] == label]

            if not subset.empty:
                subset.plot(
                    ax=ax,
                    color=COLORS[i],
                    edgecolor=CELL_EDGE_COLOR,
                    linewidth=CELL_EDGE_WIDTH,
                    alpha=CELL_ALPHA,
                    zorder=4
                )

        zone_gdf.plot(
            ax=ax,
            color=POINT_COLOR,
            markersize=POINT_SIZE,
            alpha=0.95,
            zorder=5
        )

        ax.set_title(
            f"{zone} — {band}",
            fontsize=17,
            fontweight="bold",
            pad=10
        )

        add_scale_bar(ax, 200)
        add_north_arrow(ax)
        plot_legend(ax, labels, band)

        ax.set_axis_off()
        plt.tight_layout()

        safe_zone = clean_filename(zone)

        out_png = os.path.join(
            OUT_DIR,
            f"{safe_zone}_{band}_Voronoi_Final.png"
        )

        out_tiff = os.path.join(
            OUT_DIR,
            f"{safe_zone}_{band}_Voronoi_Final.tiff"
        )

        out_pdf = os.path.join(
            OUT_DIR,
            f"{safe_zone}_{band}_Voronoi_Final.pdf"
        )

        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.savefig(out_tiff, dpi=300, bbox_inches="tight")
        plt.savefig(out_pdf, bbox_inches="tight")

        plt.close()

        print("Saved:")
        print(out_png)
        print(out_tiff)
        print(out_pdf)


print("\nDONE.")
print("Figures saved in:")
print(OUT_DIR)
