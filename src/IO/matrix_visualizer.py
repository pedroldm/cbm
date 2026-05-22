"""
Matrix Analysis Interactive Visualizer
=======================================
Usage: python matrix_visualizer.py [path/to/data.json]
       (defaults to matrix_data.json in the same directory)
"""

import json
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.widgets import Button, RadioButtons, Slider
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import LineCollection
import warnings
warnings.filterwarnings("ignore")

# ── Palette ────────────────────────────────────────────────────────────────────
BG        = "#0e1117"
PANEL     = "#161b27"
ACCENT1   = "#4fc3f7"   # sky blue
ACCENT2   = "#f06292"   # rose
ACCENT3   = "#aed581"   # lime
ACCENT4   = "#ffb74d"   # amber
ACCENT5   = "#ce93d8"   # lavender
GRID_CLR  = "#1e2535"
TEXT_CLR  = "#cdd6f4"
DIM_CLR   = "#6272a4"
CMAP_ONES = LinearSegmentedColormap.from_list("ones", ["#161b27", ACCENT1])
CMAP_DENSITY = LinearSegmentedColormap.from_list(
    "density", ["#0e1117", "#1a2a4a", ACCENT1, "#ffffff"])

ACCENTS = [ACCENT1, ACCENT2, ACCENT3, ACCENT4, ACCENT5]

plt.rcParams.update({
    "figure.facecolor":  BG,
    "axes.facecolor":    PANEL,
    "axes.edgecolor":    GRID_CLR,
    "axes.labelcolor":   TEXT_CLR,
    "axes.titlecolor":   TEXT_CLR,
    "xtick.color":       DIM_CLR,
    "ytick.color":       DIM_CLR,
    "text.color":        TEXT_CLR,
    "grid.color":        GRID_CLR,
    "grid.linewidth":    0.6,
    "font.family":       "monospace",
    "axes.titlesize":    10,
    "axes.labelsize":    8,
    "xtick.labelsize":   7,
    "ytick.labelsize":   7,
})

# ── Data loading ───────────────────────────────────────────────────────────────

def load_data(path: str) -> list[dict]:
    with open(path) as f:
        return json.load(f)


def density(d):
    total = d["num_rows"] * d["num_cols"]
    return d["total_ones"] / total if total else 0


def singleton_row_ratio(d):
    return d["singleton_rows"] / d["num_rows"] if d["num_rows"] else 0


def singleton_col_ratio(d):
    return d["singleton_cols"] / d["num_cols"] if d["num_cols"] else 0


# ── Main visualizer ────────────────────────────────────────────────────────────

class MatrixViz:

    def __init__(self, data: list[dict]):
        self.data  = data
        self.names = [d["filename"] for d in data]
        self.sel   = 0        # currently selected file index
        self.view  = "overview"  # "overview" | "detail"

        self._build_figure()
        self._draw_overview()
        plt.show()

    # ── Layout ────────────────────────────────────────────────────────────────

    def _build_figure(self):
        self.fig = plt.figure(figsize=(17, 10), facecolor=BG)
        self.fig.canvas.manager.set_window_title("Matrix Analysis Visualizer")

        # title bar
        self.fig.text(0.5, 0.972, "◈  MATRIX  ANALYSIS  DASHBOARD  ◈",
                      ha="center", va="top", fontsize=14, color=ACCENT1,
                      fontweight="bold", fontfamily="monospace",
                      path_effects=[pe.withStroke(linewidth=6, foreground=BG)])
        self.fig.text(0.5, 0.952, f"  {len(self.data)} matrices loaded",
                      ha="center", va="top", fontsize=8, color=DIM_CLR)

        # two nav buttons
        ax_ov = self.fig.add_axes([0.01, 0.945, 0.08, 0.028])
        ax_dt = self.fig.add_axes([0.10, 0.945, 0.08, 0.028])
        self.btn_ov = Button(ax_ov, "Overview", color=PANEL, hovercolor=GRID_CLR)
        self.btn_dt = Button(ax_dt, "Detail",   color=PANEL, hovercolor=GRID_CLR)
        for btn in (self.btn_ov, self.btn_dt):
            btn.label.set_color(ACCENT1)
            btn.label.set_fontfamily("monospace")
            btn.label.set_fontsize(8)
        self.btn_ov.on_clicked(lambda _: self._switch("overview"))
        self.btn_dt.on_clicked(lambda _: self._switch("detail"))

        # file selector radio (left strip)
        radio_h = min(0.55, len(self.data) * 0.055 + 0.05)
        ax_radio = self.fig.add_axes([0.005, 0.38, 0.085, radio_h])
        ax_radio.set_facecolor(PANEL)
        self.radio = RadioButtons(ax_radio, self.names, activecolor=ACCENT1)
        for lbl in self.radio.labels:
            lbl.set_color(TEXT_CLR)
            lbl.set_fontfamily("monospace")
            lbl.set_fontsize(8)
        self.radio.on_clicked(self._on_radio)

        self.fig.text(0.045, 0.38 + radio_h + 0.005, "SELECT FILE",
                      ha="center", va="bottom", fontsize=7, color=DIM_CLR)

        self.content_axes = []   # axes we clear on view switch

    def _clear_content(self):
        for ax in self.content_axes:
            self.fig.delaxes(ax)
        self.content_axes.clear()

    def _switch(self, view: str):
        self.view = view
        self._clear_content()
        if view == "overview":
            self._draw_overview()
        else:
            self._draw_detail()
        self.fig.canvas.draw_idle()

    def _on_radio(self, label: str):
        self.sel = self.names.index(label)
        if self.view == "detail":
            self._switch("detail")
        else:
            # just highlight
            self._switch("overview")

    # ── Overview view ─────────────────────────────────────────────────────────

    def _draw_overview(self):
        """6 comparative panels for all files."""
        n = len(self.data)
        xs = np.arange(n)
        labels = self.names

        # grid: 2 rows × 3 cols
        gs = gridspec.GridSpec(2, 3, figure=self.fig,
                               left=0.10, right=0.98,
                               top=0.93, bottom=0.07,
                               wspace=0.35, hspace=0.52)

        def make_ax(r, c, title):
            ax = self.fig.add_subplot(gs[r, c])
            ax.set_title(title, pad=6, color=ACCENT1, fontsize=9)
            ax.set_facecolor(PANEL)
            ax.grid(axis="y", alpha=0.3)
            ax.set_xticks(xs)
            ax.set_xticklabels(labels, rotation=30, ha="right")
            self.content_axes.append(ax)
            return ax

        # helper: highlight selected bar
        def hbar(ax, vals, cmap_or_color=ACCENT1, sel=self.sel):
            colors = [ACCENT2 if i == sel else ACCENT1 for i in range(n)]
            bars = ax.bar(xs, vals, color=colors, width=0.6, zorder=3,
                          linewidth=0, edgecolor="none")
            for bar, v in zip(bars, vals):
                ax.text(bar.get_x() + bar.get_width()/2,
                        bar.get_height() + max(vals)*0.01,
                        f"{v:g}", ha="center", va="bottom",
                        fontsize=6.5, color=TEXT_CLR)
            return bars

        # 1. Matrix size  (scatter: cols × rows, bubble = ones)
        ax1 = make_ax(0, 0, "Matrix Dimensions  (rows × cols)")
        for i, d in enumerate(self.data):
            clr = ACCENT2 if i == self.sel else ACCENT1
            sz  = d["total_ones"] / 5
            ax1.scatter(d["num_cols"], d["num_rows"], s=sz, color=clr,
                        alpha=0.85, zorder=4, edgecolors="white", linewidths=0.4)
            ax1.annotate(d["filename"],
                         (d["num_cols"], d["num_rows"]),
                         textcoords="offset points", xytext=(5, 3),
                         fontsize=7, color=clr)
        ax1.set_xlabel("num_cols")
        ax1.set_ylabel("num_rows")
        ax1.set_xticks([])
        note = mpatches.Patch(color=PANEL, label="bubble ∝ total_ones")
        ax1.legend(handles=[note], fontsize=6, loc="lower right",
                   facecolor=BG, edgecolor=DIM_CLR, labelcolor=TEXT_CLR)

        # 2. Density
        ax2 = make_ax(0, 1, "Density  (ones / cells)")
        vals = [density(d) * 100 for d in self.data]
        hbar(ax2, vals)
        ax2.set_ylabel("density %")

        # 3. Total ones
        ax3 = make_ax(0, 2, "Total Ones")
        hbar(ax3, [d["total_ones"] for d in self.data])
        ax3.set_ylabel("count")

        # 4. Singleton ratios  (stacked)
        ax4 = make_ax(1, 0, "Singleton Ratios")
        sr  = [singleton_row_ratio(d)*100 for d in self.data]
        sc  = [singleton_col_ratio(d)*100 for d in self.data]
        bar1 = ax4.bar(xs - 0.18, sr, 0.35, color=ACCENT3, label="rows %", zorder=3)
        bar2 = ax4.bar(xs + 0.18, sc, 0.35, color=ACCENT4, label="cols %", zorder=3)
        ax4.set_xticks(xs); ax4.set_xticklabels(labels, rotation=30, ha="right")
        ax4.set_ylabel("%")
        ax4.legend(fontsize=7, facecolor=BG, edgecolor=DIM_CLR, labelcolor=TEXT_CLR)
        # highlight selected
        bar1[self.sel].set_edgecolor(ACCENT2); bar1[self.sel].set_linewidth(2)
        bar2[self.sel].set_edgecolor(ACCENT2); bar2[self.sel].set_linewidth(2)

        # 5. Subset counts
        ax5 = make_ax(1, 1, "Subset Relationships")
        rs = [d["row_subsets"] for d in self.data]
        cs = [d["col_subsets"] for d in self.data]
        ax5.bar(xs - 0.18, rs, 0.35, color=ACCENT5, label="row subsets", zorder=3)
        ax5.bar(xs + 0.18, cs, 0.35, color=ACCENT2, label="col subsets", zorder=3)
        ax5.set_xticks(xs); ax5.set_xticklabels(labels, rotation=30, ha="right")
        ax5.set_ylabel("count")
        ax5.legend(fontsize=7, facecolor=BG, edgecolor=DIM_CLR, labelcolor=TEXT_CLR)
        ax5.yaxis.get_major_locator().set_params(integer=True)

        # 6. Empty rows / cols
        ax6 = make_ax(1, 2, "Empty Rows & Cols")
        er = [d["empty_rows"] for d in self.data]
        ec = [d["empty_cols"] for d in self.data]
        ax6.bar(xs - 0.18, er, 0.35, color=ACCENT3, label="empty rows", zorder=3)
        ax6.bar(xs + 0.18, ec, 0.35, color=ACCENT4, label="empty cols", zorder=3)
        ax6.set_xticks(xs); ax6.set_xticklabels(labels, rotation=30, ha="right")
        ax6.set_ylabel("count")
        ax6.legend(fontsize=7, facecolor=BG, edgecolor=DIM_CLR, labelcolor=TEXT_CLR)
        ax6.yaxis.get_major_locator().set_params(integer=True)

        # footer tip
        self.fig.text(0.54, 0.01,
                      "Click a bar to select  ·  Use radio buttons to pick file  ·  Switch to Detail for subset explorer",
                      ha="center", fontsize=7, color=DIM_CLR)

        # click on bars to select file
        def on_pick(event):
            if hasattr(event, "ind"):
                idx = event.ind[0]
                if 0 <= idx < n:
                    self.sel = idx
                    self.radio.set_active(idx)
        self.fig.canvas.mpl_connect("pick_event", on_pick)

    # ── Detail view ───────────────────────────────────────────────────────────

    def _draw_detail(self):
        d = self.data[self.sel]
        gs = gridspec.GridSpec(3, 3, figure=self.fig,
                               left=0.10, right=0.98,
                               top=0.93, bottom=0.07,
                               wspace=0.40, hspace=0.58)

        def ax_(r, c, rs=1, cs=1, title=""):
            ax = self.fig.add_subplot(gs[r:r+rs, c:c+cs])
            ax.set_facecolor(PANEL)
            if title:
                ax.set_title(title, pad=6, color=ACCENT1, fontsize=9)
            self.content_axes.append(ax)
            return ax

        # ── Info card ────────────────────────────────────────────────────────
        ax_info = ax_(0, 0, title=f"📄  {d['filename']}")
        ax_info.axis("off")
        stats = [
            ("rows",         d["num_rows"]),
            ("cols",         d["num_cols"]),
            ("total ones",   d["total_ones"]),
            ("density",      f"{density(d)*100:.3f}%"),
            ("empty rows",   d["empty_rows"]),
            ("empty cols",   d["empty_cols"]),
            ("singleton rows", d["singleton_rows"]),
            ("singleton cols", d["singleton_cols"]),
            ("row subsets",  d["row_subsets"]),
            ("col subsets",  d["col_subsets"]),
        ]
        for i, (k, v) in enumerate(stats):
            ypos = 0.95 - i * 0.095
            ax_info.text(0.02, ypos, f"  {k}", transform=ax_info.transAxes,
                         fontsize=8, color=DIM_CLR, va="top")
            ax_info.text(0.98, ypos, str(v), transform=ax_info.transAxes,
                         fontsize=8, color=ACCENT1, va="top", ha="right",
                         fontweight="bold")
            if i < len(stats)-1:
                ax_info.axhline(y=ypos * ax_info.get_ylim()[1] if ax_info.get_ylim()[1] != 1 else ypos,
                                color=GRID_CLR, linewidth=0.4,
                                transform=ax_info.transAxes)

        # ── Density heatmap (simulated row/col profile) ──────────────────────
        ax_heat = ax_(0, 1, cs=2, title="Row / Col Ones Profile  (simulated profile)")
        nr, nc = d["num_rows"], d["num_cols"]
        # approximate distribution: place ones randomly consistent with metadata
        np.random.seed(42)
        ones = d["total_ones"]
        mat  = np.zeros((nr, nc), dtype=np.uint8)
        idx  = np.random.choice(nr*nc, min(ones, nr*nc), replace=False)
        mat.flat[idx] = 1

        row_sums = mat.sum(axis=1)
        col_sums = mat.sum(axis=0)

        # show row & col profiles side by side as heatmap strips
        combined = np.vstack([
            row_sums[np.newaxis, :] / (row_sums.max() or 1),
            col_sums[np.newaxis, :nr] if nc >= nr else
               np.pad(col_sums, (0, nr-nc), constant_values=0)[np.newaxis, :nr]
        ])
        ax_heat.imshow(combined, aspect="auto", cmap=CMAP_DENSITY, vmin=0, vmax=1,
                       interpolation="nearest")
        ax_heat.set_yticks([0, 1])
        ax_heat.set_yticklabels(["rows →", "cols →"], fontsize=7)
        ax_heat.set_xlabel("index", fontsize=7)
        cbar_ax = self.fig.add_axes([
            ax_heat.get_position().x1 + 0.005,
            ax_heat.get_position().y0,
            0.008,
            ax_heat.get_position().height
        ])
        self.content_axes.append(cbar_ax)
        sm = plt.cm.ScalarMappable(cmap=CMAP_DENSITY)
        sm.set_array([])
        self.fig.colorbar(sm, cax=cbar_ax, orientation="vertical")
        cbar_ax.tick_params(labelsize=6, colors=DIM_CLR)

        # ── Row singletons + subsets bar ─────────────────────────────────────
        ax_rb = ax_(1, 0, title="Row Stats")
        categories = ["total", "singleton", "subset"]
        vals = [d["num_rows"], d["singleton_rows"], d["row_subsets"]]
        colors = [ACCENT1, ACCENT3, ACCENT5]
        bars = ax_rb.barh(categories, vals, color=colors, height=0.5, zorder=3)
        ax_rb.grid(axis="x", alpha=0.3)
        for bar, v in zip(bars, vals):
            ax_rb.text(v + max(vals)*0.01, bar.get_y()+bar.get_height()/2,
                       str(v), va="center", fontsize=8, color=TEXT_CLR)
        ax_rb.set_xlim(0, max(vals)*1.18)
        ax_rb.invert_yaxis()

        # ── Col singletons + subsets bar ─────────────────────────────────────
        ax_cb = ax_(1, 1, title="Col Stats")
        vals_c = [d["num_cols"], d["singleton_cols"], d["col_subsets"]]
        colors_c = [ACCENT1, ACCENT4, ACCENT2]
        bars_c = ax_cb.barh(categories, vals_c, color=colors_c, height=0.5, zorder=3)
        ax_cb.grid(axis="x", alpha=0.3)
        for bar, v in zip(bars_c, vals_c):
            ax_cb.text(v + max(vals_c)*0.01, bar.get_y()+bar.get_height()/2,
                       str(v), va="center", fontsize=8, color=TEXT_CLR)
        ax_cb.set_xlim(0, max(vals_c)*1.18)
        ax_cb.invert_yaxis()

        # ── Subset explorer ──────────────────────────────────────────────────
        ax_sub = ax_(1, 2, rs=2, title="Subset Relationships")
        ax_sub.axis("off")
        row_sds = d.get("row_subset_details", [])
        col_sds = d.get("column_subset_details", [])
        all_sds = [("R", sd) for sd in row_sds] + [("C", sd) for sd in col_sds]

        if not all_sds:
            ax_sub.text(0.5, 0.5, "No subset relationships",
                        ha="center", va="center", fontsize=9, color=DIM_CLR,
                        transform=ax_sub.transAxes)
        else:
            y = 0.97
            for kind, sd in all_sds:
                if kind == "R":
                    sub_k, sup_k = "subset_row", "superset_row"
                    tag, clr = "ROW", ACCENT5
                else:
                    sub_k, sup_k = "subset_column", "superset_column"
                    tag, clr = "COL", ACCENT2

                sub_id = sd[sub_k]
                sup_id = sd[sup_k]
                sub_ones = sd["subset_ones"]
                sup_ones = sd["superset_ones"]

                ax_sub.text(0.02, y, f"[{tag}]", transform=ax_sub.transAxes,
                            fontsize=7.5, color=clr, va="top", fontweight="bold")
                ax_sub.text(0.14, y,
                            f"#{sub_id}  ⊂  #{sup_id}",
                            transform=ax_sub.transAxes,
                            fontsize=7.5, color=TEXT_CLR, va="top")
                y -= 0.055
                ax_sub.text(0.06, y,
                            f"subset:   {sub_ones}",
                            transform=ax_sub.transAxes,
                            fontsize=6.5, color=DIM_CLR, va="top")
                y -= 0.047
                ax_sub.text(0.06, y,
                            f"superset: {sup_ones}",
                            transform=ax_sub.transAxes,
                            fontsize=6.5, color=DIM_CLR, va="top")
                y -= 0.06
                if y < 0.05:
                    ax_sub.text(0.5, y, f"… +{len(all_sds)-all_sds.index((kind,sd))-1} more",
                                transform=ax_sub.transAxes, ha="center",
                                fontsize=7, color=DIM_CLR)
                    break

        # ── Col subset Venn-style network ─────────────────────────────────────
        ax_net = ax_(2, 0, cs=2, title="Subset Network  (cols)")
        ax_net.set_facecolor(PANEL)
        ax_net.set_aspect("equal")
        col_sds2 = d.get("column_subset_details", [])

        if not col_sds2:
            ax_net.text(0.5, 0.5, "No column subset relationships",
                        ha="center", va="center", fontsize=9, color=DIM_CLR,
                        transform=ax_net.transAxes)
            ax_net.axis("off")
        else:
            # collect unique nodes
            nodes = {}
            edges = []
            for sd in col_sds2:
                sub, sup = sd["subset_column"], sd["superset_column"]
                if sub not in nodes: nodes[sub] = len(nodes)
                if sup not in nodes: nodes[sup] = len(nodes)
                edges.append((sub, sup, len(sd["subset_ones"]), len(sd["superset_ones"])))

            nn = len(nodes)
            angles = np.linspace(0, 2*np.pi, nn, endpoint=False)
            pos = {n: (np.cos(a), np.sin(a)) for n, a in zip(nodes, angles)}

            # edges
            for sub, sup, sz_sub, sz_sup in edges:
                x0, y0 = pos[sub]; x1, y1 = pos[sup]
                ax_net.annotate("",
                    xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=ACCENT2,
                                    lw=1.4, mutation_scale=12,
                                    connectionstyle="arc3,rad=0.15"))

            # nodes
            for col_id, idx in nodes.items():
                x, y = pos[col_id]
                sizes = [sd["subset_ones"] for sd in col_sds2 if sd["subset_column"] == col_id]
                size  = len(sizes[0]) * 80 + 120 if sizes else 200
                ax_net.scatter(x, y, s=size, color=ACCENT1, zorder=5,
                               edgecolors=BG, linewidths=1.5)
                ax_net.text(x, y + 0.12, f"c{col_id}",
                            ha="center", va="bottom", fontsize=7.5,
                            color=TEXT_CLR, fontweight="bold")

            ax_net.set_xlim(-1.5, 1.5); ax_net.set_ylim(-1.5, 1.5)
            ax_net.axis("off")

        self.fig.text(0.54, 0.01,
                      "Switch files with the radio selector · Click 'Overview' for cross-file comparison",
                      ha="center", fontsize=7, color=DIM_CLR)


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = os.path.join(os.path.dirname(__file__), "matrix_data.json")

    if not os.path.isfile(path):
        print(f"[ERROR] File not found: {path}")
        sys.exit(1)

    print(f"Loading data from: {path}")
    data = load_data(path)
    print(f"Loaded {len(data)} matrices: {[d['filename'] for d in data]}")

    MatrixViz(data)


if __name__ == "__main__":
    main()
