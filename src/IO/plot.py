import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.patches import Rectangle
import numpy as np


def plot_instance_analysis(instances):
    """
    Creates an interactive dashboard to visualize CBM instance statistics.

    Controls:
    - Next / Previous buttons to switch between instances.
    - The figure updates dynamically.

    Parameters
    ----------
    instances : list[dict]
        List of dictionaries containing instance analysis data.
    """

    class Dashboard:
        def __init__(self, instances):
            self.instances = instances
            self.idx = 0

            # Large figure with custom layout
            self.fig = plt.figure(figsize=(18, 10))
            self.fig.canvas.manager.set_window_title("CBM Instance Analysis Dashboard")

            gs = self.fig.add_gridspec(
                3,
                4,
                height_ratios=[1.1, 1.1, 1.3],
                width_ratios=[1.2, 1.2, 1.5, 1.5],
                hspace=0.45,
                wspace=0.35,
            )

            # Main axes
            self.ax_metrics = self.fig.add_subplot(gs[0, 0])
            self.ax_structure = self.fig.add_subplot(gs[0, 1])
            self.ax_singletons = self.fig.add_subplot(gs[0, 2:])
            self.ax_subset_summary = self.fig.add_subplot(gs[1, 0])
            self.ax_row_subsets = self.fig.add_subplot(gs[1, 1:3])
            self.ax_col_subsets = self.fig.add_subplot(gs[1, 3])
            self.ax_indices = self.fig.add_subplot(gs[2, :])

            # Buttons
            ax_prev = plt.axes([0.35, 0.01, 0.12, 0.05])
            ax_next = plt.axes([0.53, 0.01, 0.12, 0.05])

            self.btn_prev = Button(ax_prev, "⬅ Previous")
            self.btn_next = Button(ax_next, "Next ➡")

            self.btn_prev.on_clicked(self.prev_instance)
            self.btn_next.on_clicked(self.next_instance)

            self.update()

        def prev_instance(self, event):
            self.idx = (self.idx - 1) % len(self.instances)
            self.update()

        def next_instance(self, event):
            self.idx = (self.idx + 1) % len(self.instances)
            self.update()

        def clear_axes(self):
            for ax in [
                self.ax_metrics,
                self.ax_structure,
                self.ax_singletons,
                self.ax_subset_summary,
                self.ax_row_subsets,
                self.ax_col_subsets,
                self.ax_indices,
            ]:
                ax.clear()

        def update(self):
            self.clear_axes()

            data = self.instances[self.idx]

            filename = data["filename"]

            # =========================
            # Title
            # =========================
            self.fig.suptitle(
                f"CBM Instance Analysis — {filename} ({self.idx + 1}/{len(self.instances)})",
                fontsize=22,
                fontweight="bold",
                y=0.98,
            )

            # =========================
            # 1. Metrics Panel
            # =========================
            self.ax_metrics.axis("off")

            density = data["total_ones"] / (data["num_rows"] * data["num_cols"])

            metrics_text = (
                f"Rows: {data['num_rows']}\n"
                f"Columns: {data['num_cols']}\n"
                f"Total Ones: {data['total_ones']}\n"
                f"Density: {density:.2%}\n"
                f"Empty Rows: {data['empty_rows']}\n"
                f"Empty Cols: {data['empty_cols']}"
            )

            self.ax_metrics.text(
                0.05,
                0.95,
                metrics_text,
                va="top",
                fontsize=12,
                family="monospace",
                bbox=dict(boxstyle="round,pad=0.6", facecolor="#F7F7F7"),
            )
            self.ax_metrics.set_title("General Metrics", fontweight="bold")

            # =========================
            # 2. Matrix Structure
            # =========================
            self.ax_structure.axis("off")

            rows = data["num_rows"]
            cols = data["num_cols"]

            # Draw rectangle representing matrix shape
            rect = Rectangle((0, 0), cols, rows, facecolor="#DCEEFF", edgecolor="black")
            self.ax_structure.add_patch(rect)

            self.ax_structure.set_xlim(0, cols)
            self.ax_structure.set_ylim(0, rows)
            self.ax_structure.set_aspect("auto")
            self.ax_structure.invert_yaxis()

            self.ax_structure.text(
                cols / 2,
                rows / 2,
                f"{rows} × {cols}",
                ha="center",
                va="center",
                fontsize=14,
                fontweight="bold",
            )
            self.ax_structure.set_title("Matrix Shape", fontweight="bold")

            # =========================
            # 3. Singleton Comparison
            # =========================
            labels = ["Rows", "Columns"]
            singleton_values = [
                data["singleton_rows"],
                data["singleton_cols"],
            ]
            totals = [data["num_rows"], data["num_cols"]]
            percentages = [v / t * 100 for v, t in zip(singleton_values, totals)]

            bars = self.ax_singletons.bar(
                labels,
                percentages,
                width=0.55,
                alpha=0.85,
            )

            for bar, pct, raw in zip(bars, percentages, singleton_values):
                self.ax_singletons.text(
                    bar.get_x() + bar.get_width() / 2,
                    pct + 1,
                    f"{raw}\n({pct:.1f}%)",
                    ha="center",
                    fontsize=11,
                    fontweight="bold",
                )

            self.ax_singletons.set_ylim(0, max(percentages) * 1.35 + 1)
            self.ax_singletons.set_ylabel("% of total")
            self.ax_singletons.set_title("Singleton Rows vs Columns", fontweight="bold")
            self.ax_singletons.grid(axis="y", alpha=0.3)

            # =========================
            # 4. Subset Summary
            # =========================
            subset_values = [data["row_subsets"], data["col_subsets"]]
            subset_labels = ["Row Subsets", "Column Subsets"]

            bars = self.ax_subset_summary.bar(
                subset_labels,
                subset_values,
                width=0.55,
                alpha=0.85,
            )

            for bar, val in zip(bars, subset_values):
                self.ax_subset_summary.text(
                    bar.get_x() + bar.get_width() / 2,
                    val + 0.05,
                    str(val),
                    ha="center",
                    fontweight="bold",
                )

            self.ax_subset_summary.set_title("Subset Counts", fontweight="bold")
            self.ax_subset_summary.grid(axis="y", alpha=0.3)

            # =========================
            # 5. Row Subset Details
            # =========================
            row_details = data["row_subset_details"]
            if row_details:
                labels = [
                    f"{d['subset_row']}⊂{d['superset_row']}"
                    for d in row_details[:10]
                ]
                values = [len(d["subset_ones"]) for d in row_details[:10]]

                self.ax_row_subsets.barh(labels, values)
                self.ax_row_subsets.set_xlabel("Subset Size")
                self.ax_row_subsets.set_title(
                    "Row Subset Relationships", fontweight="bold"
                )
            else:
                self.ax_row_subsets.text(
                    0.5,
                    0.5,
                    "No row subsets",
                    ha="center",
                    va="center",
                    fontsize=12,
                    alpha=0.6,
                )
                self.ax_row_subsets.set_title(
                    "Row Subset Relationships", fontweight="bold"
                )

            # =========================
            # 6. Column Subset Details
            # =========================
            col_details = data["column_subset_details"]
            if col_details:
                labels = [
                    f"{d['subset_column']}⊂{d['superset_column']}"
                    for d in col_details[:10]
                ]
                values = [len(d["subset_ones"]) for d in col_details[:10]]

                self.ax_col_subsets.barh(labels, values)
                self.ax_col_subsets.set_xlabel("Size")
                self.ax_col_subsets.set_title(
                    "Column Subsets", fontweight="bold"
                )
            else:
                self.ax_col_subsets.text(
                    0.5,
                    0.5,
                    "No column subsets",
                    ha="center",
                    va="center",
                    fontsize=12,
                    alpha=0.6,
                )
                self.ax_col_subsets.set_title(
                    "Column Subsets", fontweight="bold"
                )

            # =========================
            # 7. Singleton Indices Distribution
            # =========================
            self.ax_indices.set_title(
                "Singleton Index Distribution", fontweight="bold"
            )

            row_idx = data["singleton_row_indices"]
            col_idx = data["singleton_col_indices"]

            if row_idx:
                self.ax_indices.scatter(
                    row_idx,
                    np.ones(len(row_idx)),
                    s=80,
                    label="Singleton Rows",
                    alpha=0.8,
                )

            if col_idx:
                self.ax_indices.scatter(
                    col_idx,
                    np.zeros(len(col_idx)),
                    s=50,
                    label="Singleton Columns",
                    alpha=0.8,
                )

            self.ax_indices.set_yticks([0, 1])
            self.ax_indices.set_yticklabels(["Columns", "Rows"])
            self.ax_indices.set_xlabel("Index")
            self.ax_indices.grid(alpha=0.3)
            self.ax_indices.legend()

            # =========================
            # Footer
            # =========================
            self.fig.text(
                0.5,
                0.065,
                "Use the buttons below to navigate between instances",
                ha="center",
                fontsize=10,
                alpha=0.7,
            )

            self.fig.canvas.draw_idle()

    Dashboard(instances)
    plt.show()