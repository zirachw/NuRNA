import math

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

from .constants import (
    NT_COLORS, NT_UNKNOWN,
    BACKBONE_COLOR, BACKBONE_COLOR_CMP, PAIR_COLOR, ARC_COLOR,
    PAIR_COLORS,
)

# Sequential colormap: warm cream to slate-blue (matches app background -> accent)
_HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "nurna_seq", ["#F7F3EE", "#4A6FA5"]
)


def _rna_layout(dot_bracket: str) -> list[tuple[float, float]]:
    """Recursive stem-loop layout. Returns (x, y) for each nucleotide."""
    n = len(dot_bracket)
    pair = [-1] * n
    stack: list[int] = []
    for k, c in enumerate(dot_bracket):
        if c == "(":
            stack.append(k)
        elif c == ")" and stack:
            j = stack.pop()
            pair[j] = k
            pair[k] = j

    xy: list[tuple[float, float] | None] = [None] * n
    STEP = 2.2
    SEP  = 3.0

    def loop_elems(lo: int, hi: int):
        result, k = [], lo
        while k <= hi:
            if 0 <= pair[k] > k and pair[k] <= hi:
                result.append((True, k, pair[k]))
                k = pair[k] + 1
            else:
                result.append((False, k, k))
                k += 1
        return result

    def place(i: int, j: int, x0: float, y0: float, angle: float) -> None:
        stem, ci, cj = [], i, j
        while ci < cj and pair[ci] == cj:
            stem.append((ci, cj))
            ci += 1
            cj -= 1

        dx, dy = math.cos(angle), math.sin(angle)
        px, py = -dy, dx
        for t, (si, sj) in enumerate(stem):
            cx = x0 + t * STEP * dx
            cy = y0 + t * STEP * dy
            xy[si] = (cx + SEP / 2 * px, cy + SEP / 2 * py)
            xy[sj] = (cx - SEP / 2 * px, cy - SEP / 2 * py)

        tip_x = x0 + len(stem) * STEP * dx
        tip_y = y0 + len(stem) * STEP * dy

        if ci > cj:
            return

        elems = loop_elems(ci, cj)
        inner_stems = [(b, c) for (a, b, c) in elems if a]
        n_inner = len(inner_stems)

        if n_inner == 0:
            # Hairpin loop: arc of unpaired nucleotides
            n_free = len(elems)
            r = max(0.8, n_free * 0.45 / math.pi)
            for t, (_, idx, _) in enumerate(elems):
                a = angle + math.pi * (t + 1) / (n_free + 1)
                xy[idx] = (tip_x + r * math.cos(a), tip_y + r * math.sin(a))

        elif n_inner == 1:
            # Interior loop or bulge: continue same direction
            for (is_s, idx, _) in elems:
                if not is_s and xy[idx] is None:
                    xy[idx] = (tip_x, tip_y)
            place(inner_stems[0][0], inner_stems[0][1], tip_x, tip_y, angle)

        else:
            # Multi-loop: fan branches symmetrically around the forward direction
            spread = min(math.pi * 0.9, n_inner * math.pi / 3)
            r = max(2.2, len(elems) * 0.35)
            for t, (si, sj) in enumerate(inner_stems):
                frac = t / (n_inner - 1) if n_inner > 1 else 0.5
                a = angle - spread / 2 + spread * frac
                place(si, sj, tip_x + r * math.cos(a), tip_y + r * math.sin(a), a)
            for (is_s, idx, _) in elems:
                if not is_s and xy[idx] is None:
                    xy[idx] = (tip_x, tip_y)

    outer = loop_elems(0, n - 1)
    outer_stems = [(b, c) for (a, b, c) in outer if a]
    n_outer = len(outer_stems)

    if n_outer == 0:
        for k in range(n):
            xy[k] = (k * STEP, 0.0)
    elif n_outer == 1:
        place(outer_stems[0][0], outer_stems[0][1], 0.0, 0.0, math.pi / 2)
        for (is_s, idx, _) in outer:
            if not is_s and xy[idx] is None:
                xy[idx] = (0.0, 0.0)
    else:
        r0 = max(6.0, len(outer) * 0.8)
        for t, (is_s, b, c) in enumerate(outer):
            a = 2 * math.pi * t / len(outer) - math.pi / 2
            x0, y0 = r0 * math.cos(a), r0 * math.sin(a)
            if is_s:
                place(b, c, x0, y0, a)
            else:
                xy[b] = (x0, y0)

    return [(v[0], v[1]) if v else (0.0, 0.0) for v in xy]


def _draw_nts(ax, seq: str, xs, ys, *, neutral: bool = False) -> None:
    for i, (x, y) in enumerate(zip(xs, ys)):
        color = NT_UNKNOWN if neutral else NT_COLORS.get(seq[i], NT_UNKNOWN)
        ax.scatter(x, y, s=60, color=color, zorder=3,
                   edgecolors="#444" if neutral else "none", linewidths=0.4)


def _nt_legend_handles():
    return [mpatches.Patch(color=NT_COLORS[nt], label=nt) for nt in "AUGC"]


def _pairs_to_db(pairs: set[tuple[int, int]], n: int) -> str:
    db = ["."] * n
    for i, j in pairs:
        db[i] = "("
        db[j] = ")"
    return "".join(db)


def plot_rna_2d(
    seq: str,
    pairs: set[tuple[int, int]],
    title: str = "",
    ax: plt.Axes | None = None,
) -> plt.Figure:
    n = len(seq)
    xy = _rna_layout(_pairs_to_db(pairs, n))
    xs = [p[0] for p in xy]
    ys = [p[1] for p in xy]

    own_fig = ax is None
    if own_fig:
        PAD = 1.5
        x_min, x_max = min(xs) - PAD, max(xs) + PAD
        y_min, y_max = min(ys) - PAD, max(ys) + PAD
        span = max(x_max - x_min, y_max - y_min)
        cx = (x_min + x_max) / 2
        cy = (y_min + y_max) / 2
        fig, ax = plt.subplots(figsize=(7, 7))
    else:
        fig = ax.figure

    ax.set_aspect("equal")
    ax.axis("off")
    if own_fig:
        ax.set_xlim(cx - span / 2, cx + span / 2)
        ax.set_ylim(cy - span / 2, cy + span / 2)

    for k in range(n - 1):
        ax.plot([xs[k], xs[k + 1]], [ys[k], ys[k + 1]],
                color=BACKBONE_COLOR, linewidth=1.2, zorder=1)
    for i, j in pairs:
        ax.plot([xs[i], xs[j]], [ys[i], ys[j]],
                color=PAIR_COLOR, linewidth=1.0, zorder=2)

    _draw_nts(ax, seq, xs, ys)
    ax.legend(handles=_nt_legend_handles(), loc="lower right",
              fontsize=14, framealpha=0.85, title="Nukleotida",
              title_fontsize=15, handlelength=2.5, handleheight=1.8,
              borderpad=1.0, labelspacing=0.8)

    if title:
        ax.set_title(title, fontsize=20, pad=12)
    if own_fig:
        fig.tight_layout()
    return fig


def plot_arc_diagram(seq: str, pairs: set[tuple[int, int]], title: str = "") -> plt.Figure:
    n = len(seq)
    fig, ax = plt.subplots(figsize=(max(8, n * 0.18), 4))
    ax.axis("off")

    ax.plot([0, n - 1], [0, 0], color=BACKBONE_COLOR, linewidth=1.5, zorder=1)

    for i, nt in enumerate(seq):
        ax.scatter(i, 0, s=80, color=NT_COLORS.get(nt, NT_UNKNOWN), zorder=3)

    for i, j in pairs:
        ax.add_patch(mpatches.Arc(
            ((i + j) / 2, 0), width=j - i, height=(j - i) / 2,
            angle=0, theta1=0, theta2=180,
            color=ARC_COLOR, linewidth=1.0, zorder=2,
        ))

    ax.set_xlim(-1, n)
    ax.set_ylim(-0.5, max((j - i) / 2 for i, j in pairs) * 1.1 if pairs else 1)
    ax.legend(handles=_nt_legend_handles(), loc="upper right",
              fontsize=8, framealpha=0.85, title="Nukleotida")
    if title:
        ax.set_title(title, fontsize=11, pad=8)
    fig.tight_layout()
    return fig


def plot_dp_matrix(matrix: list[list[int]], seq: str, vmax: int | None = None) -> plt.Figure:
    n = len(seq)
    arr = np.array(matrix, dtype=float)
    _vmax = vmax if vmax is not None else max(1, int(arr.max()))
    fig, ax = plt.subplots(figsize=(max(6, n * 0.22), max(6, n * 0.22)))
    im = ax.imshow(
        arr, cmap=_HEATMAP_CMAP, aspect="auto",
        vmin=0, vmax=_vmax,
    )
    plt.colorbar(im, ax=ax, shrink=0.7, label="gamma[i][j]")

    step = max(1, n // 20)
    ticks = list(range(0, n, step))
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels([seq[t] for t in ticks], fontsize=7, fontfamily="monospace")
    ax.set_yticklabels([seq[t] for t in ticks], fontsize=7, fontfamily="monospace")
    ax.set_xlabel("j")
    ax.set_ylabel("i")
    ax.set_title("Matriks DP Nussinov (gamma)", fontsize=11)
    fig.tight_layout()
    return fig


def plot_comparison(
    seq: str,
    predicted: set[tuple[int, int]],
    reference: set[tuple[int, int]],
) -> plt.Figure:
    n = len(seq)
    xy = _rna_layout(_pairs_to_db(reference, n))
    xs = [p[0] for p in xy]
    ys = [p[1] for p in xy]

    tp = predicted & reference
    fp = predicted - reference
    fn = reference - predicted

    PAD = 1.5
    x_min, x_max = min(xs) - PAD, max(xs) + PAD
    y_min, y_max = min(ys) - PAD, max(ys) + PAD
    span = max(x_max - x_min, y_max - y_min)
    cx = (x_min + x_max) / 2
    cy = (y_min + y_max) / 2

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(cx - span / 2, cx + span / 2)
    ax.set_ylim(cy - span / 2, cy + span / 2)

    for k in range(n - 1):
        ax.plot([xs[k], xs[k + 1]], [ys[k], ys[k + 1]],
                color=BACKBONE_COLOR_CMP, linewidth=1.2, zorder=1)

    for label, pair_set in [("fn", fn), ("fp", fp), ("tp", tp)]:
        color = PAIR_COLORS[label]
        for i, j in pair_set:
            ax.plot([xs[i], xs[j]], [ys[i], ys[j]], color=color, linewidth=2.0, zorder=2)

    _draw_nts(ax, seq, xs, ys, neutral=True)

    ax.legend(
        handles=[mpatches.Patch(color=PAIR_COLORS[k], label=f"{k.upper()} ({len(v)})")
                 for k, v in [("tp", tp), ("fp", fp), ("fn", fn)]],
        loc="lower right", fontsize=8, framealpha=0.85, title="Pasangan Basa",
    )
    ax.set_title("Perbandingan Prediksi vs Referensi", fontsize=11, pad=10)
    fig.tight_layout()
    return fig
