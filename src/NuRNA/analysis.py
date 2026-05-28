import pandas as pd
import matplotlib.pyplot as plt

from .constants import FAMILY_PALETTE


def analyze_by_family(results: list[dict]) -> pd.DataFrame:
    rows = []
    for r in results:
        m = r["metrics"]
        rows.append({
            "family":      r["family"],
            "sensitivity": m["sensitivity"],
            "ppv":         m["ppv"],
            "f1":          m["f1"],
        })
    df = pd.DataFrame(rows)
    return df.groupby("family")[["sensitivity", "ppv", "f1"]].mean().reset_index()


def plot_accuracy_vs_length(results: list[dict]) -> plt.Figure:
    lengths = [len(r["seq"]) for r in results]
    f1s     = [r["metrics"]["f1"] for r in results]
    families = [r["family"] for r in results]

    palette = _family_palette(families)

    fig, ax = plt.subplots(figsize=(8, 5))
    for fam in sorted(set(families)):
        xs = [l for l, f in zip(lengths, families) if f == fam]
        ys = [s for s, f in zip(f1s, families) if f == fam]
        ax.scatter(xs, ys, label=fam, color=palette[fam], s=20, alpha=0.7)

    ax.set_xlabel("Panjang sekuens (nt)")
    ax.set_ylabel("F1")
    ax.set_title("F1 vs Panjang Sekuens")
    ax.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    return fig


def plot_accuracy_vs_gc(results: list[dict]) -> plt.Figure:
    from .data import gc_content

    gcs      = [gc_content(r["seq"]) for r in results]
    f1s      = [r["metrics"]["f1"] for r in results]
    families = [r["family"] for r in results]

    palette = _family_palette(families)

    fig, ax = plt.subplots(figsize=(8, 5))
    for fam in sorted(set(families)):
        xs = [g for g, f in zip(gcs, families) if f == fam]
        ys = [s for s, f in zip(f1s, families) if f == fam]
        ax.scatter(xs, ys, label=fam, color=palette[fam], s=20, alpha=0.7)

    ax.set_xlabel("GC content")
    ax.set_ylabel("F1")
    ax.set_title("F1 vs GC Content")
    ax.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    return fig


def summarize(results: list[dict]) -> pd.DataFrame:
    from .data import gc_content

    rows = []
    for r in results:
        m = r["metrics"]
        rows.append({
            "id":          r.get("id", ""),
            "family":      r["family"],
            "length":      len(r["seq"]),
            "gc_content":  round(gc_content(r["seq"]), 3),
            "sensitivity": round(m["sensitivity"], 3),
            "ppv":         round(m["ppv"], 3),
            "f1":          round(m["f1"], 3),
            "tp":          m["tp"],
            "fp":          m["fp"],
            "fn":          m["fn"],
        })
    return pd.DataFrame(rows)


def _family_palette(families: list[str]) -> dict[str, str]:
    unique = sorted(set(families))
    return {fam: FAMILY_PALETTE[i % len(FAMILY_PALETTE)] for i, fam in enumerate(unique)}
