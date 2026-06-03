"""
make_paper_figure.py — 论文核心图

图1：最终干净评测(EC中和+tie修复+产物不相交)下 4 模块 × {random, ec_only, tanimoto, LTR_chem} 的 r@5。
图2：去污染证据链(每阶段 artifact → 修正后 LTR 是否仍超 Tanimoto)。
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import kio  # noqa: E402

FIG = kio.safe_output_path("modular/ltr/figures/_m").parent


def parse(v):
    return float(str(v).split("±")[0])


def main():
    kio.setup_logging()
    m = pd.read_csv(kio.OUTPUTS / "modular" / "ltr" / "clean2_metrics.csv")
    methods = ["random", "ec_only", "tanimoto", "LTR_chem"]
    mods = ["A", "B", "C", "D"]
    # 图1：grouped bar r@5
    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(mods)); w = 0.2
    colors = {"random": "#bbb", "ec_only": "#f0ad4e", "tanimoto": "#5bc0de", "LTR_chem": "#27ae60"}
    for i, meth in enumerate(methods):
        vals = [parse(m[(m.module == mod) & (m.method == meth)]["r@5"].iloc[0])
                if len(m[(m.module == mod) & (m.method == meth)]) else 0 for mod in mods]
        ax.bar(x + (i - 1.5) * w, vals, w, label=meth, color=colors[meth])
    ax.set_xticks(x); ax.set_xticklabels([f"{mm}\n{mr}" for mm, mr in
                                          zip(mods, ["small-mol", "carb", "protein", "lipid"])])
    ax.set_ylabel("recall@5 (true product in top-5)"); ax.set_ylim(0, 1.05)
    ax.set_title("Learned structure ranker (LTR_chem) vs baselines\n(EC-neutralized + product-disjoint + scaffold split, random tie-break)")
    ax.legend(); ax.grid(axis="y", alpha=0.3)
    for i, mod in enumerate(mods):
        lt = parse(m[(m.module == mod) & (m.method == "LTR_chem")]["r@5"].iloc[0])
        ax.text(x[i] + 1.5 * w, lt + 0.02, f"{lt:.2f}", ha="center", fontsize=8, color="#27ae60")
    fig.tight_layout(); fig.savefig(FIG / "fig1_clean_recall5.png", dpi=150); plt.close(fig)

    # 图2：去污染证据链
    stages = ["v2 vocab-neg", "v3 same-sub\n(junk decoy)", "v4 +product-\ndisjoint", "v5 +silver\n(EC tag leak)", "FINAL clean\n(EC-neutral+tiefix)"]
    note = ["LTR≈Tanimoto", "LTR≫(junk artifact)", "LTR≫(but...)", "LTR≫(EC leak)", "LTR>Tan A/B/C\nec_only=random"]
    statuses = [0, 1, 1, 1, 2]  # 0=tie,1=confounded,2=clean-win
    fig2, ax2 = plt.subplots(figsize=(11, 3.2))
    cmap = {0: "#5bc0de", 1: "#e74c3c", 2: "#27ae60"}
    for i, (s, n, st) in enumerate(zip(stages, note, statuses)):
        ax2.barh(0, 1, left=i, color=cmap[st], edgecolor="white")
        ax2.text(i + 0.5, 0, f"{s}\n{n}", ha="center", va="center", fontsize=8, color="white")
    ax2.set_xlim(0, len(stages)); ax2.axis("off")
    ax2.set_title("De-confounding evidence chain: every artifact removed → honest LTR>similarity (clean)")
    fig2.tight_layout(); fig2.savefig(FIG / "fig2_deconfound_chain.png", dpi=150); plt.close(fig2)
    kio.log.info("saved figures to %s", FIG)
    print("figures:", [p.name for p in FIG.glob("*.png")])


if __name__ == "__main__":
    main()
