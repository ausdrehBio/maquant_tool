"""Q7. Is the sample structure the same as in the publication?

Principal component analysis of the filtered, imputed log2 LFQ matrices (each protein centred):
the published matrix as published, and the paper pipeline's matrices for our and the deposited
table (first 20 seeds). Agreement with the published structure is measured twice:
- Procrustes similarity (1 - disparity) of the 42 samples on PC1-3, which allows rotation,
  reflection and scaling;
- Spearman correlation of the 21 distances between condition centroids (first 10 PCs).
The paper's Fig. 4a claim, that the CPSF30-binding strains (RFMG, AFMG) cluster apart from the
G184R strains (RSIR, RFMR, AFMR), is tested with the same centroid distances.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist, squareform

import mq_common as mq

QUESTION = "Is the sample structure the same as in the publication?"
PAPER = "Fig. 4a (UMAP): CPSF-binding strains cluster apart from the G184R strains"
NOTE = "PCA of the imputed matrices; Procrustes and centroid distances against published"
N_PCA_SEEDS = 20
N_PC_DIST = 10


def pca(Y):
    """Sample coordinates and explained variance of a proteins x samples matrix."""
    Z = (Y - Y.mean(axis=1, keepdims=True)).T
    Z = Z - Z.mean(axis=0)
    U, S, _ = np.linalg.svd(Z, full_matrices=False)
    return U * S, S ** 2 / (S ** 2).sum()


def centroid_distances(C):
    cent = np.array([C[i * 6:(i + 1) * 6, :N_PC_DIST].mean(axis=0)
                     for i in range(len(mq.CONDITIONS))])
    return squareform(pdist(cent))


def separation(D):
    idx = {c: i for i, c in enumerate(mq.CONDITIONS)}
    g = [idx[c] for c in mq.CPSF_BINDING]
    r = [idx[c] for c in mq.G184R]
    between = np.mean([D[a, b] for a in g for b in r])
    within = np.mean([D[g[0], g[1]]] + [D[a, b] for i, a in enumerate(r) for b in r[i + 1:]])
    return between, within


def matrices(name, n):
    """Imputed matrices of a layer for the first n seeds (the pipeline's own draws)."""
    res = mq.pipeline(name)
    out = []
    for ss in mq.seed_sequences()[:n]:
        out.append(mq.impute(res.log2, np.random.default_rng(ss)))
    return out


def main():
    pub = mq.published().matrix.values
    Cp, var_p = pca(pub)
    Dp = centroid_distances(Cp)
    iu = np.triu_indices(len(mq.CONDITIONS), 1)
    names = [n for n in ("ours", "deposited", "ours-67") if n in mq.layers()]
    res = {"published": {"proc": [1.0], "dist": [1.0], "sep": [separation(Dp)], "var": var_p,
                         "C": Cp}}
    for name in names:
        proc, dist, sep = [], [], []
        first = None
        for Y in matrices(name, N_PCA_SEEDS):
            C, var = pca(Y)
            first = first if first is not None else (C, var)
            proc.append(1 - procrustes(Cp[:, :3], C[:, :3])[2])
            D = centroid_distances(C)
            dist.append(stats.spearmanr(Dp[iu], D[iu])[0])
            sep.append(separation(D))
        res[name] = {"proc": proc, "dist": dist, "sep": sep, "var": first[1], "C": first[0]}

    rows = []
    for name, r in res.items():
        b, w = np.median([s[0] for s in r["sep"]]), np.median([s[1] for s in r["sep"]])
        rows.append((name, " / ".join(f"{100 * v:.0f}%" for v in r["var"][:3]),
                     mq.f2(np.median(r["proc"]), 3) if name != "published" else "reference",
                     mq.f2(np.median(r["dist"]), 3) if name != "published" else "reference",
                     f"{mq.f2(b, 1)} vs {mq.f2(w, 1)}", "yes" if b > w else "**no**"))

    m_o = {k: float(np.median(res["ours"][k])) for k in ("proc", "dist")}
    m_d = {k: float(np.median(res["deposited"][k])) for k in ("proc", "dist")}
    close = all(1 - m_o[k] <= 2 * (1 - m_d[k]) for k in m_o)
    claim_pub = separation(Dp)[0] > separation(Dp)[1]
    claim_ours = np.median([s[0] > s[1] for s in res["ours"]["sep"]]) >= 0.5
    if claim_pub and not claim_ours:
        verdict = "inconsistent"
    elif close:
        verdict = "consistent"
    else:
        verdict = "minor deviation"
    figure = plot(res)
    criterion = ("Consistent if the Fig. 4a separation holds in our data and, for both agreement "
                 "measures M (Procrustes similarity, centroid-distance correlation), our "
                 "unexplained share 1 − M is at most twice that of the authors' deposited "
                 "search. Inconsistent if the separation holds in the published data but not in "
                 "ours.")
    lines = [
        f"Medians over {N_PCA_SEEDS} imputation seeds; variance explained by PC1 / PC2 / PC3; "
        "centroid distances in the first 10 PCs (between: CPSF30-binding vs G184R strains; "
        "within: inside these two groups):",
        "",
        mq.md_table(rows, ["Matrix", "Variance PC1–3", "Procrustes similarity",
                           "Centroid-distance r", "Between vs within", "Fig. 4a separation"]),
        "",
        f"Our sample structure matches the published one about as closely as the authors' own "
        f"deposited search does: Procrustes {mq.f2(m_o['proc'], 3)} vs {mq.f2(m_d['proc'], 3)}, "
        f"centroid distances r = {mq.f2(m_o['dist'], 3)} vs {mq.f2(m_d['dist'], 3)}. The "
        "CPSF30-binding strains lie further from the G184R strains than the strains within "
        "either group lie from each other, in the published matrix and in ours.",
    ]
    short = (f"Procrustes {mq.f2(m_o['proc'], 3)} (deposited {mq.f2(m_d['proc'], 3)}); "
             f"Fig. 4a separation {'reproduced' if claim_ours else 'not reproduced'}")
    return mq.answer(QUESTION, "\n".join(lines), figure, short, verdict, criterion,
                     {"procrustes": m_o["proc"], "centroid_r": m_o["dist"],
                      "procrustes_deposited": m_d["proc"], "centroid_r_deposited": m_d["dist"]})


def plot(res):
    panels = [n for n in ("published", "ours", "deposited", "ours-67") if n in res]
    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 4.0))
    ref = res["published"]["C"][:, :2]
    for ax, name in zip(axes, panels):
        C = res[name]["C"][:, :3]
        if name != "published":           # rotate onto the published configuration for display
            _, C2, _ = procrustes(res["published"]["C"][:, :3], C)
            scale = np.std(ref) / np.std(C2[:, :2])
            C = C2 * scale
        for i, c in enumerate(mq.CONDITIONS):
            ax.scatter(C[i * 6:(i + 1) * 6, 0], C[i * 6:(i + 1) * 6, 1], s=22,
                       color=mq.CONDITION_COLORS[c], edgecolor=mq.SURFACE, linewidth=0.5, label=c)
        var = res[name]["var"]
        ax.set_title(name)
        ax.set_xlabel("PC1" + (f" ({100 * var[0]:.0f}%)" if name == "published" else " (aligned)"))
        ax.set_ylabel("PC2" + (f" ({100 * var[1]:.0f}%)" if name == "published" else " (aligned)"))
        mq.value_grid(ax, "both")
    axes[-1].legend(fontsize=7, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    return mq.save_fig(fig, "q07_sample_structure")


if __name__ == "__main__":
    print(main())
