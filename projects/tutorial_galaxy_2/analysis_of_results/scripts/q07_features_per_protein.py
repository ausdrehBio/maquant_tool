"""Q7. In the complete experiment, how many features has a protein on average?

The tutorial takes the mean of column 8 of ProteinLevelData, NumMeasuredFeature: the
number of features measured for a protein in a run. The distinct features per protein
over all runs (FeatureLevelData) are reported alongside.
"""
import numpy as np

import mq_common as mq

QUESTION = "In the complete experiment, how many features has a protein on average?"
TUTORIAL = "Around 6 features per protein (mean in summary statistics)"
NOTE = "mean of NumMeasuredFeature, the tutorial's column 8"


def main():
    prot = mq.load_msstats("ProteinLevelData")
    feat = mq.load_msstats("FeatureLevelData")
    measured = prot[mq.col(prot, "NumMeasuredFeature")].astype(float)
    mean_measured = measured.mean()
    distinct = feat.groupby(mq.col(feat, "PROTEIN"))[mq.col(feat, "FEATURE")].nunique()

    fig, ax = mq.plt.subplots(figsize=(7, 3.6))
    upper = int(np.percentile(measured, 99))
    ax.hist(measured.clip(upper=upper), bins=range(0, upper + 2), color=mq.NEUTRAL,
            edgecolor=mq.SURFACE, linewidth=0.8)
    ax.axvline(mean_measured, color=mq.SERIES[0], linewidth=1.5)
    ax.text(mean_measured, ax.get_ylim()[1] * 0.95, f"  mean {mean_measured:.1f}",
            color=mq.INK, fontsize=9, va="top")
    ax.set_xlabel(f"Features measured per protein and run (last bin: ≥ {upper})")
    ax.set_ylabel("Protein-run combinations")
    ax.set_title("Measured features per protein")
    mq.value_grid(ax, "y")
    figure = mq.save_fig(fig, "q07_features_per_protein")

    text = (
        f"On average **{mean_measured:.1f} features** were measured per protein and run "
        f"(mean of NumMeasuredFeature, median {measured.median():.0f}). Over all runs together, "
        f"a protein has {distinct.mean():.1f} distinct features on average (median "
        f"{distinct.median():.0f})."
    )
    return mq.answer(QUESTION, text, figure,
                     short=f"{mean_measured:.1f} per protein and run ({distinct.mean():.1f} distinct)")


if __name__ == "__main__":
    print(main())
