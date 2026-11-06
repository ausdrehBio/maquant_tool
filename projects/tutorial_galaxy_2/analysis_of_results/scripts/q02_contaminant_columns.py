"""Q2. In which columns (number) are the potential contaminants in the protein group
and evidence file respectively?

1-based position of the "Potential contaminant" column, as Galaxy counts columns.
The figure shows how many rows each file flags as potential contaminant.
"""
import mq_common as mq

QUESTION = ("In which columns (number) are the potential contaminants in the protein group "
            "and evidence file respectively?")
TUTORIAL = "Column 138 (protein groups) and 50 (evidence)"
NOTE = "positions depend on the MaxQuant version and the number of samples"


def main():
    pg, ev = mq.load_protein_groups(), mq.load_evidence()
    pos = {
        "proteinGroups.txt": list(pg.columns).index(mq.col(pg, "Potential contaminant")) + 1,
        "evidence.txt": list(ev.columns).index(mq.col(ev, "Potential contaminant")) + 1,
    }
    flagged = {"proteinGroups.txt": int(mq.is_contaminant(pg).sum()),
               "evidence.txt": int(mq.is_contaminant(ev).sum())}
    total = {"proteinGroups.txt": len(pg), "evidence.txt": len(ev)}

    fig, axes = mq.plt.subplots(1, 2, figsize=(8, 3.2))
    for ax, f in zip(axes, pos):
        b = ax.bar(["flagged", "not flagged"], [flagged[f], total[f] - flagged[f]],
                   color=[mq.NEUTRAL, mq.AXIS], width=0.6, **mq.bar_kw())
        mq.label_bars(ax, b, [mq.n(flagged[f]), mq.n(total[f] - flagged[f])])
        ax.set_title(f"{f}: column {pos[f]}")
        ax.set_ylabel("Rows")
        ax.tick_params(axis="x", length=0)
        mq.value_grid(ax, "y")
    fig.tight_layout()
    figure = mq.save_fig(fig, "q02_contaminant_columns")

    text = (
        f"**Column {pos['proteinGroups.txt']}** in proteinGroups.txt and **column "
        f"{pos['evidence.txt']}** in evidence.txt (\"Potential contaminant\", counted from 1). "
        f"They flag {mq.n(flagged['proteinGroups.txt'])} protein groups and "
        f"{mq.n(flagged['evidence.txt'])} evidence rows. Column positions change with the "
        "MaxQuant version and the number of samples, so scripts here select columns by name."
    )
    short = f"Column {pos['proteinGroups.txt']} (protein groups) and {pos['evidence.txt']} (evidence)"
    return mq.answer(QUESTION, text, figure, short=short)


if __name__ == "__main__":
    print(main())
