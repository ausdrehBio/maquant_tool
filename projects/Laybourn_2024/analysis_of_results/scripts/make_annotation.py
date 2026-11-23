"""Write annotation/go0051607_defense_response_to_virus.tsv: the pig proteins (UniProt,
organism 9823) annotated with GO:0051607 "defense response to virus" or a child term.

The paper reports the share of DEPs linked to GO:0051607 per comparison (section 3.1: 9% in the
lower trachea after either infection, 4% and 5% in the upper trachea). q06 recomputes these
shares, for Table S3 and for our DEPs, with this set. It is a snapshot of the current UniProt
annotation, not the authors' version; the file records the date of the query and is kept in
git, so the analysis does not depend on the network. Run once:
    python3 <path to>/scripts/make_annotation.py
"""
import datetime
import io
import urllib.parse
import urllib.request

import pandas as pd

import lb_common as lb

QUERY = "(organism_id:9823) AND (go:0051607)"
URL = ("https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields=accession,gene_names,"
       "protein_name&query=" + urllib.parse.quote(QUERY))


def main():
    with urllib.request.urlopen(URL, timeout=120) as fh:
        df = pd.read_csv(io.StringIO(fh.read().decode()), sep="\t", dtype=str)
    df.columns = ["accession", "gene_names", "protein_name"]
    df = df.sort_values("accession")
    lb.ANNOTATION_DIR.mkdir(exist_ok=True)
    with open(lb.GO_TABLE, "w") as out:
        out.write(f"# UniProt query {QUERY!r} on {datetime.date.today().isoformat()}; "
                  f"{len(df)} entries\n")
        df.to_csv(out, sep="\t", index=False)
    print(f"written: {lb.rel(lb.GO_TABLE)} ({len(df)} entries)")


if __name__ == "__main__":
    main()
