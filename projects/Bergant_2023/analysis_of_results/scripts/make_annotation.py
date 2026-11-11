"""Write the annotation snapshot used by the suite (run once; the output is kept in git).

annotation/ensp_gene.tsv
    Gene symbol, HGNC ID and Swiss-Prot accessions for every Ensembl protein ID that occurs in
    our, the deposited or the published protein groups (GRCh37, Ensembl BioMart). The FASTA
    headers of the authors' database (>sp|ENSP...|p) carry no gene names.
annotation/go0051607_defense_response_to_virus.tsv
    Human Swiss-Prot proteins annotated to GO:0051607 "defense response to virus" or one of its
    descendants (is_a, part_of), from QuickGO. Used for the paper's Fig. 4f claim (q08).

Both files start with comment lines naming the source and the retrieval date. Rerunning the
script fetches the current annotation, which can differ from the committed snapshot.
Usage (from any directory):  python3 <path to>/scripts/make_annotation.py
"""
from datetime import date
import io
import urllib.parse
import urllib.request

import pandas as pd

import mq_common as mq

BIOMART = "https://grch37.ensembl.org/biomart/martservice"
BIOMART_QUERY = """<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">
<Dataset name="hsapiens_gene_ensembl" interface="default">
<Attribute name="ensembl_peptide_id"/><Attribute name="external_gene_name"/>
<Attribute name="hgnc_id"/><Attribute name="uniprotswissprot"/>
</Dataset></Query>"""
QUICKGO = ("https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?goId=GO:0051607"
           "&goUsage=descendants&goUsageRelationships=is_a,part_of&taxonId=9606&taxonUsage=exact"
           "&geneProductType=protein&geneProductSubset=Swiss-Prot"
           "&selectedFields=geneProductId&selectedFields=symbol&selectedFields=goId")
CLIENT = {"User-Agent": "maquant_tool-analysis/1.0"}


def fetch(url, data=None, accept=None):
    headers = dict(CLIENT)
    if accept:
        headers["Accept"] = accept
    req = urllib.request.Request(url, data=data, headers=headers)
    with urllib.request.urlopen(req, timeout=600) as r:
        return r.read().decode()


def biomart_release():
    """Display name of the GRCh37 Ensembl mart, e.g. 'Ensembl Genes 113'."""
    text = fetch(BIOMART + "?type=registry")
    for line in text.splitlines():
        if 'name="ENSEMBL_MART_ENSEMBL"' in line:
            for part in line.split():
                if part.startswith("displayName="):
                    return line.split('displayName="')[1].split('"')[0]
    return "unknown release"


def table_ids():
    """Normalised Ensembl protein IDs of all protein groups in the three tables."""
    ids = set()
    for layer in mq.layers(include_67=False).values():
        for s in layer.pg["Protein IDs"].fillna(""):
            ids |= {t for t in mq.id_set(s) if t.startswith("ENSP")}
    for s in mq.published().ids:
        ids |= {t for t in s if t.startswith("ENSP")}
    return ids


def main():
    today = date.today().isoformat()
    ids = table_ids()
    print(f"{len(ids)} Ensembl protein IDs in the three tables")

    data = urllib.parse.urlencode({"query": BIOMART_QUERY}).encode()
    mart = pd.read_csv(io.StringIO(fetch(BIOMART, data=data)), sep="\t", dtype=str).fillna("")
    mart.columns = ["ensp", "gene", "hgnc_id", "uniprot"]
    if len(mart) < 50000:
        raise SystemExit(f"BioMart returned only {len(mart)} rows, expected a full mapping")
    mart = mart[mart["ensp"].isin(ids)]
    genes = (mart.groupby("ensp")
             .agg(gene=("gene", lambda s: next((g for g in s if g), "")),
                  hgnc_id=("hgnc_id", lambda s: next((g for g in s if g), "")),
                  uniprot=("uniprot", lambda s: ";".join(sorted({u for u in s if u}))))
             .reindex(sorted(ids)).fillna("").reset_index())
    release = biomart_release()
    mq.ANNOTATION_DIR.mkdir(parents=True, exist_ok=True)
    with open(mq.GENE_TABLE, "w") as fh:
        fh.write(f"# Ensembl protein ID -> gene symbol, HGNC ID, Swiss-Prot accessions\n"
                 f"# source: Ensembl BioMart GRCh37 ({release}), hsapiens_gene_ensembl\n"
                 f"# retrieved: {today}; IDs: all ENSP IDs in ours, deposited and published\n")
        genes.to_csv(fh, sep="\t", index=False)
    found = (genes["gene"] != "").sum()
    print(f"written: {mq.rel(mq.GENE_TABLE)} ({len(genes)} IDs, {found} with a gene symbol, "
          f"{(genes['uniprot'] != '').sum()} with a Swiss-Prot accession)")

    go = pd.read_csv(io.StringIO(fetch(QUICKGO, accept="text/tsv")), sep="\t", dtype=str)
    go.columns = [c.strip().upper() for c in go.columns]
    go = (go.groupby("GENE PRODUCT ID")
          .agg(symbol=("SYMBOL", "first"), go_terms=("GO TERM", lambda s: ";".join(sorted(set(s)))))
          .reset_index().rename(columns={"GENE PRODUCT ID": "uniprot"}).sort_values("symbol"))
    with open(mq.GO_TABLE, "w") as fh:
        fh.write("# Human Swiss-Prot proteins annotated to GO:0051607 (defense response to virus)"
                 " or a descendant (is_a, part_of)\n"
                 "# source: QuickGO annotation downloadSearch, taxon 9606, all evidence codes\n"
                 f"# retrieved: {today}\n")
        go.to_csv(fh, sep="\t", index=False)
    print(f"written: {mq.rel(mq.GO_TABLE)} ({len(go)} proteins)")


if __name__ == "__main__":
    main()
