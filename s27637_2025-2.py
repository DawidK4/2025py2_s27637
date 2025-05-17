#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import pandas as pd, matplotlib.pyplot as plt

def fetch_records(taxid, min_len, max_len, email, api_key):
    Entrez.email, Entrez.api_key = email, api_key
    h = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]", usehistory="y")
    r = Entrez.read(h)
    webenv, qk, count = r["WebEnv"], r["QueryKey"], int(r["Count"])
    h = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", retmax=500, webenv=webenv, query_key=qk)
    records = []
    for rec in SeqIO.parse(h, "genbank"):
        l = len(rec.seq)
        if min_len <= l <= max_len:
            records.append((rec.id, l, rec.description))
    return records

def save_csv(records, filename):
    df = pd.DataFrame(records, columns=["Accession", "Length", "Description"])
    df.sort_values("Length", ascending=False).to_csv(filename, index=False)
    return df

def plot_lengths(df, filename):
    df = df.sort_values("Length", ascending=False)
    plt.figure(figsize=(10,5))
    plt.plot(df["Accession"], df["Length"], marker="o")
    plt.xticks(rotation=90)
    plt.xlabel("Accession")
    plt.ylabel("Length")
    plt.tight_layout()
    plt.savefig(filename)

if __name__ == "__main__":
    email = input("NCBI Email: ")
    api_key = input("NCBI API Key: ")
    taxid = input("TaxID: ")
    min_len = int(input("Min length: "))
    max_len = int(input("Max length: "))
    records = fetch_records(taxid, min_len, max_len, email, api_key)
    if not records: exit("No records found.")
    df = save_csv(records, f"taxid_{taxid}_report.csv")
    plot_lengths(df, f"taxid_{taxid}_plot.png")
    print("Done.")
