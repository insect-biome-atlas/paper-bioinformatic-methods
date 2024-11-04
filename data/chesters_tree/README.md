# Modifications to make data compatible with EPA-NG and GAPPA

1. Trim gaps in the alignment to only keep positions with <90% gap characters

```bash
python trim_gaps.py -i chesters_new_outgroups_aligned.fasta.gz -o chesters_new_outgroups_aligned.trim0.9.fasta -f 0.9
```

2. Reformat the taxonomy file. 

Here, missing values are filled with '', columns are filtered and renamed and a
taxonomy file compatible with EPA-NG/GAPPA is written.


> [!NOTE] 
> With commit
> [85775a6](https://github.com/insect-biome-atlas/paper-bioinformatic-methods/commit/85775a6a67398b642e9dff23086f5cd6878d6755)
> the taxonomy file was updated with missing species names as well as missing
> genus for _Anoploderomorpha izumii_.

```python
import pandas as pd

# Read in the taxonomy file
df = pd.read_csv("chesters_new_outgroups_taxonomy.tsv", sep="\t", header=0, index_col=0)
df.rename(columns=lambda x: x.lower(), inplace=True)
# Fill empyty
df.fillna("", inplace=True)
# Drop the columns we don't need
df = df.loc[:, ["kingdom", "phylum", "class", "order", "family", "genus", "species"]]
df.to_csv("chesters_new_outgroups_taxonomy.updated.tsv", sep="\t")
df["lineage"] = df["kingdom"]+";"+df["phylum"]+";"+df["class"]+";"+df["order"]+";"+df["family"]+";"+df["genus"]+";"+df["species"]
df.loc[:, "lineage"].to_csv("taxonomy.tsv", sep="\t", header=False)
```