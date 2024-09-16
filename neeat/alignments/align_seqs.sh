#! /bin/bash
# Construct aligned data matrix in FASTA format for
# computing evolutionary distances

# Define order names
orders=()
while IFS= read -r line; do
  orders+=("$line")
done < orders.tsv

# Align aa sequences for each order using mafft --auto
# Convert to nucleotide alignments using pal2nal
for ord in "${orders[@]}"; do
    echo "Now processing $ord"
    mafft --auto ${ord}_aa.fasta > ${ord}_aa_aligned.fasta
    pal2nal.pl ${ord}_aa_aligned.fasta ${ord}.fasta -codontable 5 -output fasta > ${ord}_aligned.fasta
done

