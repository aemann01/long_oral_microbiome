### 1. DADA2 processing

```bash
# to run jupyter notebook with all DADA2 processing steps:
jupyter-notebook dada2_processing.ipynb
```

### 2. Taxonomic assignment with Kraken2

Once quality filter/ASV generation/etc steps are complete, can now assign taxonomy using Kraken2 (database built from Mann et al. ////)

[Code to build the rpoC database](https://github.com/aemann01/domhain/tree/main/2022-HIV_oral_microbiome/00-database_build)

```bash
kraken2 --db ~/refdb/kraken_rpoc \
	--threads 40 \
	--use-names \
	--output rep_set.kraken.out rep_set.fa \
	--unclassified-out rep_set.unclassified.out --confidence 0.01
```

Most sequences should be successfully classified (in this run, 96% were classified, ~4% unclassified)

### 3. Get taxonomic identifier for each accession number, pull ranked lineages

## TO DO: links to download appropriate ncbi databases

```bash
awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > taxids
sed "s/\t//g" ~/refdb/ncbi_taxonomy/rankedlineage.dmp > rankedlineage.dmp
sort -t "|" -k 1b,1 rankedlineage.dmp > rankedlineage_sorted
cat rankedlineage_sorted | sed 's/|\{2,\}/|/g' > rankedlineage_clean
# add unclassified to taxonomy file
sed -i '1 i\0|unclassified|' rankedlineage_clean
# check if taxids exist in ranked lineage file
cat taxids | sed 's/$/|/' | while read line; do grep -c -m 1 ^$line rankedlineage_clean | sed "s/0/$line not found/"; done > missing_check
grep "not found" missing_check # should come back with nothing if all taxids found
cat taxids | sed 's/$/|/' | while read line; do grep -m 1 ^$line rankedlineage_clean || echo $line "no lineage" ; done > lineages
grep "no lineage" lineages # should return empty
# remove taxids from file and reverse field order
sed 's/|/\t/' lineages | awk -F"\t" '{print $2}' | awk -F\| '{s=$NF;for(i=NF-1;i>=1;i--)s=s FS $i;print s}' | sed 's/^|//' | sed 's/ /_/g' | sed 's/|/;/g' > taxonomy
# merge asv ids and taxonomy
awk '{print $2}' rep_set.kraken.out > asvids
paste asvids taxonomy > taxonomy.txt
```

### 4. Filter taxonomy and fix tax strings

Here are only keeping those ASVs that were classified as bacteria at the phylum level or below

```bash
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt
```
