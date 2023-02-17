## Download bacterial assembly data from NCBI

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
cat assembly_summary.txt | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/$/\/*cds_from_genomic.fna.gz/' | sed 's/^https:/ftp:/' > query
wget -i query &
```

```bash
# make sure you got them all
grep "pattern" wget-log -B 14 | grep "2023" | grep "ftp:" | awk '{print $3}' > missing
# if genome has no cds, will not be downloaded, just a check to make sure that all possible are downloaded
wget -i missing 
# remove duplicates
rm *fna.gz.1
# unzip 
ls *fna.gz | parallel 'gzip -d {}'
# concatenate
cat *fna > all_genomes.fna
rm *genomic.fna
```

## Pull rpoC annotated sequences

```bash
grep "gene=rpoC" all_genomes.fna | grep -v "gamma" | awk '{print $1}' | sed 's/>//' > rpoc.ids
seqtk subseq all_genomes.fna rpoc.ids > rpoc.fna
```

## Importing sequences into QIIME2 to pull amplicons and dereplicate

```bash
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path rpoc.fna \
  --output-path rpoc.qza
qiime feature-classifier extract-reads \
    --i-sequences rpoc.qza \
    --p-f-primer MAYGARAARMGNATGYTNCARGA \
    --p-r-primer GMCATYTGRTCNCCRTCRAA \
    --p-n-jobs 8 \
    --p-min-length 400 \
    --p-max-length 600 \
    --o-reads rpoc_extracted.qza
# export fasta file
qiime tools export --input-path rpoc_extracted.qza --output-path .
```

## Get taxonomic identifier for each accession number
From here modified following QIIME2_RefDB_development (ITS2) script (Dubois et al. 2022. A detailed workflow to develop QIIME2-formatted reference databases for taxonomic analysis of DNA metabarcoding data. BMC Genomic Data)

```bash
grep ">" dna-sequences.fasta | cut -d ">" -f 2 | cut -d " " -f 1 | awk -F"|" '{print $2}' | sed 's/_.*//' > AccessionNumbers
paste <(cat AccessionNumbers) <(sed '/^>/d' rpoc.fna) > AccessionNumbers_seqs_linking_table
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
# gzip -d nucl_gb.accession2taxid.gz
fgrep -w -f AccessionNumbers nucl_gb.accession2taxid > AccessionNumbers_taxids_linking_table
wc -l AccessionNumbers_taxids_linking_table 
# retrieving accessions not found in tax table
awk -F '\t' '{print $2}' AccessionNumbers_taxids_linking_table > AccessionNumbers_found_in_accession2taxid
cat AccessionNumbers AccessionNumbers_found_in_accession2taxid | sort | uniq -u > AccessionNumbers_not_found
url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=xml&id="
url+=$(paste -s -d "," AccessionNumbers_not_found)
curl $url > AccessionNumbers_not_found.xml
# extract nodes 
paste \
 <(xmllint --xpath '//TSeq_accver/node()' AccessionNumbers_not_found.xml) \
 <(xmllint --xpath '//TSeq_taxid/node()' AccessionNumbers_not_found.xml) > Missing_taxids
# merge files into single accession number/taxids linking table
awk 'BEGIN {FS=OFS="\t"} {print $2,$3}' AccessionNumbers_taxids_linking_table > AccessionNumbers_taxids_linking_table_extracted
cat AccessionNumbers_taxids_linking_table_extracted Missing_taxids > AccessionNumbers_taxids_linking_table_final
# Clean up temporary files
rm AccessionNumbers AccessionNumbers_found_in_accession2taxid AccessionNumbers_not_found AccessionNumbers_not_found.xml AccessionNumbers_taxids_linking_table AccessionNumbers_taxids_linking_table_extracted Missing_taxids
```

## Geting taxonomic lineages for each taxid

```bash
awk -F '\t' '{print $2}' AccessionNumbers_taxids_linking_table_final | sort | uniq > Taxids_uniq
wc -l Taxids_uniq 
# mkdir taxdump
# wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
# mv new_taxdump.tar.gz taxdump/
# tar -xvzf taxdump/new_taxdump.tar.gz -C taxdump
sed -i "s/\t//g" taxdump/rankedlineage.dmp
sort -t "|" -k 1b,1 taxdump/rankedlineage.dmp > taxdump/rankedlineage_sorted

# join taxonomic lineages to taxids
join -t "|" -1 1 -2 1 -a 1 Taxids_uniq taxdump/rankedlineage_sorted > Taxids_taxonomic_lineages_linking_table
wc -l Taxids_taxonomic_lineages_linking_table
# recover missing taxonomic lineages from NCBI
awk -F '|' '$2=="" {print $0}' Taxids_taxonomic_lineages_linking_table > Taxids_not_found
url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&rettype=xml&id="
url+=$(paste -s -d "," Taxids_not_found)
curl $url > Taxids_not_found.xml

# extracting missing nodes
paste \
 <(cat Taxids_not_found) \
 <(paste -d "," \
   <(xmllint --xpath "//TaxaSet/Taxon/TaxId/node()" Taxids_not_found.xml) \
   <(xmllint --xpath "//TaxaSet/Taxon/ScientificName/node()" Taxids_not_found.xml) \
   <(for i in $(cat Taxids_not_found);do paste -d "" <(echo "ZZ") <(xmllint --xpath "//TaxaSet/Taxon[./TaxId='$i']/LineageEx/Taxon[./Rank='kingdom']/ScientificName/node()" Taxids_not_found.xml);done | sed "s/ZZ//g") \
   <(for i in $(cat Taxids_not_found);do paste -d "" <(echo "ZZ") <(xmllint --xpath "//TaxaSet/Taxon[./TaxId='$i']/LineageEx/Taxon[./Rank='phylum']/ScientificName/node()" Taxids_not_found.xml);done | sed "s/ZZ//g") \
   <(for i in $(cat Taxids_not_found);do paste -d "" <(echo "ZZ") <(xmllint --xpath "//TaxaSet/Taxon[./TaxId='$i']/LineageEx/Taxon[./Rank='class']/ScientificName/node()" Taxids_not_found.xml);done | sed "s/ZZ//g") \
   <(for i in $(cat Taxids_not_found);do paste -d "" <(echo "ZZ") <(xmllint --xpath "//TaxaSet/Taxon[./TaxId='$i']/LineageEx/Taxon[./Rank='order']/ScientificName/node()" Taxids_not_found.xml);done | sed "s/ZZ//g") \
   <(for i in $(cat Taxids_not_found);do paste -d "" <(echo "ZZ") <(xmllint --xpath "//TaxaSet/Taxon[./TaxId='$i']/LineageEx/Taxon[./Rank='family']/ScientificName/node()" Taxids_not_found.xml);done | sed "s/ZZ//g") \
   <(for i in $(cat Taxids_not_found);do paste -d "" <(echo "ZZ") <(xmllint --xpath "//TaxaSet/Taxon[./TaxId='$i']/LineageEx/Taxon[./Rank='genus']/ScientificName/node()" Taxids_not_found.xml);done | sed "s/ZZ//g")) > Missing_taxonomic_lineages
# reformat tables
paste \
  <(awk -F "|" '$2!="" {print $1}' Taxids_taxonomic_lineages_linking_table) \
  <(awk -F "|" 'BEGIN {OFS=""} $2!="" {print "k__",$9,"; p__",$8,"; c__",$7,"; o__",$6,"; f__",$5,"; g__",$4,"; s__",$2}' Taxids_taxonomic_lineages_linking_table) > Taxids_taxonomic_lineages_linking_table_reformatted
paste \
  <(awk -F '\t' '{print $1}' Missing_taxonomic_lineages) \
  <(awk -F '\t' '{print $2}' Missing_taxonomic_lineages | awk -F ',' 'BEGIN {OFS=""} {print "k__",$3,"; p__",$4,"; c__",$5,"; o__",$6,"; f__",$7,"; g__",$8,"; s__",$2}') > Missing_taxonomic_lineages_reformatted
# merge files
cat Taxids_taxonomic_lineages_linking_table_reformatted Missing_taxonomic_lineages_reformatted > Taxids_taxonomic_lineages_linking_table_merged
# remove genus names in species rank annotations
paste -d "" \
  <(awk -F 's__' '{OFS=""} {print $1,"s__"}' Taxids_taxonomic_lineages_linking_table_merged) \
  <(awk -F 's__' '{print $2}' Taxids_taxonomic_lineages_linking_table_merged | cut -d " " -f 2-) > Taxids_taxonomic_lineages_linking_table_final
```

## Creating a global table gathering accession numbers, taxids, taxonomic lineages, and nucleotide sequences

```bash
join -t $'\t' -1 2 -2 1 -a 1 \
    <(sort -t $'\t' -n -k 2 AccessionNumbers_taxids_linking_table_final) \
    <(sort -t $'\t' -n -k 1 Taxids_taxonomic_lineages_linking_table_final) | \
    awk 'BEGIN {FS=OFS="\t"} {print $2, $1, $3}' > AccessionNumbers_taxids_Taxonomic_lineages_linking_table
join -t $'\t' -1 1 -2 1 -a 1 \
      <(sort -t $'\t' -k 1b,1 AccessionNumbers_taxids_Taxonomic_lineages_linking_table)\
      <(sort -t $'\t' -k 1b,1 AccessionNumbers_seqs_linking_table) > Global_table
# check that columns are complete -- all should return 0
awk -F '\t' '{print $1}' Global_table | grep -c "^$"       
awk -F '\t' '{print $2}' Global_table | grep -c "^$"       
awk -F '\t' '{print $3}' Global_table | grep -c "^$"       
awk -F '\t' '{print $4}' Global_table | grep -c "^$" 
# clean up temporary files
rm Missing_taxids Missing_taxonomic_lineages Missing_taxonomic_lineages_reformatted Taxids_taxonomic_lineages_linking_table_merged Taxids_not_found Taxids_not_found.xml AccessionNumbers_taxids_linking_table_final AccessionNumbers_taxids_Taxonomic_lineages_linking_table AccessionNumbers_seqs_linking_table Taxids_taxonomic_lineages_linking_table Taxids_taxonomic_lineages_linking_table_final Taxids_taxonomic_lineages_linking_table_reformatted Taxids_uniq
```

## Creating QIIME2 formatted fasta and taxonomic lineages files and importation into QIIME2

```bash
awk -F '\t' 'BEGIN {OFS=""} {print ">",$1,"\n",$4}' Global_table | sed 's/-//g' > Fasta_file
awk 'BEGIN {FS=OFS="\t"} {print $1,$3}' Global_table > Taxonomic_lineages
# add numbers to duplicate accession numbers
perl -i -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' Fasta_file
# replace duplicate entries in taxonomic lineages file
awk '{print $1}' Taxonomic_lineages | awk '{print $0 (/^/ ? "_" (++c[$1]) : "")}' | sed 's/_1//' > temp
awk -F"\t" '{print $2}' Taxonomic_lineages > temp2
paste temp temp2 > Taxonomic_lineages
# import sequences into qiime2
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path Fasta_file \
  --output-path rpoc_qiime2_db_seq_2023-2-16.qza
```

## Import taxonomy into QIIME2

```bash
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path Taxonomic_lineages \
  --output-path rpoc_qiime2_db_tax_2023-2-16.qza
```

## Train classifier

```bash
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads rpoc_qiime2_db_seq_2023-2-16.qza \
  --i-reference-taxonomy rpoc_qiime2_db_tax_2023-2-16.qza \
  --o-classifier rpoc_qiime2_classifier_2023-2-16.qza
# cleanup
rm -r assembly_summary.txt query wget-log missing all_genomes.fna Fasta_file* nucl_gb.accession2taxid duplicate_rpoc_loci rpoc.ids rpoc.fna rpoc.qza rpoc_extracted.qza dna-sequences.fasta Taxids Global_table Fasta_file temp* Taxonomic_lineages taxdump/
```



