# Treponema revision analyses

Comment from reviewer: How specific is rpoC at distinguishing different Treponema species? It is curious that a non-oral treponeme comes up with high importance.

I want to pull all of the treponema rpoC genes from our database and do some ANI and phylogenetic analyses -- maybe rpoC isn't particularly good at distinguishing this genus?

```bash
# working on Hillary
cd /home/allie/long_oral_microbiome/07-Revisions-10.7.24
grep "Treponema" ~/domhain/2022-HIV_oral_microbiome/00-database_build/rpoc_ref.fa -A 1 > treponemes.fa
# alignment
mafft --adjustdirection treponemes.fa > treponemes.align.fa
# fix headers for raxml
sed -i 's/|/ /' treponemes.align.fa
rm *tre
# build tree
raxmlHPC-PTHREADS-SSE3 -T 40 \                                       
    -m GTRCAT \
    -c 25 \
    -e 0.001 \
    -p 31514 \
    -f a \
    -N 100 \
    -x 02938 \
    -n treponeme.tre -s treponemes.align.fa
# visualize in figtree -- do species cluster in expected ways?
grep ">" treponemes.align.fa | sed 's/>//' | awk '{print $1, $3 "_" $4}' | sed 's/ /\t/' | sed '1i accession\tspecies' > ann
otations.txt
```

The tree looks really good -- let's try to place our T. phagedenis ASV to make sure that it properly clusters using EPA

```bash
cat treponemes.fa t.phagendenis.fa > query.fa
mafft --adjustdirection query.fa > query.align.fa
sed -i 's/|/ /' query.align.fa

raxmlHPC-PTHREADS-SSE3 -f v -G 0.2 -m GTRCAT -n treponeme.epa.tre -s query.align.fa -t RAxML_bestTree.treponeme.tre -T 10

sed 's/QUERY___//g' RAxML_labelledTree.treponeme.epa.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placementTree.treponeme.epa.tre
```
