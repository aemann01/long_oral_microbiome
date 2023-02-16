## QIIME2 denoising with DADA2

### QIIME2 conda environment installation and setup
Running on pickles

```bash
# Download and install QIIME2 conda environment
conda update -n base -c defaults conda
wget https://data.qiime2.org/distro/core/qiime2-2022.11-py38-linux-conda.yml
conda env create -n qiime2-2022.11 --file qiime2-2022.11-py38-linux-conda.yml
rm qiime2-2022.11-py38-linux-conda.yml

# Activate environment and test installation
conda activate qiime2-2022.11
qiime --help

# Enable tab completion
source tab-qiime
```

Import raw sequence data into qiime format

```bash
qiime tools import \
	--type 'SampleData[PairedEndSequencesWithQuality]' \
	--input-path sample_manifest.txt \
	--output-path sample_manifest.qza \
	--input-format PairedEndFastqManifestPhred33V2
```

### Trim primers

```bash
qiime cutadapt trim-paired \
	--i-demultiplexed-sequences sample_manifest.qza \
	--p-cores 8 \
	--p-front-f MAYGARAARMGNATGYTNCARGA \
	--p-front-r GMCATYTGRTCNCCRTCRAA \
	--verbose \
	--o-trimmed-sequences primer-trimmed.qza | tee primer_trimming.stats 

# summary statistics and quality after primer trimming
qiime demux summarize --i-data primer-trimmed.qza --o-visualization primer-trimmed-quality.qzv
```

### Denoising

```bash
time qiime dada2 denoise-paired \
  --i-demultiplexed-seqs primer-trimmed.qza \
  --p-trunc-len-f 275\
  --p-trunc-len-r 275\
  --p-max-ee-f 4\
  --p-max-ee-r 6\
  --p-min-overlap 10 \
  --p-pooling-method 'pseudo' \
  --p-chimera-method 'pooled' \
  --p-n-threads 8 \
  --o-representative-sequences rep_seq.qza \
  --o-table feature-table.qza \
  --o-denoising-stats dada2-stats.qza 1>dada2.out 2>dada2.err

# run statistics
qiime metadata tabulate \
  --m-input-file dada2-stats.qza \
  --o-visualization dada2-stats-summ.qzv

# summary
qiime feature-table summarize \
  --i-table feature-table.qza \
  --m-sample-metadata-file metadata-1.26.23.txt \
  --o-visualization feature-table-summ.qzv

qiime feature-table tabulate-seqs \
  --i-data rep_seq.qza \
  --o-visualization asv-sequences-summ.qzv
  ```

