## importing all the demuxed reads back into qiime2
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ../metadata/manifest.tsv \
--input-format PairedEndFastqManifestPhred33V2 --output-path demux-paired-end.qza

## visualize demuxed reads
qiime demux summarize --i-data demux-paired-end.qza --o-visualization demux-paired-end.qzv

## DADA2
qiime dada2 denoise-paired --i-demultiplexed-seqs demux-paired-end.qza \
--p-n-threads 12 --o-table table.qza \
--p-trim-left-f 6 --p-trunc-len-f 247 \
--p-trim-left-r 6 --p-trunc-len-r 247 \
--o-representative-sequences rep-seqs.qza --output-dir DADA2

## classifying
qiime feature-classifier classify-sklearn --i-classifier /s/angus/h/nobackup/databases/QIIME2_classifiers/qiime2-2020.11/GG/gg_13_8$
--i-reads rep-seqs.qza --o-classification taxonomy.qza

## making the no-chloro-no-mito directory
mkdir no-chloro-no-mito

## copy the taxonomy file into it
cp taxonomy.qza no-chloro-no-mito/

## filtering ASVs
qiime taxa filter-table --i-table table.qza --i-taxonomy taxonomy.qza \
--p-exclude mitochondria,chloroplast --o-filtered-table no-chloro-no-mito/table.qza

## filtering rep-seqs
qiime taxa filter-seqs --i-sequences rep-seqs.qza --i-taxonomy taxonomy.qza \
--p-exclude mitochondria,chloroplast --o-filtered-sequences no-chloro-no-mito/rep-seqs.qza

qiime alignment mafft --i-sequences no-chloro-no-mito/rep-seqs.qza \
--o-alignment no-chloro-no-mito/aligned-rep-seqs.qza
qiime alignment mask --i-alignment no-chloro-no-mito/aligned-rep-seqs.qza \
--o-masked-alignment no-chloro-no-mito/masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment no-chloro-no-mito/masked-aligned-rep-seqs.qza \
--o-tree no-chloro-no-mito/unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree no-chloro-no-mito/unrooted-tree.qza \
--o-rooted-tree no-chloro-no-mito/rooted-tree.qza

## change into no-chloro-no-mito
cd no-chloro-no-mito

## make phyloseq dir
mkdir exported

## export to phyloseq dir
qiime tools export --input-path rep-seqs.qza --output-path exported/
qiime tools export --input-path taxonomy.qza --output-path exported/
qiime tools export --input-path rooted-tree.qza --output-path exported/
qiime tools export --input-path table.qza --output-path exported/

## need to to alter the headers of taxonomy then amend to the biom file for import into phyloseq
