#### amending our exported feature-table.biom file with the taxonomy
#### before running this you needs to edit the headers of your taxonomy.tsv files
#### to look like this: #OTUID  taxonomy  confidence
#### those are tab separated

### use the biom convert command that comes as part of QIIME2
biom add-metadata -i feature-table.biom -o table-with-taxonomy.biom \
--observation-metadata-fp biom-taxonomy.tsv --sc-separated taxonomy
