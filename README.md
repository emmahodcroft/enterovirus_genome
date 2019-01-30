# Enterovirus D68 Full-Genome Nextstrain Analysis

Performs a full Nextstrain analysis on full-genome Enterovirus D68.

## Quickstart 
### Setup
Download in _tab delimited format_ all samples that are Enterovirus -> Enterovirus D -> Enterovirus D68 using ViPR's search function, with sequence length min:6400, max:8000. (There should be ~520+.)
_(Using the 'full genome' tick-box will result in fewer sequences)_

Place this file in `enterovirus_genome/data` (you may need to create `data`), and include the name of that file on line 9 of the Snakefile (replacing `data/genomeEntero-30Jan18.tsv` or similar). 

Place sequences and metadata from Swedish sequences in the `data` folder, and ensure the filenames on lines 12 and 13 of the Snakefile match your own. 

If you have other sequences & metadata (not on GenBank) you'd like to add, you can include these on lines 15 and 16.

### Regions
This script will allow you to look at sequences by region as well as country. The Snakefile is already set up for this kind of analysis, and region will be automatically generated for all downloaded sequences.

*However*, you should ensure the Swedish metadata file, and any additional 'manual' files, have an additional column called 'region' with an entry for each sample ('europe' - lowercase). Otherwise, no Swedish/manual sequences will have a region. 

### Running
Navigate to the `enterovirus_genome` folder and run `snakemake "auspice/enterovirus_d68_genome_tree.json"` to do a full-genome build. Initial runs may take some time, as downloading all sequences from GenBank is slow.

All accession numbers are compared, so a sequence already included in 'Swedish' or 'manual' files will not be downloaded from GenBank.

## Reruns
This Snakefile is written to make adding new data from ViPR easier. Simply download the latest full collection of samples from ViPR (using the same instructions as above), place the new file in `data`, and replace the filename on line 9 of the Snakefile. Run `snakemake`, and the script should automatically only download and BLAST sequences with accesssion numbers that have not previously been checked (even if they were not included in the analysis). 

After adding any new sequences, the a new full Nextstrain analysis will proceed. 


# Technical Notes
## Strain names
Unlike the VP1 build, strain names are not modified in this pipeline.

## Blasting
Because only whole-genome sequences are used, no Blasting is done in this pipeline (unlike VP1)

## Reruns
This Snakefile saves a copy of the most recently run parsed, downloaded ViPR file, and uses this to decide whether an accession number is 'new.' f you delete or modify the files in the 'genbank' folder that's created, then you may trigger a completely new run.
