# project-covid19_subgenomic
Scripts for identification of SARS-CoV-2 subgenomic RNA reads from BAM alignment files and subsequent statistical analyses.

## Description

These scripts contain the methods used for identification of SARS-CoV-2 subgenomic RNA reads from BAM alignment files and perform the subsequent statistical analyses in the study: Comparative subgenomic mRNA profiles of SARS-CoV-2 Alpha, Delta and Omicron BA.1, BA.2 and BA.5 sub-lineages using Danish COVID-19 genomic surveillance data. 

Note that this source code is shared for information only with the assumption that the user has the required input WGS and clinical metadata in hand.
Raw reads have been deposited in the GISAID's EpiCoV database under the EPI ISL identifiers listed in Supplementarytable3.csv
Alignments were performed against the SARS-CoV-2 reference with GenBank accession number MN908947.3 
Scripts have been tested using R 3.6.3

## File manifest

extract_mate_reads.R          : Script to collect all subgenomic RNA reads, outputs them into a table by position  

get_readstats.R               : Collects the number of mapped reads

collect_site_occ_w_pseudoct.R : Aggregates reads by subgenomic RNA open reading frames (sites) and by lineage groups (variant) 

do_statistical_analysis.R     : Scripts used to perform the analysis and to generate the figures (Note: the code has been provided for information only, in a prototyping shape and without metadata inputs

Supplementarytable1.csv       : Sample list with accession ids.
