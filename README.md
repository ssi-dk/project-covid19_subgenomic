# project-covid19_subgenomic
Scripts for identification of SARS-CoV-2 subgenomic RNA reads from BAM alignment files and subsequent statistical analyses.

## Description

These scripts contain the methods used for identification of SARS-CoV-2 subgenomic RNA reads from BAM alignment files and perform the subsequent statistical analyses in the study: Comparative subgenomic mRNA profiles of SARS-CoV-2 Alpha, Delta and Omicron BA.1, BA.2 and BA.5 sub-lineages using Danish COVID-19 genomic surveillance data. 

Note that this source code is shared for information only with the assumption that the user has the required input WGS and clinical metadata in hand.
Raw reads have been deposited in the GISAID's EpiCoV database under the EPI ISL identifiers listed in Supplementarytable1.csv
Alignments were performed against the SARS-CoV-2 reference with GenBank accession number MN908947.3 
Scripts have been tested using R 3.6.3

## File manifest

Scripts have been tested using R 3.6.3 with the following package versions:
Rsamtools 2.8.0, IRanges 2.26.0, glmnet 4.1-6, caret 6.0-93, ggplot2 3.4.0, ggpubr 0.5.0, peakPick 0.11, ggfortify 0.4.14, cluster 2.1.4

extract_mate_reads.R          : Script to collect all subgenomic RNA reads, outputs them into a table by position  
                                Input: BAM alignement files, metadata file with sampleids and path to BAM files, WGS metadata with mean coverage information 
                                Output: Table of subgenomic RNA reads with the following columns qname|flag|rname|strand|pos|cigar|mpos|isize|seq|sample|leaderpos|mappedreads|substr_match 
 
collect_site_occ_w_pseudoct.R : Aggregates reads by subgenomic RNA open reading frames (sites) and by lineage groups (variant) 
                                Input: Table of subgenomic RNA reads
                                Output: Table of occurrences with the folling columns: sample|pos|count|leader_depth|lineage|mapped

get_readstats.R               : Collects the number of mapped reads
                                Input: BAM alignement files, metadata file with sampleids and path to BAM files, WGS metadata with mean coverage information 
                                Output: Table with the number of mapped reads in each sample

do_statistical_analysis.R     : Scripts used to perform the analysis and to generate the figures (Note: the script does not include any example patient metadata but describes the expected data columns to be imported.
                                
Supplementarytable1.csv       : Sample list with accession ids.
