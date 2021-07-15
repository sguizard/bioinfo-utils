#! /usr/bin/env Rscript

args<-commandArgs(TRUE)

## DESCRIPTION
# This script link data beetween isoseq3 cluster report, and TAMA collapse/merge reports.
# This way, it's possible to know which ccs is in which transcript and which transcript eis in which gene model
# Options are this following
# args[1] ==> *.ccs.clustered.cluster_report.csv
# args[2] ==> .tc_trans_read_up.bed
# args[3] ==> merged_*.subreads_merge.txt
# args[4] ==> prefix

## _up.bed generation (We add a suplemental to TAMA collpase output for easier linking with TAMA merge)
# mkdir 8-link_ccs_to_models
# cd 8-link_ccs_to_models
# scp -r "eddie.ecdf.ed.ac.uk:~/gp_vol/Projects/isoseq_cluster/8-link_ccs_to_models/data/*" .
# ls|grep -v 'merged_'|tr '.' ' '|cut -d ' ' -f 1|uniq|xargs mkdir
# mv merged_m64036_210418_013328.subreads_merge.txt m64036_210418_013328/
# mv merged_m64036_210415_040208.subreads_merge.txt m64036_210415_040208/
# mv merged_m64036_210413_213334.subreads_merge.txt m64036_210413_213334/
# mv merged_m64036_210320_212651.subreads_merge.txt m64036_210320_212651/
# mv merged_m64036_210412_151836.subreads_merge.txt m64036_210412_151836/
# mv merged_m64036_210322_035446.subreads_merge.txt m64036_210322_035446/
# mv merged_m64036_210319_151149.subreads_merge.txt m64036_210319_151149/
# mv merged_m64036_210315_182820.subreads_merge.txt m64036_210315_182820/
#
# for f in (ls)
#     cd $f
#     for i in *_tc_trans_read.bed
#             set n (string replace "_trans_read" "" $i)
#             set o (string replace ".bed" "_up.bed" $i)
#             set cmd "perl -F'\t' -nae 'chomp; \$F[3] =~ /(G\d+\.\d+);/; print \"\$_\t"$n"_\$1\n\"' $i > $o"
#             eval $cmd
#     end
#     set out $f"_tc_trans_read_up.bed"
#     cat *_up.bed > $out
#     cd /home/sguizard/Work/Projects/isoseq/cluster/8-link_ccs_to_models
# end

library(tidyverse)

read_bed13 <- function(file) {
    require(readr)
    read_tsv(
        file,
        col_names = c(
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "blockCount",
            "blockSizes",
            "blockStarts",
            "tc_id"),
        col_types = cols(
            chrom       = col_character(),
            chromStart  = col_integer(),
            chromEnd    = col_integer(),
            name        = col_character(),
            score       = col_double(),
            strand      = col_character(),
            thickStart  = col_integer(),
            thickEnd    = col_integer(),
            itemRgb     = col_character(),
            blockCount  = col_integer(),
            blockSizes  = col_character(),
            blockStarts = col_character(),
            tc_id       = col_character()))
}

read_bed12 <- function(file) {
    require(readr)
    read_tsv(
        file,
        col_names = c(
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "blockCount",
            "blockSizes",
            "blockStarts"),
        col_types = cols(
            chrom       = col_character(),
            chromStart  = col_integer(),
            chromEnd    = col_integer(),
            name        = col_character(),
            score       = col_double(),
            strand      = col_character(),
            thickStart  = col_integer(),
            thickEnd    = col_integer(),
            itemRgb     = col_character(),
            blockCount  = col_integer(),
            blockSizes  = col_character(),
            blockStarts = col_character()))
}

ccs <- read_csv(args[1])
bed <- 
    read_bed13(args[2]) %>% 
    select(name, tc_id) %>%  
    separate(col = name, into = c("gene_id", "cluster_id"), sep = ";")
mtc <- read_bed12(args[3])

r2g <- bed %>%   
    right_join(ccs) %>%  
    select(read_id, cluster_id, tc_id) %>%  
    left_join( 
        mtc %>%  
        separate(col = name, into = c("gene_merge", "tc_id"), sep = ";", remove = FALSE) %>%  
        select(chrom, chromStart, chromEnd, name, gene_merge, tc_id))

out <- paste0(args[4], "_ccs2g.tsv", sep = "")
r2g %>% write_tsv(out)
