#!/bin/bash

#nextflow get-gemma-de-results.nf --level class -resume --outdir /space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma-de-sc-pipeline-class --filter_genes false
#nextflow get-gemma-de-results.nf --level subclass -resume --outdir /space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma-de-sc-pipeline-subclass --filter_genes false
#nextflow get-gemma-de-results.nf --author_submitted true -resume --outdir /space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma-de-author --filter_genes false
nextflow gemma-vs-gemma-pseudobulk-de.nf -params-file params.class.json -resume --outdir /space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma-vs-gemma-class
nextflow gemma-vs-gemma-pseudobulk-de.nf -params-file params.subclass.json -resume --outdir /space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma-vs-gemma-subclass
