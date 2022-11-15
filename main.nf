#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/**
===============================
IKMB GWAS Association Testing Pipeline
===============================

This Pipeline performs Association Testing for GWAS on QC'ed (and imputed) chromosome-wise genotyping data in vcf.gz format.

### Homepage / git
git@github.com:ikmb/gwas-assoc.git

**/

// Pipeline version

params.version = workflow.manifest.version

run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"




include { assoc } from './workflows/main' params(params)

workflow {

	assoc()

}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}

