#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/**
===============================
Pipeline
===============================

This Pipeline performs ....

### Homepage / git
git@github.com:ikmb/pipeline.git

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

