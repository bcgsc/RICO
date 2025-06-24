#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {STAR_RSEM} from './rico_processes.nf'
include {immunedeconv as quantiseq} from './rico_processes.nf'
include {immunedeconv as TIMER} from './rico_processes.nf'
include {immunedeconv as MCPCounter} from './rico_processes.nf'
include {immunedeconv as epic} from './rico_processes.nf'
include {immunedeconv as xcell} from './rico_processes.nf'
include {prep_files_for_immunedeconv} from './rico_processes.nf'
include {immunedeconv_cibersort} from './rico_processes.nf'
include {extract_matrix} from './rico_processes.nf'
include {run_cibersort} from './rico_processes.nf'
include {update_LM22} from './rico_processes.nf'
include {IPASS} from './rico_processes.nf'
include {M1M2} from './rico_processes.nf'

/*
* pipeline input parameters
*/

    params.release = "0.4.0"

	log.info """\
    RICO Analysis Pipeline
    ===================================
    version         : ${params.release}
    samples_file    : ${params.samples_file}
    out_dir         : ${params.out_dir}
    """
	.stripIndent()

    // Using explicit paths that can be found by immunedeconv container
    immunedeconv_script = "$workflow.projectDir/scripts/immunedeconv.R"

//main workflow
workflow {

    //Load in the samples file
    samples = Channel
        .fromPath(params.samples_file)
        .splitCsv(header:true)
        .map{ row-> tuple(row.library_id, row.fastq1, row.fastq2) }

    //Will copy the matrix out of 1 container to make the data available to the other
    extract_matrix("tmp")

    //update the orignal LM22 from cibersort to match ens100 ids
    if (params.cibersort_script && params.cibersort_matrix) {
        update_LM22(params.cibersort_matrix)
    }

    STAR_RSEM(samples)
    prep_files_for_immunedeconv(STAR_RSEM.out.genes)

    //run the immunedeconv cell populaton estimators
    quantiseq(immunedeconv_script, prep_files_for_immunedeconv.out, "quantiseq" , extract_matrix.out)
    MCPCounter(immunedeconv_script, prep_files_for_immunedeconv.out, "mcp_counter" , extract_matrix.out)
    epic(immunedeconv_script, prep_files_for_immunedeconv.out, "epic" , extract_matrix.out)
    xcell(immunedeconv_script, prep_files_for_immunedeconv.out, "xcell" , extract_matrix.out)

    //run cibersort outside of immunedeconv (results differ between standalone and immunedeconv)
    if (params.cibersort_script && params.cibersort_matrix) {
        run_cibersort(prep_files_for_immunedeconv.out, params.cibersort_script, update_LM22.out.updated_LM22, extract_matrix.out)
    }

    //run immune specific scoring tools
    M1M2(prep_files_for_immunedeconv.out, extract_matrix.out)
    IPASS(prep_files_for_immunedeconv.out, extract_matrix.out)
}

//will run at the end of the analysis.
workflow.onComplete {
	//final_results_ch.view{ println "[complete] $it"}
	println "Pipeline $workflow.scriptName completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
