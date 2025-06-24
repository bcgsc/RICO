#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
binaries
*/

//Run STAR/RSEM pipeline to generate TPM values for hg38_no_alt Ens100 annotations
process STAR_RSEM {
	tag "$library_id"
	publishDir "${params.out_dir}/${library_id}/STAR_RSEM/", mode: 'copy'

	input:
		tuple val(library_id), path(fastq1), path(fastq2)

	output:
		// tuple val("${library_id}"), path("${library_id}.STAR.genome.bam"), emit: bam
		tuple val("${library_id}"), path("${library_id}.isoforms.results"), emit: isoform
		tuple val("${library_id}"), path("${library_id}.genes.results"), emit: genes

	script:
		reference_name = "/rico/ref/hg38_no_alt"
		STAR_path = "/rico/STAR-2.5.2b/bin/Linux_x86_64/"
	"""
        /rico/RSEM-1.3.0/rsem-calculate-expression \
			--star \
			--no-bam-output \
			--star-path ${STAR_path} \
			--star-gzipped-read-file \
			--star-output-genome-bam \
			--estimate-rspd \
			--paired-end \
			--seed 12345 \
			-p 48 \
			--forward-prob 0 \
			${fastq1} ${fastq2} \
			${reference_name} \
			${library_id}
  	"""
}

process prep_files_for_immunedeconv {
    tag "$library_id"
	publishDir "${params.out_dir}/${library_id}/immunedeconv/", mode: 'copy'

	input:
		tuple val(library_id), path(gene_counts)

	output:
		tuple val(library_id), path("hgnc_tpm.csv")

	script:
		ens_HGNC_file = "/rico/ref/ens_HGNC_production_sorted.txt"
	"""
		#need to put this somewhere that travels with the code
		awk 'NR>1 { print \$1, \$(NF-1) }' ${gene_counts} | sort > ensg_tpm.csv
		echo "HGNC TPM" > header.txt
		join ensg_tpm.csv ${ens_HGNC_file} | awk '{ print \$3, \$2 }' | cat header.txt - > hgnc_tpm.csv
	"""
}

process extract_matrix {
	tag "$library_id"

	input:
		val(tmp)

	output:
		path("MATRIX.csv")

	script:
	"""
		cp /rico/ref/MATRIX.csv .
	"""
}

process immunedeconv {
	tag "$library_id"
	publishDir "${params.out_dir}/${library_id}/immunedeconv/", mode: 'copy'

	input:
		path(immunedeconv_script)
		tuple val(library_id), path(gene_counts)
		val(method)
		path(matrix_file)

	output:
		tuple val(library_id), path("immunedeconv_out_${method}.csv")

	script:
	"""
		Rscript --vanilla ${immunedeconv_script} ${gene_counts} ${method} ${matrix_file}
		mv immunedeconv_out.csv immunedeconv_out_${method}.csv
	"""
}

// Currently Un-used
process immunedeconv_cibersort {
	tag "$library_id"
	publishDir "${params.out_dir}/${library_id}/immunedeconv/", mode: 'copy'

	input:
		path(immunedeconv_cibersort_script)
		tuple val(library_id), path(gene_counts)
		path(cibersort_binary)
		path(cibersort_path)
		path(matrix_file)

	output:
        tuple val(library_id), path("immunedeconv_out_cibersort.csv")

	script:
	"""
		Rscript --vanilla ${immunedeconv_cibersort_script} ${gene_counts} ${cibersort_binary} ${cibersort_path} ${matrix_file}
		mv immunedeconv_out.csv immunedeconv_out_cibersort.csv
	"""
}

// Runs cibersort outside of immunedeconv
// immunedeconv approach could not reproduce production
// https://github.com/omnideconv/immunedeconv/discussions/138
process run_cibersort {
	tag "$library_id"
	publishDir "${params.out_dir}/${library_id}/cibersort/", mode: 'copy'

	input:
		tuple val(library_id), path(gene_counts)
		path(cibersort_binary)
		path(cibersort_LM22)
		path(matrix_file)

	output:
        tuple val(library_id), path("abundance_results.txt")

	script:
		cibersort_runner = "/rico/run_cibersort.R"
	"""
		Rscript --vanilla ${cibersort_runner} -s ${cibersort_LM22} -m ${gene_counts} -p 500 -q -a -t ${cibersort_binary}  -M ${matrix_file}
	"""
}

// Update the original LM22.txt file to use gene names mapping the ens100 reference
process update_LM22 {
	cpus 1
	memory = "1 GB"

	input:
		path(LM22_file)

	output:
        path("LM22_updated.txt"), emit: updated_LM22

	script:
		gene_mappings = "/rico/ref/Gene_name_mapping.txt"
	"""
		awk 'NR>1{ printf \"sed -i 's/%s/%s/' $LM22_file\\n\", \$1, \$2}' $gene_mappings | bash
		cp ${LM22_file} LM22_updated.txt
	"""
}


// Run IPASS from Chelsey Mayoh's group
process IPASS {
	tag "$library_id"
	publishDir "${params.out_dir}/${library_id}/IPASS/", mode: 'copy'

	input:
		tuple val(library_id), path(gene_counts)
		path(matrix_file)

	output:
		tuple val(library_id), path("IPASS_Results.txt")

	script:
		IPASS_script = "/rico/IPASS_Calculation.R"
		IPASS_RData = "/rico/ref/IPASS_hg38.RData"
	"""
		Rscript --vanilla ${IPASS_script} ${gene_counts} ${matrix_file} ${IPASS_RData}
	"""
}

// Run the M1M2 score from BCGSC
process M1M2 {
	tag "$library_id"
	publishDir "${params.out_dir}/${library_id}/M1M2/", mode: 'copy'

	input:
		tuple val(library_id), path(gene_counts)
		path(matrix_file)

	output:
		tuple val(library_id), path("hgnc_tpm.csv_M1M2.txt")

	script:
		M1M2_script = "/rico/M1M2.R"
	"""
		Rscript --vanilla ${M1M2_script} ${gene_counts} ${matrix_file}
	"""
}
