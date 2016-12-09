#!/usr/bin/env nextflow

params.pair1 = "${PWD}/*_R1*.fastq.gz*"
params.pair2 = "${PWD}/*_R2*.fastq.gz*"
params.output = "${PWD}"
params.threads = 15

process print_log {

"""
#!/usr/bin/env bash
echo "
AmrPlusPlus - NF ~ version 1.0.0
================================
AMR Database       : \$amr_db
Kraken Database    : \$kraken_db
Host Genome        : \$genome
Forward Reads      : ${params.pair1}
Reverse Reads      : ${params.pair2}
Number of Threads  : ${params.threads}

Program paths used:
Java               : `whereis java`
BWA:               : `whereis bwa-mem`
CoverageSampler    : `whereis csa`
Kraken             : `whereis kraken`
"
"""

}

threads = params.threads

forward_reads = Channel
		.fromPath(params.pair1)
		.map { path -> [ path.toString().replace('_R1', '_RX'), path ] }

reverse_reads = Channel
		.fromPath(params.pair2)
		.map { path -> [ path.toString().replace('_R2', '_RX'), path ] }

params.read_pairs = forward_reads
	     .phase(reverse_reads)
             	.map { pair1, pair2 -> [ extractSampleName(pair1[1]), pair1[1], pair2[1] ] }

params.read_pairs.into {
	read_files_trimmed
}

process trimmomatic_qc {
	maxForks 2
	publishDir "${params.output}/trimmomatic_output"

        input:
        set dataset_id, file(forward), file(reverse) from read_files_trimmed

        output:
        set dataset_id, file("${dataset_id}_trimmed_forward.fastq"), file("${dataset_id}_trimmed_reverse.fastq") into read_files_genome, read_files_amr, read_files_kraken
	set dataset_id, file("${dataset_id}_unpaired_forward.fastq"), file("${dataset_id}_unpaired_reverse.fastq") into unpaired_files_genome, unpaired_files_amr

        """
        java -jar ${TRIM_PATH}/trimmomatic-0.32.jar PE -threads ${threads} -phred33 ${forward} ${reverse} ${dataset_id}_trimmed_forward.fastq ${dataset_id}_unpaired_forward.fastq ${dataset_id}_trimmed_reverse.fastq ${dataset_id}_unpaired_reverse.fastq ILLUMINACLIP:${TRIM_PATH}/adapters/TruSeq3-PE.fa:2:30:10:3:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """
}

process bwa_genome_alignment {
	maxForks 1
	publishDir "${params.output}/sam_bam_output"
	
	input:
	set dataset_id, file(forward), file(reverse) from read_files_genome

	output:
	set dataset_id, file("${dataset_id}_host_alignment.sam") into genome_sam_files
	
	"""
	bwa mem -t ${threads} ${GENOME} ${forward} ${reverse} > ${dataset_id}_host_alignment.sam
	"""
}

process samtools_view_to_fastq {
	maxForks 15
	publishDir "${params.output}/sam_bam_output"

	input:
	set dataset_id, file(host_alignment) from genome_sam_files

	output:
	set dataset_id, file("${dataset_id}_nonhost_sorted.bam") into nonhost_bam_files
	set dataset_id, file("${dataset_id}_nonhost_forward.fastq"), file("${dataset_id}_nonhost_reverse.fastq") into read_files_nonhost_amr, read_files_nonhost_kraken

	"""
	samtools view -hbS ${host_alignment} > ${dataset_id}_host_alignment.bam
	samtools sort --threads ${threads} ${dataset_id}_host_alignment.bam > ${dataset_id}_host_sorted.bam
	samtools view -h -f 4 -b ${dataset_id}_host_sorted.bam > ${dataset_id}_nonhost_sorted.bam
	bamToFastq -i ${dataset_id}_nonhost_sorted.bam -fq ${dataset_id}_nonhost_forward.fastq -fq2 ${dataset_id}_nonhost_reverse.fastq
	"""
}

process bwa_amr_alignment {
	maxForks 1
	publishDir "${params.output}/amr_output"
	
	input:
	set dataset_id, file(forward), file(reverse) from read_files_nonhost_amr

	output:
	set dataset_id, file("${dataset_id}_amr_alignment.sam") into amr_sam_files
	
	"""
	bwa mem -t ${threads} ${AMR_DB} ${forward} ${reverse} > ${dataset_id}_amr_alignment.sam
	"""
}

process amr_coverage_sampler {
	maxForks 4
	publishDir "${params.output}/amr_output"
	
	input:
	set dataset_id, file(amr_sam_alignment) from amr_sam_files

	output:
	set dataset_id, file("${dataset_id}_coverage_sampler_amr.tab") into amr_csa_files

	"""
	csa -ref_fp ${AMR_DB} -sam_fp ${amr_sam_alignment} -min 100 -max 100 -skip 5 -t 80 -samples 1 -out_fp ${dataset_id}_coverage_sampler_amr.tab
	"""
}

process kraken_classification {
	maxForks 1
	publishDir "${params.output}/kraken_output"

	input:
	set dataset_id, file(forward), file(reverse) from read_files_nonhost_kraken

	output:
	set dataset_id, file("${dataset_id}_kraken.raw") into kraken_raw
	set dataset_id, file("${dataset_id}_kraken.report") into kraken_report

	"""
	kraken --preload --threads ${threads} --fastq-input --paired ${forward} ${reverse} > ${dataset_id}_kraken.raw
	kraken-report ${dataset_id}_kraken.raw > ${dataset_id}_kraken.report
	"""
}

process amr_aggregate_results {
	maxForks 1
	cache false
	publishDir "${params.output}"
	
	input:
	set dataset_id, file(amr_file) from amr_csa_files.toList()
	
	output:
	file("${params.output}/AMR_aggregated_output.csv") into aggregated_output
	
	"""
	nextflow_aggregate_results.py AMR ${AMR_ANNOT} ${amr_file} >> ${params.output}/AMR_aggregated_output.csv
	"""
}

process kraken_aggregate_results {
	maxForks 1
	cache false
	publishDir "${params.output}"
	
	input:
	set dataset_id, file(kraken_report_file) from kraken_report.toList()
	
	output:
	file("${params.output}/kraken_aggregated_output.csv") into aggregated_output
	
	"""
	nextflow_aggregate_results.py kraken ${kraken_report_file} >> ${params.output}/kraken_aggregated_output.csv
	"""
}

def extractSampleName(s) {
	ret = s =~ /\/(.+)_R/;
	basepath = ~/.+\//
  	return ret[0][1] - basepath;
}

