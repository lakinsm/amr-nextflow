#!/usr/bin/env nextflow

params.pair1 = "${PWD}/*_R1*.fastq"
params.pair2 = "${PWD}/*_R2*.fastq"
params.genome = "/s/angus/index/databases/bwa_indexes/mod_bos_taurus/mod_bos_taurus.fna"
params.amr_db = "/s/angus/index/databases/MEGARes/megares_database.fasta"
params.kraken_db = "/s/angus/index/databases/kraken_databases/kraken_refcustom_Oct2016/database"
params.threads = 1

log.info "AmrPlusPlus - NF ~ version 1.0.0"
log.info "================================"
log.info "AMR Database       : ${params.amr_db}"
log.info "Kraken Database    : ${params.kraken_db}"
log.info "Host Genome        : ${params.genome}"
log.info "Forward Reads      : ${params.pair1}"
log.info "Reverse Reads      : ${params.pair2}"
log.info "Number of Threads  : ${params.threads}"
log.info "\n"

genome = file(params.genome)
amr_db = file(params.amr_db)
kraken_db = file(params.kraken_db)
threads = params.threads

if(!genome.exists()) {
	exit 1, "Unable to find genome file: {params.genome}"
}
if(!amr_db.exists()) {
        exit 1, "Unable to find genome file: {params.amr_db}"
}
if(!kraken_db.exists()) {
        exit 1, "Unable to find virulence file: {params.kraken_db}"
}

forward_reads = Channel
		.fromPath(params.pair1)
		.map { path -> [ path.toString().replace('_R1', '_RX'), path ] }

reverse_reads = Channel
		.fromPath(params.pair2)
		.map { path -> [ path.toString().replace('_R2', '_RX'), path ] }

params.read_pairs = forward_reads
	     .phase(reverse_reads)
             .map { pair1, pair2 -> [ pathToDatasetID(pair1[1]), pair1[1], pair2[1] ] }

params.read_pairs.into {
	read_files_genome;
	read_files_amr;
	read_files_kraken
}

process bowtie2_genome_alignment {
	input:
	set dataset_id, file(forward), file(reverse) from read_files_genome
	file genome

	output:
	set dataset_id, file('host_alignment.sam') into genome_sam_files
	set dataset_id, file('nonhost_alignment.bam') into genome_sam_files
	set dataset_id, file('nonhost_forward.fastq') into read_files_nonhost
	set dataset_id, file('nonhost_reverse.fastq') into read_files_nonhost

	"""
	bowtie2 -p ${threads} -x ${genome} -1 ${forward} -2 ${reverse} -S host_alignment.sam
	samtools view -h -f 4 -b host_alignment.sam > nonhost_alignment.bam
	bamToFastq -i nonhost_alignment.bam -fq nonhost_forward.fastq -fq2 nonhost_reverse.fastq
	"""
}

process bowtie2_amr_alignment {
	input:
	set dataset_id, file(forward), file(reverse) from read_files_nonhost
	file index from amr_db.first()

	output:
	set dataset_id, file('amr_alignment.sam') into amr_sam_files

	"""
	bowtie2 -p ${threads} -x amr.index -1 $forward -2 $reverse -S amr_alignment.sam
	"""
}

process amr_coverage_sampler {
	input:
	set dataset_id, file(amr_sam_alignment) from amr_sam_files
	file amrdb from amr_db

	output:
	set dataset_id, file('coverage_sampler_amr.tab') into amr_csa_files

	"""
	csa -ref_fp $amrdb -sam_fp $amr_sam_alignment -min 100 -max 100 -skip 5 -t 80 -samples 1 -out_fp coverage_sampler_amr.tab
	"""
}

process kraken_classification {
	input:
	set dataset_id, file(forward), file(reverse) from read_files_nonhost
	file kdb from kraken_db

	output:
	set dataset_id, file('kraken.raw') into kraken_files
	set dataset_id, file('kraken.report') into kraken_files

	"""
	kraken --preload --db ${kdb} --threads ${threads} --fastq-input --paired ${forward} ${reverse} > kraken.raw
	kraken-report -db ${kdb} kraken.raw > kraken.report
	"""
}

def pathToDatasetID(path) {
  	return path.getParent().toString();
}
