#!/usr/bin/env nextflow

params.pair1 = "${PWD}/*_R1*.fastq"
params.pair2 = "${PWD}/*_R2*.fastq"
params.threads = 1

process print_log {

'''
echo "
AmrPlusPlus - NF ~ version 1.0.0
================================
AMR Database       : \${amr_db}
Kraken Database    : \${kraken_db}
Host Genome        : \${genome}
Forward Reads      : ${params.pair1}
Reverse Reads      : ${params.pair2}
Number of Threads  : ${params.threads}
\n"

Program paths used:
Java               : `whereis java`
Bowtie2:           : `whereis bowtie2`
CoverageSampler    : `whereis csa`
Kraken             : `whereis kraken`
\n"
'''

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
             	.map { pair1, pair2 -> [ pathToDatasetID(pair1[1]), pair1[1], pair2[1] ] }

params.read_pairs.into {
	read_files_trimmed
}

process trimmomatic_qc {
        input:
        set dataset_id, file(forward), file(reverse) from read_files_trimmed

        output:
        set dataset_id, file('trimmed_forward.fastq'), file('trimmed_reverse.fastq') into read_files_genome, read_files_amr, read_files_kraken

        """
        java -jar /s/angus/index/common/tools/Trimmomatic-0-1.32/trimmomatic-0.32.jar PE -threads ${threads} -phred33 ${forward} ${reverse} trimmed_forward.fastq 1U.fastq trimmed_reverse.fastq 2U.fastq ILLUMINACLIP:/s/angus/index/common/tools/Trimmomatic-0-1.32/adapters/TruSeq3-PE.fa:2:30:10:3:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """
}

process bowtie2_genome_alignment {
	publishDir 'sam_bam_output'
	
	input:
	set dataset_id, file(forward), file(reverse) from read_files_genome
	file genome

	output:
	set dataset_id, file('host_alignment.sam') into genome_sam_files
	set dataset_id, file('nonhost_alignment.bam') into nonhost_bam_files
	set dataset_id, file('nonhost_forward.fastq'), file('nonhost_reverse.fastq') into read_files_nonhost_amr, read_files_nonhost_kraken

	"""
	bowtie2 -p ${threads} -x \${genome} -1 ${forward} -2 ${reverse} -S host_alignment.sam
	samtools view -h -f 4 -b host_alignment.sam > nonhost_alignment.bam
	bamToFastq -i nonhost_alignment.bam -fq nonhost_forward.fastq -fq2 nonhost_reverse.fastq
	"""
}

process bowtie2_amr_alignment {
	publishDir 'amr_output'
	
	input:
	set dataset_id, file(forward), file(reverse) from read_files_nonhost_amr
	file index from amr_db.first()

	output:
	set dataset_id, file('amr_alignment.sam') into amr_sam_files

	"""
	bowtie2 -p ${threads} -x \${amr_db} -1 $forward -2 $reverse -S amr_alignment.sam
	"""
}

process amr_coverage_sampler {
	publishDir 'amr_output'
	
	input:
	set dataset_id, file(amr_sam_alignment) from amr_sam_files
	file amrdb from amr_db

	output:
	set dataset_id, file('coverage_sampler_amr.tab') into amr_csa_files

	"""
	csa -ref_fp \$amrdb -sam_fp $amr_sam_alignment -min 100 -max 100 -skip 5 -t 80 -samples 1 -out_fp coverage_sampler_amr.tab
	"""
}

process kraken_classification {
	publishDir 'kraken_output'

	input:
	set dataset_id, file(forward), file(reverse) from read_files_nonhost_kraken
	file kdb from kraken_db

	output:
	set dataset_id, file('kraken.raw') into kraken_raw
	set dataset_id, file('kraken.report') into kraken_report

	"""
	kraken --preload --db \${kraken_db} --threads ${threads} --fastq-input --paired ${forward} ${reverse} > kraken.raw
	kraken-report -db \${kraken_db} kraken.raw > kraken.report
	"""
}

def pathToDatasetID(path) {
  	return path.getParent().toString();
}
