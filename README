# $${\color{lightgreen}ENA-submitter}$$

Welcome to the ENA-submitter app.

A Python app allowing for a concurrent submission of many paired-end FASTQ reads or BAM/CRAM alignments to the European Nucleotide Archive Database (ENA); requires the webin-cli.jar app to be installed. 

The following options are available:

	-u, --user: Webin user name
	-p --password: Webin password
	-S --study_id: Webin study ID
	-n --study_name: Webin study NAME
	-I --instrument: NGS instrument name, e.g., 'Illumina NovaSeq 6000', 'Illumina MiSeq', 'Illumina iSeq 100'
	-i, --insert_size: The number of sequenced nucleotide in two reads forming a pair (optional)
	-s, --source:, NGS library source, e.g., 'GENOMIC', 'GENOMIC SINGLE CELL', 'TRANSCRIPTOMIC', 'TRANSCRIPTOMIC SINGLE CELL'
	-c, --library_selection: Method of the NGS library selection (default: RANDOM)
	-t, --library_strategy: NGS library type, e.g., 'WGS' (genome), 'WXS' (exome), 'RNA-Seq' (transcriptome), 'Targeted-Capture','miRNA-Seq'
	-f, --reads_dir: defines the path to the directory containing FASTQ/BAM/CRAM files.
	-r, --reads_type: determines the type of NGS reads. It should be: 'FASTQ', 'BAM' or 'CRAM'.
	-a, --Webin_accessions: defines the path to a tab-delimited txt file with sample accession numbers (samples aliases should be the same as sample prefixes in the corresponding FASTQ/BAM/CRAM file names). The sample aliases/prefixes must not contain underscore characters '_'.
	-T, --threads: number of CPU threads to use (default: 1)
	-h, --help: prints this help message.
	-v, --version: prints the version of this program.

INFO: In order to run this program one has to provide a path to the ENA webin-cli app (the bash wrapper script for the webin_cli.jar file).

INFO: This program expects all the fastq.gz reads for a given sample from the same end to be merged into a single fastq.gz file. Thus, all the first reads (R1) of the pair should be merged into one fastq.gz file, while the second fastq.gz file should contain all the second reads (R2) of the pair-end reads. Single-end reads are currently unsupported by this app.

INFO: Sample-determining prefixes in the FASTQ/BAM/CRAM file names should be placed at the beginning of the file name. These aliases/prefixes must not contain underscore characters '_'.
