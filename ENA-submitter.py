#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Łukasz Szafron
@email: lukszafron@gmail.com
"""
appname = "ENA-submitter"
app_version = "1.0"

import pandas as pd
import time, datetime, subprocess, re, os, sys, getopt, multiprocessing, pathlib

start_time = time.time()
def usage():
    print(
        "\nWelcome to the " + appname + " app.\n\n"
        "The following options are available:\n\n"
        "\t-u, --user: Webin user name\n"
        "\t-p --password: Webin password\n"
        "\t-S --study_id: Webin study ACCESSION no.\n"
        "\t-n --study_name: Webin study NAME\n"
        "\t-I --instrument: NGS instrument name, e.g., 'Illumina NovaSeq 6000', 'Illumina MiSeq', 'Illumina iSeq 100'\n"
        "\t-i, --insert_size: The number of sequenced nucleotide in two reads forming a pair (optional)\n"
        "\t-s, --source:, NGS library source, e.g., 'GENOMIC', 'GENOMIC SINGLE CELL', 'TRANSCRIPTOMIC', 'TRANSCRIPTOMIC SINGLE CELL'\n"
        "\t-c, --library_selection: Method of the NGS library selection (default: RANDOM)\n"
        "\t-t, --library_strategy: NGS library type, e.g., 'WGS' (genome), 'WXS' (exome), 'RNA-Seq' (transcriptome), 'Targeted-Capture', 'miRNA-Seq'\n"
        "\t-f, --reads_dir: defines the path to the directory containing FASTQ/BAM/CRAM files.\n"
        "\t-r, --reads_type: determines the type of NGS reads. It should be: 'FASTQ', 'BAM' or 'CRAM'.\n"
        "\t-F, --FASTQ_type: defines whether the submitted FASTQ files contain paired-end 'paired' or single reads 'single'. The default value is: 'paired'.\n"
        "\t-x, --suffix: a suffix determining the type of the BAM/CRAM files used. This is a character string separating sample names from the file extension (.bam or .cram).\n"
        "\t-a, --Webin_accessions: defines the path to a tab-delimited txt file with sample accession numbers (samples aliases should be the same as sample prefixes in the corresponding FASTQ/BAM/CRAM file names). The sample aliases/prefixes must not contain underscore characters '_'.\n"
        "\t-T, --threads: number of CPU threads to use (default: 1)\n"
        "\t-h, --help: prints this help message.\n"
        "\t-v, --version: prints the version of this program.\n"
        )
# The next two lines are for debugging purposes only and should be commented in the final program.
# option_list = ["-u", "******", "-p", "******", "-S", "PRJEB61419", "-n", "44genes", "-I", "Illumina MiSeq", "-s", "GENOMIC", "-t", "WXS","-a", "/workspace/lukasz/NGS-all-in-one/RUNS/44GENESVEP100/MAPPINGS_TRIMMED/GRCh38.chr.only.genome/HISAT2/CRNDE-44g_capture_targets_hg38.bed_and_SureSelect_All_Exon_V7_hg38_Padded.bed.subset/Webin-accessions-2023-04-18T09 37 05.538+01 00.txt", "-f", "/workspace/lukasz/NGS-all-in-one/RUNS/44GENESVEP100/MAPPINGS_TRIMMED/GRCh38.chr.only.genome/HISAT2/CRNDE-44g_capture_targets_hg38.bed_and_SureSelect_All_Exon_V7_hg38_Padded.bed.subset/", "-r", "BAM", "-x" , "_sorted.no_dups", "-T", "5"]
# opts, args = getopt.getopt(option_list, "u:p:S:n:I:i:s:c:t:f:r:F:x:a:T:hv", ["user=","password=","study_id=","study_name=","instrument=","insert_size=","source=","library_selection=","library_strategy=","reads_dir=", "reads_type=", "FASTQ_type=", "suffix=", "Webin_accessions=","threads=","help","version"])

opts, args = getopt.getopt(sys.argv[1:], "u:p:S:n:I:i:s:c:t:f:r:F:x:a:T:hv", ["user=","password=","study_id=","study_name=","instrument=","insert_size=","source=","library_selection=","library_strategy=","reads_dir=", "reads_type=", "FASTQ_type=", "suffix=", "Webin_accessions=","threads=","help","version"])

try:
        opts,args
        
        if len(opts) == 0:
                usage()
                sys.exit()
        for o, a in opts:
            if o in ("-u", "--user"):
                user = a
            elif o in ("-p", "--password"):
                password = a
            elif o in ("-S", "--study_id"):
                study_id = a
            elif o in ("-n", "--study_name"):
                study_name = a
            elif o in ("-I", "--instrument"):
                instrument = a
            elif o in ("-i", "--insert_size"):
                insert_size = a
            elif o in ("-s", "--source"):
                source = a
            elif o in ("-c", "--library_selection"):
                lib_selection = a
            elif o in ("-t", "--library_strategy"):
                lib_strategy = a
            elif o in ("-f", "--reads_dir"):
                reads_dir = a
            elif o in ("-r", "--reads_type"):
                reads_type = a
            elif o in ("-F", "--FASTQ_type"):
                fastq_type = a
            elif o in ("-x", "--suffix"):
                suffix = a
            elif o in ("-a", "--Webin_accessions"):
                accessions_path = a
            elif o in ("-T", "--threads"):
                threads = a
            elif o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-v", "--version"):
                print(appname, " version: ", app_version, sep = "")
                sys.exit()
            else:
                assert False, "Unhandled option: "+o
except getopt.GetoptError as err:
    # print help information and exit:
    print("\n"+str(err)) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)
    
p = pathlib.Path("/programs/miniconda3/envs")
p = list(p.rglob("*webin*.jar"))
if(len(p) != 1):
    raise(Exception("ENA Webin-CLI executable was not found or the path to it is ambiguous."))

webin_cli_path = str(p[0])
    
try:
    user
except:
    raise(Exception("The user name is missing."))
try:
    password
except:
    raise(Exception("The password is missing."))
try:
    study_id
except:
    raise(Exception("The study accession no. is missing."))
try:
    study_name
except:
    raise(Exception("The study name is missing."))
try:
    instrument
except:
    raise(Exception("The instrument name is missing."))
try:
    insert_size
except:
    insert_size = None
try:
    source
except:
    raise(Exception("The library source argument is missing."))
try:
    lib_selection
except:
    lib_selection = "RANDOM"
try:
    lib_strategy
except:
    raise(Exception("The library strategy argument is missing."))
try:
    reads_dir
except:
    raise(Exception("The path to the FASTQ/BAM/CRAM-containing directory is missing."))
try:
    reads_type
except:
    raise(Exception("The type of the reads to be submitted to the ENA is missing (possible values are: 'FASTQ', 'BAM', 'CRAM')."))
if reads_type in ("BAM", "CRAM"):
    try:
        suffix
    except:
        raise(Exception("For BAM and CRAM files, the file name suffix must be provided."))
elif reads_type == "FASTQ":
    try:
        fastq_type
    except:
        fastq_type = "paired"
try:
    accessions_path
except:
    raise(Exception("The path to the CSV file with sample accession numbers is missing."))
try:
    threads
except:
    threads = 1
finally:
    threads = int(threads)

with open(file=accessions_path, mode="rt") as acc:
    accessions = pd.read_csv(acc, delimiter = "\t", skipfooter=1, engine="python")
    
nrows,ncols = accessions.shape
aliases = [(accessions.at[i,"ACCESSION"],accessions.at[i,"ALIAS"]) for i in range(nrows)]

class SAMPLE:
    def __init__(self, alias):
        print("Submitting the {} sample, accession_id: {} ({} out of {})...\n".format(alias[1], alias[0], aliases.index((alias[0],alias[1]))+1, len(aliases)), flush = True)
        self.alias = alias
        if reads_type == "FASTQ":
            if fastq_type == "paired":
                self.matching_reads = [reads for reads in os.listdir(reads_dir) if re.search(string=reads, pattern="^"+str(alias[1]).replace(".", "\.")+"_.*\.fastq\.gz$", flags=re.I)]
                if(len(self.matching_reads) != 2):
                    raise Exception("There is an incorrect number of merged {0} files ({1}) for the sample: {2}. All R1 reads should be provided in one merged fastq file, and all R2 reads in another merged fastq file.\n".format("FASTQ", len(self.matching_reads), self.alias[1]))
            elif fastq_type == "single":
                self.matching_reads = [reads for reads in os.listdir(reads_dir) if re.search(string=reads, pattern="^"+str(alias[1]).replace(".", "\.")+"_.*\.fastq\.gz$", flags=re.I)]
                if(len(self.matching_reads) != 1):
                    raise Exception("There is an incorrect number of merged {0} files ({1}) for the sample: {2}. All reads should be provided in one merged fastq file.\n".format("FASTQ", len(self.matching_reads), self.alias[1]))    
        elif reads_type == "BAM":
            self.matching_reads = [reads for reads in os.listdir(reads_dir) if re.search(string=reads, pattern="^"+str(alias[1]).replace(".", "\.")+suffix+"\.bam$", flags=re.I)]
            if(len(self.matching_reads) != 1):
                raise Exception("There is an incorrect number of {0} files ({1}) for the sample: {2}.\n".format("BAM", len(self.matching_reads), self.alias[1]))
        elif reads_type == "CRAM":
            self.matching_reads = [reads for reads in os.listdir(reads_dir) if re.search(string=reads, pattern="^"+str(alias[1]).replace(".", "\.")+suffix+"\.cram$", flags=re.I)]
            if(len(self.matching_reads) != 1):
                raise Exception("There is an incorrect number of {0} files ({1}) for the sample: {2}.\n".format("CRAM", len(self.matching_reads), self.alias[1]))
    def manifest(self):
        # print("Generating the MANIFEST file...\n", flush=True)
        self.manifest_dict = {
            "SAMPLE":self.alias[0], 
            "STUDY":study_id, 
            "NAME":'_'.join([study_name, str(self.alias[1])]), 
            "INSTRUMENT":instrument, 
            "INSERT_SIZE":str(insert_size), 
            "LIBRARY_SOURCE":source, 
            "LIBRARY_SELECTION":lib_selection, 
            "LIBRARY_STRATEGY":lib_strategy,
            "SUBMISSION_TOOL":appname,
            "SUBMISSION_TOOL_VERSION":app_version
            }
        self.manifest_dict = {k:v for k,v in self.manifest_dict.items() if v != "None"}
        self.items = [v for v in self.manifest_dict.items()]
        self.tsv = []
        for i,v in enumerate(self.items):
            self.tsv.append('\t'.join(self.items[i]))
        self.tsv = self.tsv + ['\t'.join([reads_type, reads]) for reads in self.matching_reads]
        self.file_alias = os.path.join(reads_dir, str(self.alias[1])+"_Manifest.txt")
        with open(file = self.file_alias, mode = "wt") as out:
            out.write('\n'.join(self.tsv))
    def submission(self):
#        print("Submitting the FASTQ/BAM/CRAM files and metadata to the European Nucleotide Archive (ENA) database...\n", flush=True)
        self.res = subprocess.run(["java", "-jar", webin_cli_path, "-context", "reads", "-inputDir", reads_dir, "-outputDir", reads_dir, "-manifest", self.file_alias, "-userName", user, "-password", password, "-submit"], text = True)
        if self.res.returncode == 0:
            print("{} sample submission is done.\n".format(self.alias[1]), flush=True)
        else:
            print("{} SAMPLE SUBMISSION HAS FAILED!!!\n".format(self.alias[1]), flush=True)

def exec_func(alias):
    sample = SAMPLE(alias)
    sample.manifest()
    sample.submission()
    return sample.res.returncode

# protect the entry point
if __name__ == '__main__':
    with multiprocessing.Pool(threads) as pool:
        returncodes = pool.map(exec_func, aliases)
        
if len([v for v in returncodes if v == 0]) == len(aliases):
    print("All samples have been submitted to the ENA successfully.\n")
else:
    print("WARNING: Some errors occurred during the submission process.\n")
    erroneous_submissions = [aliases[i][1] for i,v in enumerate(returncodes) if v != 0]
    print("{} sample(s): {} has/have not been submitted to the ENA because of errors.\n".format(len(erroneous_submissions), ', '.join(erroneous_submissions)))

print("--- The sample submission process took %s hh:mm:ss. ---\n" % str(datetime.timedelta(seconds = round(time.time() - start_time, 0))))
