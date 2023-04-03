#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Łukasz Szafron
@email: lukszafron@gmail.com
"""
appname = "ENA-submitter"
app_version = "1.0"
webin_cli_path = "/programs/webin-cli"

import pandas as pd
import time, datetime, subprocess, re, os, sys, getopt, multiprocessing

start_time = time.time()
def usage():
    print(
        "\nWelcome to the " + appname + " app.\n\n"
        "The following options are available:\n\n"
        "\t-u, --user: Webin user name\n"
        "\t-p --password: Webin password\n"
        "\t-S --study_id: Webin study ID\n"
        "\t-n --study_name: Webin study NAME\n"
        "\t-I --instrument: NGS instrument name, e.g., 'Illumina NovaSeq 6000', 'Illumina iSeq 100'\n"
        "\t-i, --insert_size: The number of sequenced nucleotide in two reads forming a pair\n"
        "\t-s, --source:, NGS library source, e.g., 'GENOMIC', 'GENOMIC SINGLE CELL', 'TRANSCRIPTOMIC', 'TRANSCRIPTOMIC SINGLE CELL'\n"
        "\t-c, --library_selection: Method of the NGS library selection (default: RANDOM)\n"
        "\t-t, --library_strategy: NGS library type, e.g., 'WGS' (genome), 'WXS' (exome), 'RNA-Seq', 'miRNA-Seq'\n"
        "\t-f, --fastq_dir: defines the path to the directory containing FASTQ files.\n"
        "\t-a, --Webin_accessions: defines the path to a tab-delimited txt file with sample accession numbers (samples aliases should be the same as sample prefixes in the corresponding FASTQ files).\n"
        "\t-T, --threads: number of CPU threads to use (default: 1)\n"
        "\t-h, --help: prints this help message.\n"
        "\t-v, --version: prints the version of this program.\n"
        )
# The next two lines are for debugging purposes only and should be commented in the final program.
# option_list = ["-u", "********", "-p", "********", "-S", "PRJEB61011", "-n", "CRNDE_RNA-seq", "-I", "Illumina NovaSeq 6000", "-i", "200", "-s", "TRANSCRIPTOMIC", "-t", "RNA-Seq","-a", "/workspace/lukasz/NGS-all-in-one/RUNS/CRNDEAB1/FASTQ/Webin-accessions-2023-03-29T12 59 28.917+01 00.txt", "-f", "/workspace/lukasz/NGS-all-in-one/RUNS/CRNDEAB1/FASTQ", "-T", "5"]
# opts, args = getopt.getopt(option_list, "u:p:S:n:I:i:s:c:t:f:a:T:hv", ["user=","password=","study_id=","study_name=","instrument=","insert_size=","source=","library_selection=","library_strategy=","fastq_dir=","Webin_accessions=","threads=","help","version"])

opts, args = getopt.getopt(sys.argv[1:], "u:p:S:n:I:i:s:c:t:f:a:T:hv", ["user=","password=","study_id=","study_name=","instrument=","insert_size=","source=","library_selection=","library_strategy=","fastq_dir=","Webin_accessions=","threads=","help","version"])

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
            elif o in ("-f", "--fastq_dir"):
                fastq_dir = a
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
    print("\n"+str(err),"red") # will print something like "option -a not recognized"
    usage()
    sys.exit(2)
    
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
    raise(Exception("The study id is missing."))
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
    raise(Exception("The insert size is missing."))
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
    fastq_dir
except:
    raise(Exception("The path to the FASTQ-containing directory is missing."))

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
        self.matching_fastqs = [fastq for fastq in os.listdir(fastq_dir) if re.search(string=fastq, pattern="^"+alias[1].replace(".", "\.")+".*\.fastq\.gz$")]
    def fastq_exists(self):
        if(len(self.matching_fastqs) != 2):
            raise Exception("There is an incorrect number of merged FASTQ files ({1}) for the sample: {0}.\n".format(self.alias[1], len(self.matching_fastqs)))
    def run_tests(self):
        # print("Testing if exactly two merged FASTQ files exist for the analyzed sample...\n", flush=True)
        self.fastq_exists()
    def manifest(self):
        # print("Generating the MANIFEST file...\n", flush=True)
        self.manifest_dict = {"SAMPLE":self.alias[0], "STUDY":study_id, "NAME":'_'.join([study_name, self.alias[1]]), "INSTRUMENT":instrument, "INSERT_SIZE":insert_size, "LIBRARY_SOURCE":source, "LIBRARY_SELECTION":lib_selection, "LIBRARY_STRATEGY":lib_strategy}
        self.items = [v for v in self.manifest_dict.items()]
        self.tsv = []
        for i,v in enumerate(self.items):
            self.tsv.append('\t'.join(self.items[i]))
        self.tsv = self.tsv + ['\t'.join(["FASTQ", fastq]) for fastq in self.matching_fastqs]
        self.file_alias = os.path.join(fastq_dir, self.alias[1]+"_Manifest.txt")
        with open(file = self.file_alias, mode = "wt") as out:
            out.write('\n'.join(self.tsv))
    def submission(self):
#        print("Submitting the FASTQ files and metadata to the European Nucleotide Archive (ENA) database...\n", flush=True)
        self.res = subprocess.run([webin_cli_path, "-context reads", "-inputDir", fastq_dir, "-outputDir", fastq_dir, "-manifest", self.file_alias, "-userName", user, "-password", password, "-submit"], text = True)
        if self.res.returncode == 0:
            print("{} sample submission is done.\n".format(self.alias[1]), flush=True)
        else:
            print("{} SAMPLE SUBMISSION HAS FAILED!!!\n".format(self.alias[1]), flush=True)

def exec_func(alias):
    sample = SAMPLE(alias)
    sample.run_tests()
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

print("--- The sample submission process took %s hh:mm:ss. ---\n" % str(datetime.timedelta(seconds = round(time.time() - start_time, 0))))
