#!/usr/bin/env python3

# fastq-hdfinder gets the number of heteroduplexes from a set of FASTQ CCS Reads by searching for reads that have
# low -confidence basecalls.

__author__ = 'Ted Verhey'
__email__ = 'verheytb@gmail.com'

import argparse
import csv
import glob
import os
import sys
from datetime import datetime

from Bio import SeqIO

# parses the command-line arguments and displays usage information
parser = argparse.ArgumentParser(
    description="fastq-hdfinder gets the number of heteroduplexes from a set of FASTQ CCS" +
                " reads by searching for those that have low-confidence basecalls.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--input_dir", "-i",
                    required=True,
                    type=str,
                    help="Directory containing fasta files for alignment.")
parser.add_argument("--qvcutoff", "-q",
                    type=int,
                    default=5,
                    help="An integer QV score in the interval [0,93], including and below which to call a "
                         "heteroduplex.")
args = parser.parse_args()
inputdir = os.path.normpath(args.input_dir)

# imports file names from input directory, ignoring empty files
if not os.path.exists(inputdir):
    sys.exit("Error: Input directory doesn't exist")
filelist = [fn for fn in glob.iglob(os.path.join(inputdir, "*.fastq")) if os.stat(fn).st_size > 0]
filelist.sort()
if len(filelist) == 0:
    sys.exit("Error: No non-empty FASTQ files in input directory.")

# opens output report file
reportpath = os.path.join(inputdir, "fastq-hdfinder_report_" + datetime.now().strftime("%Y.%m.%d-%H.%M.%S") + ".csv")
with open(reportpath, 'a', newline='') as outhandle:
    writer = csv.writer(outhandle, delimiter=',')
    writer.writerow(["Phred cutoff:", str(args.qcutoff)])
    writer.writerow([""])
    writer.writerow(["Filename", "Total Reads", "Heteroduplexes", "Heteroduplex reads per read"])
    for fn in filelist:
        with open(fn, "r") as fn_handle:
            num_records = len(SeqIO.index(fn, "fastq"))
            num_heteroduplexes = 0
            for record in SeqIO.parse(fn_handle, "fastq"):
                if min(record.letter_annotations["phred_quality"]) <= args.qcutoff:
                    num_heteroduplexes += 1
        writer.writerow([os.path.basename(fn), num_records, num_heteroduplexes, num_heteroduplexes / num_records])
