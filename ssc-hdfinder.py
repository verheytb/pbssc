#!/usr/bin/env python3

# ssc-hdfinder takes single-stranded reads (specifically, those in the format outputted by pbssc) and identifies
# the number that have different sequences in either strand.

__author__ = 'Ted Verhey'
__email__ = 'verheytb@gmail.com'

import argparse, os, sys, csv, glob
from datetime import datetime
import swalign
from pyfaidx import Fasta
from Bio import SeqIO

scoring = swalign.NucleotideScoringMatrix(2, -1)

def get_opposite_strand(readname):
    if readname.endswith("/0/ssc"):
        return readname[:readname.rfind("/0/ssc")] + "/1/ssc"
    elif readname.endswith("/1/ssc"):
        return readname[:readname.rfind("/1/ssc")] + "/0/ssc"
    else:
        sys.exit("Error: Report.csv is not properly formatted.")


def is_similar(seq1, seq2, window):
    """

    :param seq1:
    :param seq2:
    :return: True if there are no differences or if the differences are 20bp-separated single indels.
    """
    if seq1 == seq2:
        return True
    sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...
    al = sw.align(seq1, seq2)
    seq1g = seq1
    seq2g = seq2
    assert al.q_pos == al.r_pos
    position = 0
    for c in al.cigar:
        if c[1] == "M":
            position += c[0]
        elif c[1] == "I":
            seq1g = seq1g[0:position] + "-"*c[0] + seq1g[position:]
        elif c[1] == "D":
            seq2g = seq2g[0:position] + "-"*c[0] + seq2g[position:]
    assert len(seq1g) == len(seq2g)
    for x in range(len(seq1g)):
        if seq1g[x] != seq2g[x]:
            if seq1g[x] in ("A", "C", "T", "G") and seq2g[x] in ("A", "C", "T", "G"):  # substitution
                return False
            elif seq1g[x] == "-" or seq2g[x] == "-":  # indels
                right_side = max(x-window, 0)
                left_side = min(x+window, len(seq1g))
                flankseqs = [seq1g[right_side:x], seq1g[x+1:left_side],
                             seq2g[right_side:x], seq2g[x+1:left_side]]
                if any([("-" in s) for s in flankseqs]):
                    return False
            else:
                sys.exit("Error: DNA sequences must be in capital letters.")
    return True

# Displays usage information and parses arguments.
parser = argparse.ArgumentParser(description="SSC-HDFINDER takes single-stranded reads (specifically, those in the " +
                                             "format outputted by pbssc) and identifies the number that have " +
                                             "different sequences in either strand.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-j", "--jobDir",
                    required=True,
                    type=str,
                    help='The directory containing output from pbssc. Must contain a "Demultiplexed" folder ' +
                         'and a Report.csv file. (eg. "~/Desktop/Project1")')
parser.add_argument("-r", "--reference",
                    required=False,
                    type=str,
                    help="A FASTA file containing the reference sequence by which to predict the expected level of " +
                         "heteroduplexes. It must be trimmed to where the sequences are trimmed.")
parser.add_argument("-s", "--ignore_single_indels",
                    required=False,
                    action="store_true",
                    help="If specified, heteroduplexes are not counted if both strands differ from each other by " +
                         "only a single indel.")
parser.add_argument("-w", "--window_length",
                    required=False,
                    default=20,
                    type=int,
                    help="Ignoring equencing errors (single indels) that are separated from each other by greater " +
                         "than this distance.")
args = parser.parse_args()

# Checks input files
if not os.path.exists(args.jobDir):
    sys.exit("Error: specified job directory does not exist.")
if not os.path.exists(os.path.join(args.jobDir, "Report.csv")):
    sys.exit("Error: Report.csv file could not be found in the specified job directory.")
if len(glob.glob(os.path.join(args.jobDir, "Demultiplexed/*.fasta"))) == 0:
    sys.exit('Error: No FASTA files could be found in the "Demultiplexed" directory.')

# Reads the CSV to memory for quick access.
with open(os.path.join(args.jobDir, "Report.csv"), "r") as r:
    csvReader = csv.reader(r)
    next(csvReader)
    reportDict = {}
    for row in csvReader:
        reportDict[row[0]] = row[1:]

# Creates a new dictionary mapping the barcodes to good reads. Discards low quality reads and those with only one
# strand.
barcodeForPrefix = {}
for readName, stats in reportDict.items():
    prefix = readName[:-6]
    ignoreflag = False
    if not prefix in barcodeForPrefix:
        if not get_opposite_strand(readName) in reportDict:  # Checks that the other strand exists
            continue
        if not reportDict[readName][0] == reportDict[get_opposite_strand(readName)][
            0]:  # Ensure strands have same barcodes
            continue
        for strand in (readName, get_opposite_strand(readName)):
            # if reportDict[strand][0] == "--":  # Check that both strands have non-null barcodes
            #    ignoreflag = True
            # if not float(reportDict[strand][3]) >= 0:  # Checks that they both exceed 0 estimated accuracy
            #    ignoreflag = True
            if not float(reportDict[strand][4]) >= 0.9:  # Checks that they both exceed 0.9 minimum basecall confidence
                ignoreflag = True
            if not reportDict[strand][5] == "False" and not reportDict[strand][6] == "False":
                ignoreflag = True
        if not ignoreflag:
            barcodeForPrefix[prefix] = reportDict[readName][0]

fastas = {}
hdcount = {}
numsequences = {}

for prefix, barcode in barcodeForPrefix.items():
    if barcode not in fastas:
        fastas[barcode] = Fasta(os.path.join(args.jobDir, "Demultiplexed/" + barcode + ".fasta"), as_raw=True)
        hdcount[barcode] = 0
        numsequences[barcode] = 0
    numsequences[barcode] += 1
    fstrand = fastas[barcode][prefix + "/0/ssc"][:]
    rstrand = fastas[barcode][prefix + "/1/ssc"][:]
    if args.ignore_single_indels:
        if not is_similar(fstrand, rstrand, args.window_length):
            hdcount[barcode] += 1
    else:
        if not fstrand == rstrand:
            hdcount[barcode] += 1

totalsequences = {}
differentsequences = {}

if args.reference:
    with open(args.reference, 'r') as r:
        reference = str(SeqIO.read(r, 'fasta').seq)
    for barcode, file in fastas.items():
        if barcode not in totalsequences:
            totalsequences[barcode] = 0
            differentsequences[barcode] = 0
        for entry in file:
            totalsequences[barcode] += 1
            if args.ignore_single_indels:
                if not is_similar(entry[:], reference, args.window_length):
                    differentsequences[barcode] += 1
            else:
                if entry[:] != reference:
                    differentsequences[barcode] += 1

csvOut = open(os.path.join(args.jobDir, "ssc-hdfinder_report_" + datetime.now().strftime("%Y.%m.%d-%H.%M.%S") + ".csv"),
              "w")
csvOut.write(
    "Barcode,Number_of_Sequences,Heteroduplexes,HD_Freq,Total_Strands,Switched_Strands,Switched_Strands_Freq\n")
for barcode in fastas:
    if args.reference:
        printarray = (barcode, str(numsequences[barcode]), str(hdcount[barcode]),
                      str(hdcount[barcode] / numsequences[barcode]), str(totalsequences[barcode]),
                      str(differentsequences[barcode]), str(differentsequences[barcode] / totalsequences[barcode]))
    else:
        printarray = (barcode, str(numsequences[barcode]), str(hdcount[barcode]),
                      str(hdcount[barcode] / numsequences[barcode]), '', '', '')
    csvOut.write("%s,%s,%s,%s,%s,%s,%s\n" % printarray)
