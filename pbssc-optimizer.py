#!/usr/bin/env python3

# Optimize looks at the data output from pbssc and tests the different parameters to see how many reads they
# include and the error rates.

import os
import argparse
import pickle
import csv
import multiprocessing
from datetime import datetime
from statistics import mean
from time import time
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
pairwise2.MAX_ALIGNMENTS = 100

__author__ = "Theodore (Ted) B. Verhey"
__version__ = "3.0"
__email__ = "verheytb@gmail.com"
__status__ = "Development"

# these are the barcodes that constitute the control bins to be used in this analysis:

# B31 SMRT Cells
barcodes = ("0005_Forward--0018_Reverse", "0006_Forward--0018_Reverse", "0007_Forward--0018_Reverse",
            "0002_Forward--0019_Reverse", "0003_Forward--0019_Reverse")
# JD1 SMRT Cells
# barcodes = ("0001_Forward--0010_Reverse", "0002_For ward--0010_Reverse", "0003_Forward--0010_Reverse")


# parse arguments and display help information
parser = argparse.ArgumentParser(description="Optimize looks at the data output from SSCONSENSUS and tests the " +
                                             "different parameters to plot the resulting dataset size vs error rates.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input_dir",
                    required=True,
                    type=str,
                    help='The directory containing the SSCONSENSUS files (eg. "~/Desktop/Project1")')
parser.add_argument("-c", "--cpus",
                    default=multiprocessing.cpu_count(),
                    type=int,
                    help="Number of CPUs to use")
parser.add_argument("-z", "--nozeros",
                    action="store_true",
                    help="Exclude parameter sets with no reads from the report.")
parser.add_argument("--reference", "-r",
                    required=True,
                    default=None,
                    type=str,
                    help="FASTA file containing a single reference sequence to compare to.")
parser.add_argument("--cassettes", "-s",
                    required=True,
                    default=None,
                    type=str,
                    help="A .pyc file from a cassette_aligner.py run that contains the best cassette alignments.")
args = parser.parse_args()


class Counter(object):
    # a multiprocess-safe counter
    def __init__(self):
        self.val = multiprocessing.Value('i', 0)

    def increment(self, n=1):
        with self.val.get_lock():
            self.val.value += n

    @property
    def value(self):
        return self.val.value


class Alignment(object):
    """
    Corresponds to one way of aligning a read to a reference.
    """

    def __init__(self, gappy_r, gappy_q):
        self.transform = get_transform(gappy_q=gappy_q, gappy_r=gappy_r)
        self.cigar = get_cigar(gappy_r=gappy_r, gappy_q=gappy_q)
        if self.cigar[0][0] == 2:
            self.begin = self.cigar[0][1]
            self.cigar.pop(0)
        else:
            self.begin = 0
        self.snps = count_snps(self.transform)


class Read(object):
    """
    This is a single-stranded consensus read, which may be mapped or unmapped. When mapped, all of the equivalently
    highest-scoring alignments are retained.
    """
    def __init__(self, query, reference):
        self.reference = reference
        self.query = query
        self.alignments = None
        self.is_aligned = False

    def align(self, cassettes):
        self.alignments = align(reference=self.reference, query=self.query)
        self.is_aligned = True

        # choose alignments that are closest to the ensemble of cassettes.
        alignment_distances = [sum(map_distance(aln, cas.alignments) for cas in cassettes) for aln in self.alignments]
        min_distance = min(alignment_distances)
        self.alignments = [self.alignments[x] for x, score in enumerate(alignment_distances) if score == min_distance]


class SscRead(object):
    # a class that carries some metrics for each single-stranded consensus
    def __init__(self, name, barcode, num_passes, coverage,
                 predicted_accuracy, min_confidence, trim_fail, mapping_fail):
        self.name = name
        self.barcode = barcode
        self.numPasses = int(num_passes)
        self.coverage = int(coverage)
        self.predictedAccuracy = float(predicted_accuracy)
        self.minConfidence = float(min_confidence)
        self.trimFail = trim_fail == "True"  # to convert string to boolean
        self.mappingFail = mapping_fail == "True"


# General-purpose functions


def print_message(contents, ontop=False):
    """

    :param contents: a string to print as a message
    :param ontop: if True, does not end with a newline
    :return:
    """
    if not ontop:
        print(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "   " + contents)
    else:
        print(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "   " + contents, end="\r")


def get_cigar(gappy_q, gappy_r):
    """

    :param gappy_q: gapped query sequence
    :param gappy_r: gapped reference sequence
    :return: returns a list of lists, where each sublist is [operation, length]
    """
    assert len(gappy_q) == len(gappy_r)
    cigar = []
    for q, r in zip(gappy_q, gappy_r):
        if q == "-":  # deletion
            if cigar and cigar[-1][0] == 2:
                cigar[-1][1] += 1
            else:
                cigar.append([2, 1])
        elif r == "-":  # insertion
            if cigar and cigar[-1][0] == 1:
                cigar[-1][1] += 1
            else:
                cigar.append([1, 1])
        else:
            if cigar and cigar[-1][0] == 0:
                cigar[-1][1] += 1
            else:
                cigar.append([0, 1])
    return cigar


def get_transform(gappy_q, gappy_r):
    """

    :return: tuple of tuples
    """
    transform = []
    rgaps = 0  # counts number of reference gaps
    for x, r in enumerate(str(gappy_r)):
        x_ref = x-rgaps
        q = str(gappy_q)[x]
        if q == r:  # match
            continue
        elif q == "-":  # deletion
            transform.append([x_ref, "D"])
        elif r == "-":  # insertion
            if transform and transform[-1][1] == "I" and transform[-1][0] == x_ref:
                transform[-1][2] += q
            else:
                transform.append([x_ref, "I", q])
            rgaps += 1
        else:  # mismatch
            transform.append([x_ref, "S", q])
    return tuple(tuple(x) for x in transform)  # convert to tuple of tuples


def transform(reference, transformer):
    query = list(reference.seq)
    transformer = sorted(transformer, key=lambda x: x[0])
    offset = 0
    for op in transformer:
        pos = op[0] + offset
        if op[1] == "S":
            query[pos] = op[2]
        elif op[1] == "I":
            query[pos:pos] = list(op[2])
            offset += len(op[2])
        elif op[1] == "D":
            query.pop(pos)
            offset -= 1
    return "".join(query)


def count_snps(transform, start=None, stop=None):
    """
    Counts the number of SNPs represented in a transform within the [start, stop] interval
    :param transform:
    :param start:
    :param stop:
    :return:
    """
    snp_count = 0
    for op in transform:
        if start and op[0] < start:
            continue
        if stop and op[0] > stop:
            continue
        if op[1] in "SD":
            snp_count += 1
        elif op[1] in "I":
            snp_count += len(op[2])
    return snp_count


def map_distance(aln1, aln2):
    """
    Measures the distance of any two reads from each based on a mapping to a common reference
    :param aln1: A mapping of a read to a reference
    :param aln2: A mapping of another read to the same reference
    :return: an integer score corresponding to the number of different bases that differ between the two mappings.
    """
    score = 0
    # substitutions and deletions are per base, so we just count the number of different entries in the transforms
    aln1_sd = {x for x in aln1.transform if x[1] in ("S", "D")}
    aln2_sd = {x for x in aln2.transform if x[1] in ("S", "D")}
    score += len(aln1_sd ^ aln2_sd)

    # insertions that are non-identical must be compared at the subsequence level in order to count the number of
    # differences
    aln1_i = {x for x in aln1.transform if x[1] == "I"}
    aln2_i = {x for x in aln2.transform if x[1] == "I"}
    inserts = aln1_i ^ aln2_i
    while inserts:
        position1, _, inseq1 = inserts.pop()
        others = {x for x in inserts if x[0] == position1}
        if others:
            assert len(others) == 1
            other = others.pop()
            inserts.discard(other)
            position2, _, inseq2 = other
            if inseq1.endswith(inseq2) or inseq1.startswith(inseq2):
                score += len(inseq1) - len(inseq2)
            elif inseq2.endswith(inseq1) or inseq2.startswith(inseq1):
                score += len(inseq2) - len(inseq1)
            else:
                score += max(len(inseq1), len(inseq2))
        else:
            score += len(inseq1)
    return score


def align(reference, query):
    """
    do a pairwise alignment of the query to the reference, outputting up to 10000 of the highest-scoring alignments.
    :param reference: a SeqRecord object containing the reference
    :param query: a SeqRecord object containing the query sequence
    :return: a list of up to 10000 Alignment objects
    """
    alns = pairwise2.align.localms(reference.seq, query.seq, 1, -1, -2, -1)
    alignments = []
    for aln in alns:
        al1, al2, score, begin, end = aln
        alignments.append(Alignment(gappy_r=al1, gappy_q=al2))
    return alignments


def align_worker(read):
    read.align(cassettes)
    return read


if __name__ == "__main__":

    # makes temporary directory and report file
    tempdir = "/tmp/optimize_" + datetime.now().strftime("%Y%m%d-%H%M%S")
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)
    reportpath = os.path.join(args.input_dir, "optimize_report_" + datetime.now().strftime("%Y%m%d-%H%M%S") + ".csv")

    # defines some global objects
    startTime = time()

    # Load the read statistics from CSV to memory for the whole dataset
    print_message("Importing data files.", ontop=False)
    with open(os.path.join(args.input_dir, "Report.csv")) as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        next(csvreader)  # skip the header
        read_stats = [SscRead(*row) for row in csvreader]

    # extracts the reference sequences
    with open(args.reference, "r") as h:
        reference = SeqIO.read(h, "fasta")

    # Loads the control read sequences
    control_reads = []
    for x in barcodes:
        with open(os.path.join(args.input_dir, "Demultiplexed/" + x + ".fasta"), "r") as h:
            for record in SeqIO.parse(h, format="fasta"):
                control_reads.append(Read(query=record, reference=reference))

    print_message("There are %d reads to align." % len(control_reads))

    # imports cassettes
    with open(args.cassettes, "rb") as h:
        prealigned_cassettes = pickle.load(h)
    cassettes = []
    for name, seq, aln_obj in prealigned_cassettes:
        cassette = Read(query=SeqRecord(Seq(seq), id=name), reference=reference)
        cassette.alignments = aln_obj
        cassette.is_aligned = True
        cassettes.append(cassette)

    print_message("Aligning the reads")
    with multiprocessing.Pool(args.cpus) as P:
        control_reads = P.map(align_worker, control_reads)
    # starts counting
    print_message("Comparing parameter sets")
    with open(reportpath, 'a', newline='') as rhandle:
        writer = csv.writer(rhandle)
        writer.writerow(["Coverage", "Predicted_Accuracy", "Min_Confidence", "Dataset Size (Reads)", "Errors (bp per read)"])
        for coverage in range(0, 41, 1):
            for predicted_accuracy in (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.80, 0.825, 0.85, 0.875, 0.9, 0.925,
                                       0.95, 0.975, 0.99, 0.995, 0.9975, 0.999, 0.9995, 0.99975, 0.9999, 0.99999):
                for min_confidence in (0, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 0.9125,
                                       0.925, 0.95, 0.975, 0.99, 0.9925, 0.995, 0.9975, 0.999):
                    selected_read_stats = [x for x in read_stats if x.coverage >= coverage and
                                           x.predictedAccuracy >= predicted_accuracy and
                                           x.minConfidence >= min_confidence]  # a subset of reads from the whole dataset
                    selected_control_reads = [read for read in control_reads if read.query.id in
                                              [x.name for x in selected_read_stats]]
                    error_rate = mean(mean(aln.snps for aln in read.alignments) for read in selected_control_reads)
                    writer.writerow([coverage, predicted_accuracy, min_confidence, len(selected_read_stats), error_rate])

    elapsedTime = datetime.fromtimestamp(time()) - datetime.fromtimestamp(startTime)
    print_message("\nFinished! Elapsed time is %s" % str(elapsedTime), ontop=False)
