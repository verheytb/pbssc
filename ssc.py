#!/usr/bin/env python

# SSCONSENSUS computes single-stranded consensus reads from an aligned_reads.cmp.h5 file that has been filtered for
# barcode score ratios using pbbarcode. It outputs to FASTA or FASTQ formats and can also output a CSV file with
# detailed statistics for each read. SSCONSENSUS must be called from the SMRT Shell.

from __future__ import print_function

__author__ = 'Ted Verhey'
__email__ = 'verheytb@gmail.com'

import os
import sys
import argparse
import multiprocessing
from datetime import datetime
from time import time, sleep
import numpy as np
import ConsensusCore as cc
from GenomicConsensus.quiver.model import loadQuiverConfig
from GenomicConsensus.quiver.utils import refineConsensus, refineDinucleotideRepeats,\
    consensusConfidence, QuiverConsensus
from pbcore.io.align.CmpH5IO import CmpH5Reader
from pbcore.io.FastaIO import IndexedFastaReader, FastaWriter
from pbcore.io.FastqIO import FastqWriter


def consensusForAlignments(refWindow, refSequence, alns, quiverConfig):
    """
    This function overrides the one in the GenomicConsensus module in order to start using the reference instead
    of a POA consensus.
    :param refWindow:
    :param refSequence:
    :param alns:
    :param quiverConfig:
    :return:
    """
    _, refStart, refEnd = refWindow
    # No POA --- just use the reference as a starting point.
    poaCss = refSequence
    # Extract reads into ConsensusCore-compatible objects, and map them into the
    # coordinates relative to the POA consensus
    mappedReads = [ quiverConfig.extractMappedRead(aln, refStart) for aln in alns ]
    # Load the mapped reads into the mutation scorer, and iterate
    # until convergence.
    configTbl = quiverConfig.ccQuiverConfigTbl
    mms = cc.SparseSseQvMultiReadMutationScorer(configTbl, poaCss)
    for mr in mappedReads:
        mms.AddRead(mr)
    # Iterate until convergence
    _, quiverConverged = refineConsensus(mms, quiverConfig)
    if quiverConverged:
        if quiverConfig.refineDinucleotideRepeats:
            refineDinucleotideRepeats(mms)
        quiverCss = mms.Template()
        if quiverConfig.computeConfidence:
            confidence = consensusConfidence(mms)
        else:
            confidence = np.zeros(shape=len(quiverCss), dtype=int)
        return QuiverConsensus(refWindow, quiverCss, confidence, mms)
    else:
        return QuiverConsensus.noCallConsensus(quiverConfig.noEvidenceConsensus, refWindow, refSequence)


# parse arguments and display help information
parser = argparse.ArgumentParser(description="SSCONSENSUS computes single-stranded consensus reads from an " +
                                             "aligned_reads.cmp.h5 file that has been filtered for barcode score " +
                                             "ratios using pbbarcode. It outputs to FASTA or FASTQ formats and can " +
                                             "also output a CSV file with detailed statistics for each read. " +
                                             "SSCONSENSUS must be called from the SMRT Shell.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-j", "--jobDir", required=True, type=str,
                    help='The directory containing the SMRT Analysis job for processing ' +
                         '(eg. "~/Programs/smrtanalysis/userdata/jobs/016/016450/")')
parser.add_argument("-r", "--reference", required=True, type=str,
                    help='The reference FASTA file used by SMRT Analysis(eg. ' +
                         '"~/Programs/smrtanalysis/userdata/references/JD1_vlsE_1/sequence/JD1_vlsE_1.fasta")')
parser.add_argument("-o", "--outDir", required=True, type=str,
                    help='The directory for output files (eg. "~/Desktop/Project1")')
parser.add_argument("-a", "--minAvgConfidence", type=float, default=None,
                    help="Ignores strand-specific consensuses below an average confidence value (0,1).")
parser.add_argument("-m", "--minConfidence", type=float, default=None,
                    help="Ignores strand-specific consensuses where any base is below a confidence value (0,1).")
parser.add_argument("-p", "--minCoverage", type=int, default=None,
                    help="Ignores strand-specific consensuses that have fewer than the specified number of passes, " +
                         "including partial passes.")
parser.add_argument("-t", "--trim", type=str, default=None,
                    help="Trims all sequences to the specified tuple of sequences (eg. ACAGCTG, CGGCGAAT). These " +
                         "sequences must be in the same strand as the reference.")
parser.add_argument("-q", "--fastq", action='store_true', help="Outputs FASTQ files instead of FASTA files.")
parser.add_argument("-c", "--cpus", default=multiprocessing.cpu_count() - 1, type=int, help="Number of CPUs to use")
args = parser.parse_args()


def printmessage(contents, ontop):
    if not ontop:
        print(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "   " + contents)
    else:
        print(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "   " + contents, end="\r")


def intersect(i1, i2):
    s1, e1 = i1
    s2, e2 = i2
    return not((s1 >= e2) or (s2 >= e1))


def hull(i1, i2):
    s1, e1 = i1
    s2, e2 = i2
    return (min(s1, s2), max(e1, e2))


def hullMany(intervals):
    h = (float("inf"),float("-inf"))
    for v in intervals:
        h = hull(h, v)
    return h


def within(needle, haystack):
    ns, ne = needle
    hs, he = haystack
    return hs <= ns <= ne <= he


def ilen(interval):
    return interval[1]-interval[0]


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


class ConsensusSequence(object):
    # a class that carries a consensus sequence and some metrics
    def __init__(self, name, seq, qual, numPasses, coverage, predictedAccuracy, barcode):
        self.name = name
        self.seq = seq
        self.qual = qual
        self.numPasses = numPasses
        self.coverage = coverage
        self.predictedAccuracy = predictedAccuracy
        self.minConfidence = 0
        self.barcode = barcode
        self.minNumPassesFail = False
        self.mappingFail = False
        self.trimFail = False
        self.minCoverageFail = False
        self.minAvgConfidenceFail = False
        self.minConfidenceFail = False


def mappingsAreConcordant(templateIntervals):
    # Are the mappings concordant?
    FUDGE = 30
    s, e = max(templateIntervals, key=ilen)
    s -= FUDGE
    e += FUDGE
    return all(within(i, (s, e)) for i in templateIntervals)


def checkMapping(alnSubreads):
    # Require unique mapping for now.
    if not all(sr.MapQV == 254 for sr in alnSubreads): return False
    if not len(set(sr.referenceId for sr in alnSubreads)) == 1: return False
    templateIntervals = [(a.tStart, a.tEnd) for a in alnSubreads]
    if not mappingsAreConcordant(templateIntervals): return False
    return True


def unphred(q):
    return 10 ** (-q / 10)


def phred(p):
    return -10 * np.log10(p)


def estimateAccuracy(confidence):
    errorProbs = unphred(np.array(confidence, dtype=float))
    errorRate = np.mean(errorProbs)
    errorRate = (errorRate * len(confidence) + 1.) / (len(confidence) + 1)
    return 1. - errorRate


def trim(cssObj, leftSeq, rightSeq):
    # trims a cssObj to the specified sequences
    if all(s in cssObj.seq for s in (leftSeq, rightSeq)):
        a = cssObj.seq.find(leftSeq)
        b = cssObj.seq.rfind(rightSeq) + len(rightSeq)
        cssObj.seq = cssObj.seq[a:b]
        cssObj.qual = cssObj.qual[a:b]
    else:
        cssObj.trimFail = True
    return cssObj


# The worker function
def workerProcess(inQueue, refFile, quiverConfig):
    countFromReset = 0
    with CmpH5Reader(os.path.join(args.jobDir, "data/aligned_reads.cmp.h5")) as d:
        while not inQueue.empty() and countFromReset < 10000:
            referenceSeq = IndexedFastaReader(refFile)[0].sequence
            while not inQueue.empty():
                countFromReset += 1
                movieID, holeNumber, rcrefstrand = inQueue.get()
                alns = d[((d.MovieID == movieID) & (d.HoleNumber == holeNumber) & (d.RCRefStrand == rcrefstrand))]
                cssName = "/".join(alns[0].readName.split("/")[:-1]) + "/" + str(alns[0].RCRefStrand) + "/ssc"
                cssObj = ConsensusSequence(cssName, "", "", len(alns), 0, 0, alns[0].barcodeName)
                if not cssObj.numPasses >= args.minCoverage:
                    cssObj.minNumPassesFail = True
                if not checkMapping(alns):
                    cssObj.mappingFail = True
                if not cssObj.minNumPassesFail and not cssObj.mappingFail:
                    refId = alns[0].referenceId
                    v = hullMany([(a.tStart, a.tEnd) for a in alns])
                    window = (refId, v[0], v[1])
                    windowLen = v[1] - v[0]
                    refSeqInWindow = referenceSeq[v[0]:v[1]]
                    css = consensusForAlignments(window, refSeqInWindow, alns, quiverConfig)
                    cssObj.seq = css.sequence
                    cssObj.qual = css.confidence
                    cssObj.coverage = sum((a.referenceSpan > 0.8 * windowLen) for a in alns)
                    if not cssObj.coverage >= args.minCoverage:
                        cssObj.minCoverageFail = True
                    cssObj.predictedAccuracy = estimateAccuracy(css.confidence)
                    if not cssObj.predictedAccuracy >= args.minAvgConfidence:
                        cssObj.minAvgConfidenceFail = True
                    if args.trim:
                        cssObj = trim(cssObj, lseq, rseq)
                    cssObj.minConfidence = 1 - unphred(np.amin(np.array(cssObj.qual, dtype=float)))
                    if not cssObj.minConfidence >= args.minConfidence:
                        cssObj.minConfidenceFail = True
                resultQueue.put(cssObj)
                counter.increment()


def writerProcess(outDir):
    # makes output directories
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    fastOutDir = os.path.join(outDir, "Demultiplexed/")
    if not os.path.exists(fastOutDir):
        os.makedirs(fastOutDir)

    # opens files
    csvOut = open(os.path.join(outDir, "Report.csv"), "w")
    csvOut.write("Name,Barcode,NumPasses,Coverage,AvgConfidence,MinConfidence,TrimFail,MappingFail\n")
    writers = {}
    while counter.value < totalNumber:
        result = resultQueue.get()
        csvOut.write("%s,%s,%d,%d,%0.6f,%0.6f,%s,%s\n" % (
            result.name, result.barcode, result.numPasses, result.coverage, result.predictedAccuracy,
            result.minConfidence, result.trimFail, result.mappingFail))
        if not result.barcode in writers:
            if args.fastq:
                writers[result.barcode] = FastqWriter(os.path.join(fastOutDir, result.barcode + ".fastq"))
            else:
                writers[result.barcode] = FastaWriter(os.path.join(fastOutDir, result.barcode + ".fasta"))
        if not any((result.minNumPassesFail, result.mappingFail, result.trimFail, result.minCoverageFail,
                    result.minAvgConfidenceFail, result.minConfidenceFail)):
            if args.fastq:
                writers[result.barcode].writeRecord(result.name, result.seq, result.qual)
            else:
                writers[result.barcode].writeRecord(result.name, result.seq)


# defines some global objects
taskQueue = multiprocessing.Queue()
resultQueue = multiprocessing.Queue()
counter = Counter()
startTime = time()
quiverErrorModel = loadQuiverConfig("P6-C4.AllQVsMergingByChannelModel")

if args.trim:
    dupleSeq = args.trim.upper()
    if not all([x in ("A", "C", "G", "T", ",") for x in dupleSeq]) or not dupleSeq.count(",") == 1:
        sys.exit("-t requires a duple of sequences separated by a comma. (eg. ACTAGGA,CTACGAG)")
    lseq = dupleSeq[:args.trim.find(",")]
    rseq = args.trim[args.trim.find(",") + 1:]

# get the list of reads to extract, and populate the task queue
printmessage("Importing data files", ontop=False)
with CmpH5Reader(os.path.join(args.jobDir, "data/aligned_reads.cmp.h5")) as c:
    c.attach(os.path.join(args.jobDir, "input.fofn"))
    readSet = set(zip(c.MovieID, c.HoleNumber, c.RCRefStrand))
totalNumber = len(readSet)
for i in readSet:
    taskQueue.put(i)

# starts the processes
processList = [multiprocessing.Process(target=workerProcess, args=(taskQueue, args.reference, quiverErrorModel)) for i
               in range(args.cpus)]
for i in processList:
    i.start()
printmessage("%d processes started" % args.cpus, ontop=False)

# starts a process to write the results queue to disk.
writerProcess = multiprocessing.Process(target=writerProcess, args=(args.outDir,))
writerProcess.start()

# every second, updates terminal and checks for dead processes
while counter.value < totalNumber:
    printmessage("Generating single-stranded consensuses: %d of %d (Write backlog: %d)" % (
        counter.value, totalNumber, resultQueue.qsize()), ontop=True)
    sys.stdout.flush()
    for x, i in enumerate(processList):
        if not i.is_alive():
            i.join()
            processList[x] = None
            processList[x] = multiprocessing.Process(target=workerProcess,
                                                     args=(taskQueue, args.reference, quiverErrorModel))
            processList[x].start()
    sleep(0.1)

# waits for processes to finish
for i in processList:
    i.join()
writerProcess.join()

elapsedTime = datetime.fromtimestamp(time()) - datetime.fromtimestamp(startTime)
printmessage("\nFinished! Elapsed time is %s" % str(elapsedTime), ontop=False)
