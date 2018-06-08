[![DOI](https://zenodo.org/badge/DOI/10.1111/mmi.13873.svg)](https://doi.org/10.1111/mmi.13873)

# pbssc
Generates single-stranded CCS reads from PacBio .bax.h5 files of multiplexed amplicon sequences from an RSII instrument.

## Introduction
PacBio's CCS method in SMRT Analysis 2.3 cannot distinguish between the strands of heteroduplex amplicons, so pbssc is a workaround method that is compatible with SMRT Analysis 2.3 and the .bax.h5 files generated from RSII instruments. 

## Dependencies

`pbssc.py` is compatible with Python 2.7, and requires the following modules to run:
- [ConsensusCore](https://github.com/PacificBiosciences/ConsensusCore)
- [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus)
- [pbcore](https://github.com/PacificBiosciences/pbcore)

In order to do the first few steps of the analysis, you will also need [SMRT Analysis 2.3](http://www.pacb.com/support/software-downloads/).

## Instructions
1. Using SMRT Analysis, set up a job with the RS_Resequencing_Barcode.1 method to align subreads to the reference sequence. Use an appropriate reference sequence for the target amplicon, and ensure that the barcode FASTA file meets the [specifications](http://www.pacb.com/wp-content/uploads/2015/09/Shared-Protocol-PacBio-Barcodes-for-SMRT-Sequencing.pdf) and is appropriate for your sample. Since the reads must align to the reference, pbssc only generates single-stranded reads for variants of a known amplicon sequence.

2. Take note of the job number and directory; it should be in your `userdata` directory that was specified in the SMRT Analysis installation. An example job directory would be `/opt/smrtanalysis/userdata/jobs/016/016531`. In a terminal, assign that directory to a variable, as well as the location of your SMRT Analysis installation.
  ```sh
  $ JOBDIR="/opt/smrtanalysis/userdata/jobs/016/016531"
  $ SMRT_HOME="/opt/smrtanalysis"
  ```

3. Run `pbbarcode labelAlignments` from the SMRT Shell to filter reads for those that were high-scoring during demultiplexing. Failure to do this step will result in incorrect barcode sorting.
  ```sh
  $ $SMRT_HOME/smrtcmds/bin/smrtshell
  $ pbbarcode labelAlignments --minScoreRatio 1.1 $JOBDIR/data/barcode.fofn $JOBDIR/data/aligned_reads.cmp.h5
  $ exit
  ```

4. Ensure that you have run the exit command to leave the SMRT Shell. Then run `pbssc.py`, using the reference originally used in the SMRT Analysis Resequencing Job. The reference must have an associated FASTA index file (references stored in SMRT Portal's 'userdata/references/' directory are already indexed).
  ```sh
  $ python pbssc.py -j $JOBDIR -r reference_amplicon.fasta -o outdir -t ATCTTCGATCGA,TGTAACTGAAGA
  ```

## Detecting Heteroduplexes
Heteroduplexes can be measured in demultiplexed samples using two methods, implemented in the following included scripts:

### FASTQ method
Using `fastq-hdfinder.py`, you can supply a directory of fastq files generated using PacBio's CCS generating software. If trimmed to ensure that the ends of the reads are high quality, heteroduplex reads can be detected if single nucleotides have very low QVs in a read that is otherwise high quality. By specifying a QV cutoff with `-q`, fastq-hdfinder will report the number of reads in each file that have any bases below the QV cutoff. 

While the absolute values are not very precise, it should give values that are relatively comparable. This method can be used to identify whether sequences from PCR reactions with high template diversity have more heteroduplexes than sequences from PCR reactions with low template diversity.

### SSC method
The SSC method is much more precise. `ssc-hdfinder.py` looks through the output of pbssc, and identifies all the reads for which both strands have available post-filter sequences. It then compares the sequences of each strand, and identifies heteroduplexes and homoduplexes.

By default, ssc-hdfinder considers molecules that contain only single indels separated by at least 20bp as homoduplexes. This is to prevent sequencing errors from dominating the reported frequencies of heteroduplexes. This can be bypassed using the `-s` flag or by setting the window size using `-w`.

## Notes

### Trimming
- If amplicon sequences need to be trimmed, it is recommended to use pbssc for this task, especially if using the filtering options. Because the ends of the sequences may have lower consensus accuracy, pbssc trims before filtering.
