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

3. Run `pbbarcode labelZmws` from the SMRT Shell to filter reads for those that were high-scoring during demultiplexing. Failure to do this step will result in incorrect barcode sorting.
  ```sh
  $ $SMRT_HOME/smrtcmds/bin/smrtshell
  $ pbbarcode labelAlignments --minScoreRatio 1.1 $JOBDIR/data/barcode.fofn $JOBDIR/data/aligned_reads.cmp.h5
  exit
  ```

4. Run `pbssc.py`, using the reference originally used in the SMRT Analysis Resequencing Job, and trim sequences if needed.
  ```sh
  $ python pbssc.py -j $JOBDIR -r reference_amplicon.fasta -o outdir -t ATCTTCGATCGA,TGTAACTGAAGA
  ```

## Notes

### Trimming
- If amplicon sequences need to be trimmed, it is recommended to use this program to do so, especially if using the filtering options. Because the ends of the sequences may have lower consensus accuracy, pbssc trims before filtering.
