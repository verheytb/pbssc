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
1. Using SMRT Analysis, set up a job with the RS_Resequencing_Barcode.1 method to align subreads to the reference sequence. Use an appropriate reference sequence for the target amplicon, and ensure that the barcode FASTA file meets PacBio specifications. Since the reads must align to the reference, pbssc only generates single-stranded reads for variants of a known amplicon sequence.

2. Take note of the job number and directory; it should be in your specified `userdata` directory that was specified in the SMRT Analysis installation. An example job directory would be: `/opt/smrtanalysis/userdata/jobs/016/016531`

2. Run `pbbarcode labelZmws` from the SMRT Shell to filter reads by (replace `/opt/smrtanalysis` with your SMRT Analysis install directory):
```sh
/opt/smrtanalysis/smrtcmds/bin/smrtshell
pbbarcode labelAlignments 
```
3. Run `pbssc.py` using the SMRT Analysis job as the input directory (eg. "~/userdata/jobs/016/016531").

```sh
python pbssc.py -j ~/userdata/jobs/016/016592 -r reference_amplicon.fasta -o outdir -t ATCTTCGATCGA,TGTAACTGAAGA
```
