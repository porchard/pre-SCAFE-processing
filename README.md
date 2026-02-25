# pre-SCAFE processing

This is code that can be run to prepare ONT-sequenced 10x Genomics 3' GEX datasets for SCAFE. It assumes one has already run the official ONT pipeline on the dataset.

## Steps

The pipeline consists of two steps.

The first step consists of several substeps, which together ensure the BAM file meets various expectations of SCAFE:

1. The CIGAR string is updated to remove any hardclipping.
2. The amount of softclipping at the beginning of the read is determined, and the read is skipped if there are more than 4 bps of softclipping. The official ONT pipeline leaves the 'GGG' trinucleotide from the last 3 basepairs of the TSO sequence at the 5' of the read (which will be trimmed below), and the mRNA cap is expected to result in a fourth 'G' (which SCAFE will look for). Therefore, many reads have 4 bps of softclipping. Softclipping is not well-handled within SCAFE, so we filter reads with extensive softclipping.
3. The 5’ end of the read is expected to begin with ‘GGG’ (the last 3 basepairs of the TSO sequence; the rest of the TSO sequence is trimmed earlier in the pipeline). If it does not (allowing for one mismatch) the read is filtered out; if it does, the three bps at the 5’ end of the read are trimmed.
4. The read is trimmed from the 3’ end until it is no longer than 100 bps.

In the second step -- which is optional but recommended -- we handle cases where multiple reads from the same nucleus and with the same UMI and gene assignment suggest different 5' ends. These likely represent technical artifacts, and can inflate the supposed support each TSS cluster gets in the SCAFE filtering step, as well as inflate SCAFE tCRE counts. We've found this to noticeably impact SCAFE results obtained from some libraries, and to have minimal impact for other libraries. To do this, we identify groups of reads with the same cell barcode (CB tag), UMI (UB tag), and gene assignment (GX tag), and check the 5' end positions. If they do not all nominate the same 5' end position, we determine the 5' end position with the most reads supporting it and keep only the corresponding reads. If two 5' end positions have equal support, we keep the most upstream one.


## Example data

ONT has published an example ONT-sequenced 10x Genomics 3' GEX dataset here: https://epi2me.nanoporetech.com/sc-gemx-2025.02/#dataset

If you have `aws` CLI tool on your machine, you can download a portion of their (processed) dataset, which includes the BAM file, as follows:

```
aws s3 sync --no-sign-request s3://ont-open-data/sc_gemx_2025.02/analysis/workflow_outputs/293t/293t_1_1/PBC16207_20250123_293t_1_1 PBC16207_20250123_293t_1_1
```

## Running

There are two ways to run the code -- as stand-alone python scripts, or as a small NextFlow pipeline.

### Running as stand-alone python scripts

If you wish to simply run the python scripts directly, you'll need python3 and samtools installed, as well as the following python packages: `pysam`, `numpy`, `pandas`, `matplotlib`, and `seaborn`.

Using the example ONT BAM file above, you can run the two steps as follows:

```
# step 1
# this uses minimal memory, and should run fine on a laptop
# should take ~2.5 hours on the example dataset
python3 /path/to/bin/trim-for-scafe.py --input-bam /path/to/PBC16207_20250123_293t_1_1/PBC16207_20250123_293t_1_1.tagged.bam --output-bam trimmed.unsorted.bam --trim-to 100 --ggg-mismatches-allowed 1 --max-softclipping 4 2> log.txt
samtools sort -m 3G -o trimmed.bam trimmed.unsorted.bam
samtools index trimmed.bam

# step 2
# this uses considerable memory (~70-75GB on the example dataset) so is best run on a server
# should take ~45 minutes on the example dataset
python3 /path/to/bin/filter-bam-to-most-supported-5prime-ends.py --bam-in trimmed.bam --bam-out filtered.bam --prefix filtering.
```

`filtered.bam` is the file that should be used for SCAFE input.

### Running as a NextFlow pipeline

If you have NextFlow and Singularity installed (and NextFlow is appropriately configured for your system), you can instead run the NextFlow pipeline. In this case, the pipeline utilizes Singularity containers for all other dependencies (there is no need to install python, python packages, or samtools).

The NextFlow pipeline needs only the ONT BAM file as input. With the example dataset, the preprocessing can be run as follows:

```
nextflow run -resume --results results --bam_glob '/path/to/PBC16207_20250123_293t_1_1/PBC16207_20250123_293t_1_1.tagged.bam' /path/to/main.nf
```

There will be two output subdirectories in the `results` directory:
1. `trimmed`. This is the output of step 1 above.
2. `filtered`. This the the output of step 2 above. These BAM files should be used for SCAFE input.
