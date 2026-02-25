# pre-SCAFE processing

This is a small NextFlow pipeline that can be run to prepare ONT-sequenced 10x Genomics 3' GEX datasets for SCAFE. It assumes one has already run the official ONT pipeline on the dataset.

## Steps

The pipeline consists of two steps.

The first step consists of several substeps, which together ensure the BAM file meets various expectations of SCAFE.

1. The CIGAR string is updated to remove any hardclipping.
2. The amount of softclipping at the beginning of the read is determined, and the read is skipped if there are more than 4 bps of softclipping. The official ONT pipeline leaves the 'GGG' trinucleotide from the last 3 basepairs of the TSO sequence at the 5' of the read (which will be trimmed below), and the mRNA cap is expected to result in a fourth 'G' (which SCAFE will look for). Therefore, many reads have 4 bps of softclipping. Softclipping is not well-handled within SCAFE, so we filter reads with extensive softclipping.
3. The 5’ end of the read is expected to begin with ‘GGG’ (the last 3 basepairs of the TSO sequence; the rest of the TSO sequence is trimmed earlier in the pipeline). If it does not (allowing for one mismatch) the read is filtered out; if it does, the three bps at the 5’ end of the read are trimmed.
4. The read is trimmed from the 3’ end until it is no longer than 100 bps.

In the second step, we handle cases where multiple reads from the same nucleus and with the same UMI and gene assignment suggest different 5' ends. These likely represent technical artifacts, and can inflate the supposed support each TSS cluster gets in the SCAFE filtering step, as well as inflate SCAFE tCRE counts. We've found this to noticeably impact SCAFE results obtained from some libraries, and to have minimal impact for other libraries. To do this, we identify groups of reads with the same cell barcode (CB tag), UMI (UB tag), and gene assignment (GX tag), and check the 5' end positions. If they do not all nominate the same 5' end position, we determine the 5' end position with the most reads supporting it and keep only the corresponding reads. If two 5' end positions have equal support, we keep the most upstream one.


## Example data

ONT has published an example ONT-sequenced 10x Genomics 3' GEX dataset here: https://epi2me.nanoporetech.com/sc-gemx-2025.02/#dataset

If you have `aws` CLI tool on your machine, you can download a portion of their (processed) dataset, which includes the BAM file, as follows:

```
aws s3 sync --no-sign-request s3://ont-open-data/sc_gemx_2025.02/analysis/workflow_outputs/293t/293t_1_1/PBC16207_20250123_293t_1_1 PBC16207_20250123_293t_1_1
```

## Dependencies

You must have NextFlow and Singularity installed. NextFlow should be configured as appropriate for your system. The pipeline utilizes Singularity containers for all other dependencies.

## Running

The NextFlow pipeline needs only the ONT BAM file as input. With the example dataset, the preprocessing can be run as follows:

```
nextflow run -resume --results results --bam_glob '/path/to/PBC16207_20250123_293t_1_1/PBC16207_20250123_293t_1_1.tagged.bam' /path/to/main.nf
```

## Output

There will be two subdirectories output:
1. `trimmed`. This is the output of step 1 above.
2. `filtered`. This the the output of step 2 above. These BAM files should be used for SCAFE input.
