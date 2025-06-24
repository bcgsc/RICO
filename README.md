# RICO : Rna-seq Immune COntainer

Immune profiling from RNA-seq data.
Uses STAR/RSEM, and immunedeconv to provide estimates of cell populations in short read RNA-Seq data.

## Requirements

Nextflow version 20.12.0 or later is required.

Container building/testing was done with Apptainer version 1.3.1-1, but Singularity (v3.5.2+) should peform equally well.

During your first run in a new storage location an internet connection is required for the 
associated containers to be pulled down.  

To run the required cibersort component of the workflow you need to get permission from the [cibersort
team](https://cibersort.stanford.edu/). Note that this is the original version of CIBERSORT, **not
CIBERSERTx**. Once permission is granted (it may take multiple days to get approval) you can download
the required scripts: `CIBERSORT.R` and `LM22.txt`.  Add the path to your copies of these scripts
to the `nextflow.config` file (we recommend placing them in the `scripts` directory of this repository).
If these script paths are not provided, RICO will skip this step automatically.

### Space

The container itself requires approximately 28GB of space fully built.

The amount of storage required to analyze a sample is a function of the number and length of reads
provided for analysis.

We have seen up to 160Gb of space required to analyze a sample with 750 million pairs of 150bp reads.
More typically we see space requirement of 30-40Gb for up to 200 milllion read pairs.

### CPU and RAM

STAR genome indices generation, during container build step, can require a significant amount of RAM (32GB+).
It has been configured to use 16 threads, but can be adjusted in the `rico.recipe` file.

Typical CPU usage is up to 300 CPU hours.

Peak RAM usage is 48Gb during the STAR/RSEM process.

### I/0

A typical sample will have a total network and disk traffic of up to 600Gb (500Gb of read operations and 100Gb of write).

## Setup

A **one-time** build of the workflow container is required. This involves the installation of necessary dependencies and the transferring of scripts and references (`rico.recipe`).

Firstly, a number of external genome references are required. Please download these files and place them in the`ref/` directory.
```
hg38_no_alt.chrlist
hg38_no_alt.fa
hg38_no_alt.fa.fai
hg38_no_alt.grp
hg38_no_alt.idx.fa
hg38_no_alt.n2g.idx.fa
hg38_no_alt.seq
hg38_no_alt.ti
hg38_no_alt.transcripts.fa
Homo_sapiens.GRCh38.100.remapped.gtf
```

Then, run the following command...
```
apptainer build </path/to/rico_container.sif> rico.recipe
```
*depending on your server permissions, you may need to include `--no-https` or `--fakeroot` parameters in the command*

Finally, update `params.rico_container` in `nextflow.config` to the path to your container.

## How To Run

Samples to analyze need an ID, and paths to the paried end fastq files.  
The **Header Line is Required!** and must be `library_id,fastq1,fastq2`.
Multiple samples can be provided in this file.
For example...

```
cat samples.csv

library_id,fastq1,fastq2
PROABC,AHF52WDRXX_294571.fastq.gz,AHF52WDRXX_294571.fastq.gz
PROBCD,20_3310_R1.fastq.gz,20_3310_R2.fastq.gz
```

Then, to run the analysis...
`nextflow run rico_workflow.nf --out_dir output_folder --samples_file samples.csv`.

If nextflow is properly installed, the above command will pull down the analysis container and run on the local machine.
For most first time users there is no need to download the container manually.
If you want to run the analysis on a cluster you can modify the lines in the `nextflow.config` file to provide the name
of your executor (see [Nextflow Executors](https://www.nextflow.io/docs/latest/executor.html) ).  Similar options
are available for [AWS](https://www.nextflow.io/docs/latest/awscloud.html), [Google Cloud](https://www.nextflow.io/docs/latest/google.html),
and [Azure](https://www.nextflow.io/docs/latest/azure.html).

## Output

A Separate folder will be created for each sample listed in the .csv file.   These folders will be in the
folder provided to the `--out_dir` parameter at run time.    Inside each output folder,  files
are organized under folders called `immunedeconv`, `cibersort`, `M1M2`, and `IPASS`.  Each analysis will output
a .csv file corresponding to each tool that is run.  

Each output .csv file will have 3 columns simila to the following:

| cell_type | Score | ptile |
| --------- | ----- | ----- |
| T cell | 6.57 | 22 |
| B cell | 0.23 | 2 |
| ... | ... | ... |

The rows reported are unique to each of the tools as are the way the scores are calcualated.   No post processing
is perfomed on the scores after they are reported from each tool.  The `ptile` column reports the percentile
of each score when compared to the provided matrix of control samples.   Currently, there is a matrix of
222 RNA libraries included in the container for this purpose.

## Reference Details

This workflow uses a version of hg38 that doesn't include either alt contigs or multiple HLA alleles.
Ensembl version 100 gene models are used during alignment and are provided to the cell estimation
scripts in immunedeconv.

Once the main container is built you can access the reference files inside the container at
`/rico/ref`.

## Known issues

None at this time
