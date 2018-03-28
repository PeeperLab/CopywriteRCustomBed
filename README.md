# CopywriteRCustomBed
## This is an adapted version of CopywriteR (https://github.com/PeeperLab/CopywriteR). 

CopywriteRCustomBed is an adapted version of CopywriteR developed for a specific project. 

CopywriteRCustomBed allows users to provide a custom bed file for each samples with the 
peak calling to identify the on-target sequence reads. All other functionality is identical 
to CopywriteR. 




Current methods for detection of copy number variants and aberrations (CNV and
CNA) from targeted sequencing data are based on the depth of coverage of
captured exons. Accurate CNA determination is complicated by uneven genomic
distribution and non-uniform capture efficiency of targeted exons. Here we
present CopywriteR which eludes these problems by exploiting ‘off-target’
sequence reads. CopywriteR allows for extracting uniformly distributed copy
number information, can be used without reference and can be applied to
sequencing data obtained from various techniques including chromatin
immunoprecipitation and target enrichment on small gene panels. CopywriteR
outperforms existing methods and constitutes a widely applicable alternative to
available tools.

CopywriteR has been described in
[Kuilman et al., 2015](http://genomebiology.com/2015/16/1/49/abstract). 

## Installation (not via Bioconductor)

CopywriteR can be installed without Bioconductor, which is useful when you have
R 3.1 and Bioconductor 3.0 installed. Installation should be performed as
follows:

    > source("http://bioconductor.org/biocLite.R")
    > biocLite(c("matrixStats", "gtools", "data.table", "S4Vectors", "chipseq",
                 "IRanges", "Rsamtools", "DNAcopy", "GenomicAlignments",
                 "GenomicRanges", "GenomeInfoDb", "BiocParallel",
                 "futile.logger"))

In addition, CopywriteRCustomBed requires the CopyhelpeR package, which can be downloaded
as tarball (.tar.gz-file) from the
[CopyhelpeR releases webpage](https://github.com/PeeperLab/CopyhelpeR/releases).
Subsequently, the CopyhelpeR package can be installed from the command line
using the following command:

    $ R CMD INSTALL CopyhelpeR*.tar.gz

As the last step in the installation process, the latest CopywriteR package can
be downloaded from the
[CopywriteR releases webpage](https://github.com/PeeperLab/CopywriteR/releases)
and installed using the following command:

    $ R CMD INSTALL CopywriteRCustomBed*.tar.gz

Now you are all set to start your analysis.

## CopywriteR usage:

Load the CopywriteR package in R using:

    > library("CopywriteRCustomBed")

CopywriteR contains three main functions:

`preCopywriteR` will generate a GRanges object containing mappability and
GC-content information, and one containing 'blacklisted' regions that contain
are subject to copy number variation. These 'helper' files can be created for
any specified bin size that is a multiple of 1000 bp, and for any of the
available reference genomes (hg18, hg19, hg38, mm9 and mm10). The helper files
can be re-used and need to be created only once for every combination of
reference genome and bin size. `preCopywriteR` uses information stored in
pre-assembled 1kb bin mappability and GC-content GRanges objects to create the
custom bin size helper files. These objects are stored in the CopyhelpeR
annotation package.

preCopywriteR can be run as follows:

    > preCopywriteR(output.folder, bin.size, ref.genome, prefix = "")

`CopywriteR` will generate separate tables with compensated read counts and
log2-transformed normalized compensated read counts after correction for
GC-content, mappability and upon removal of blacklisted regions.

    > CopywriteRCustomBed(sample.control, Custom.bed, destination.folder, reference.folder,
                       bp.param, capture.regions.file,
                       keep.intermediary.files = FALSE, bpWidth = 0)

`plotCNA` performs segmentation using the DNAcopy Bioconductor package, and
plotting of copy number profiles.

    > plotCNA(destinationFolder, set.nchrom)

For more details please refer to the CopywriteR vignette. Alternatively, one of
the following commands can be used to show help files for the corresponding
function:

    > ?preCopywriteR
    > ?CopywriteRCustomBed
    > ?plotCNA

## Quick start

A typical analysis using CopywriteR could be as follows. First, CopywriteR needs
to be loaded:

    > library(CopywriteRCustomBed)

Then, preCopywriteR can be run using the command:

    > preCopywriteR(output.folder = file.path("./path/to/output/folder"),
                    bin.size = 20000,
                    ref.genome = "mm10",
                    prefix = "")

Next, we need to specify the settings for parallel computing. We have
implemented use of the
[BiocParallel](http://bioconductor.org/packages/release/bioc/html/BiocParallel.html)
package, which supports several different types of environments. For every
environment, a BiocParallelParam can be specified that defines how parallel
computation is executed. Below, we use a SnowParam instance of
BiocParallelParam, which is based on the
[snow](http://cran.r-project.org/web/packages/snow/index.html) package. Please
refer to the
[BiocParallel](http://bioconductor.org/packages/release/bioc/html/BiocParallel.html)
package for more information. A SnowParam using 12 CPUs can be defined as
follows:

    > bp.param <- SnowParam(workers = 12, type = "SOCK")

Next, we need to specify which samples, controls and custom Peak files correspond to each other
using the sample.control variable. This information will then be used to
remove on-target reads in the corresponding sample. We specify the
sample.control variable as follows:

    > path <- "./path/to/bam"
    > samples <- list.files(path = path, pattern = "bam$", full.names = TRUE)
    > controls <- samples[c(1:6, 1:6)]
    > sample.control <- data.frame(samples, controls)

This might result in the following variable:

    > sample.control
          samples   controls	CustomPeakfiles
    1  ./C003.bam ./C003.bam	./C003Peaks.bed
    2  ./C016.bam ./C016.bam	./C016Peaks.bed
    3  ./C024.bam ./C024.bam	./C024Peaks.bed
    4  ./M003.bam ./C003.bam	./M003Peaks.bed
    5  ./M016.bam ./C016.bam	./M016Peaks.bed
    6  ./M024.bam ./C024.bam	./M024Peaks.bed
    
Next, we need a file with the 'on-target' peaks in bed format. 

	 > Bedfile
	CustomPeakfiles
	1  ./C003Peaks.bed
	2  ./C016Peaks.bed
	3  ./C024Peaks.bed
	4  ./M003Peaks.bed
    5  ./M016Peaks.bed
    6  ./M024Peaks.bed
  
The 'on-target' regions can be extended on both sides with n bps. 

	>   bpWidth = 100
	

We are now set for running CopywriteR:

    > CopywriteRCustomBed(sample.control = sample.control, Custom.bed = Bedfile
                 destination.folder = file.path("./path/to/destination/folder"),
                 reference.folder = file.path("./path/to/reference/folder", "mm10_20kb"),
                 bp.param = bp.param, bpWidth = 0)

After running CopywriteRCustomBed the workflow is identical to the workflow of CopywriteR. 


## Troubleshooting
Please be aware that this version of CopywriteR was generated for a specific project and  
could contain bugs and issues when run on samples other than the samples in this specific project. 

The manuscript for which this version was generated is currently under review and a link to 
the paper will be added here as soon as is has been accepted for publication.  

If you want to use the code and have any questions feel free to email:

- [Oscar Krijgsman](mailto:o.krijgsman@nki.nl)
- [Thomas Kuilman](mailto:t.kuilman@nki.nl)










