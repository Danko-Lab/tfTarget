tfTarget
========

Transcription factors (TFs) regulate complex programs of gene transcription by binding to short DNA sequence motifs within transcription regulatory elements (TRE). Here we introduce tfTarget, a unified framework that identifies the "TF -> TRE -> target gene" regulatory network differential regulated between two conditions, e.g. experimental vs. control, using PRO-seq/GRO-seq/ChRO-seq data as input.

Workflow of tfTarget
<img src="img/img1.png">

Sources of Position Weight Matrices (PWMs).
-------------------------------------------


Requires
--------

* R packages:
	
	rphast, rtfbdbs, grid, cluster, apcluster, DESeq2, gplots.
	
	rtfbsdb (https://github.com/Danko-Lab/rtfbs_db) 
	
	bigWig  (https://github.com/andrelmartins/bigWig.git)

* bioinformatics tools or exterior command:
	
	awk, sort: Unix commands
	
	mergeBed, bedtools (http://bedtools.readthedocs.org/en/latest/)
	
	sort-bed (http://bedops.readthedocs.org/en/latest/index.html)
	
	twoBitToFa, faToTwoBit (http://hgdownload.cse.ucsc.edu/admin/exe/)

* 2bit files for your genome of interest.  Find links to these here: 
    
	http://hgdownload.cse.ucsc.edu/downloads.html

* TREs regions identified by dREG, or equivalent tools, in bed format. 
	https://github.com/Danko-Lab/dREG

* Gene annotation file in bed6 format. Can be prepared from gencode or Refseq gtf files. We recommend to use gene ID and gene name for the 4th and 5th columns. The information will show up in the output.
	https://www.gencodegenes.org/releases/current.html

Installation
--------

* If all dependent packages and commands have been installed, please use the following codes to install/update the package. 

```````
library("devtools");
install_github("Danko-Lab/tfTarget/tfTarget")
```````

Documents
----------

* R vignette:
 (Coming soon)

* R manual:
 (Coming soon)

