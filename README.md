tfTarget
========

Transcription factors (TFs) regulate complex programs of gene transcription by binding to short DNA sequence motifs. Here we introduce tfTarget, a unified framework that identifies the "transcription factor -> regulatory elements -> target gene" regulatory network that are differential regulated between two conditions using PRO-seq/GRO-seq/ChRO-seq data as input.

integrates a database of more than 65,000 TF binding motifs with tools to easily and efficiently scan target genome sequences. Rtfbsdb clusters motifs with similar DNA sequence specificities and optionally integrates RNA-seq or PRO-seq data to restrict analyses to motifs recognized by TFs expressed in the cell type of interest.  Our package allows common analyses to be performed rapidly in an integrated environment.  

Uses: Parse TF motifs from public databases, read into R, and scan using 'rtfbs'.

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
https://github.com/Danko-Lab/rtfbs_db/blob/master/rtfbsdb-vignette.pdf (Coming soon)

* R manual:
https://github.com/Danko-Lab/rtfbs_db/blob/master/rtfbsdb-manual.pdf 	 (Coming soon)

