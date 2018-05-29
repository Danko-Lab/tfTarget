tfTarget
========

Transcription factors (TFs) regulate complex programs of gene transcription by binding to short DNA sequence motifs. Here we introduce rtfbsdb, a unified framework that integrates a database of more than 65,000 TF binding motifs with tools to easily and efficiently scan target genome sequences. Rtfbsdb clusters motifs with similar DNA sequence specificities and optionally integrates RNA-seq or PRO-seq data to restrict analyses to motifs recognized by TFs expressed in the cell type of interest.  Our package allows common analyses to be performed rapidly in an integrated environment.  

Uses: Parse TF motifs from public databases, read into R, and scan using 'rtfbs'.

<img src="img/FIG1.png">

Sources of Position Weight Matrices (PWMs).
-------------------------------------------


Requires
--------

* rtfbdbs: In R, type: "install.packages('rtfbs')"

* bedops:
	* Get the latest version of the bedops binaries here: https://bedops.readthedocs.org/en/latest/
	* Install, and add them to your path.

* The twoBitToFa program from the Kent libraries.  Download it here: http://hgdownload.cse.ucsc.edu/admin/exe/

* 2bit files for your genome of interest.  Find links to these here: http://hgdownload.cse.ucsc.edu/downloads.html

Installation
--------

* If the package has been installed in R, please use the following codes to update the package. 

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

