
library(tfTarget)

#system variables


Half.size = 150
Min.size = 2500
Exp.cutoff = 2
Sites.num.cutoff = 10

#default variables

BigWig.path = "./"
Ncores = 1
MTH = 7
Run.repeats = 2
Pval.cutoff.up = 0.01
Pval.cutoff.down = 0.1
Fdr.cutoff = 0.01
Closest.N = 2
Distance.cutoff = 50E3
Gene.pval.cutoff = 0.05

deseq.only = F
rtfbsdb.only = F
Tfs.path =system.file("extdata", "tfs.rdata", package="tfTarget")

args<-commandArgs(TRUE)

tag.pos <- which(sapply(args,FUN=function(x) substr(x,1,1)=="-" ))


#check variables not understand,e.g. misspelling

all.args<-c("-query","-control","-prefix","-TRE.path","-gene.path",
		"-2bit.path","-bigWig.path","-ncores","-deseq","-rtfbsdb",
		"-tfs.path","-mTH","-cycles","-pval.up","-pval.down","-fdr.cutoff",
		"-dist","-closest.N","-pval.gene")
all.tags<-args[tag.pos]

unmatched.tags <- all.tags[!all.tags %in% all.args]

if(length(unmatched.tags)>0) {
	print("Error! The following tags are not understood. Please check the spelling.")
	stop(unmatched.tags)
}

tag.pos <- c(tag.pos,length(args)+1)


#assign required arguments

Files.query = args[(which(args=="-query")+1):(tag.pos[which(tag.pos==which(args=="-query"))+1]-1)]
Files.control = args[(which(args=="-control")+1):(tag.pos[which(tag.pos==which(args=="-control"))+1]-1)]
Plus.files.query = Files.query[seq(1,length(Files.query),by=2)]
Minus.files.query = Files.query[seq(2,length(Files.query),by=2)]
Plus.files.control = Files.control[seq(1,length(Files.control),by=2)]
Minus.files.control = Files.control[seq(2,length(Files.control),by=2)]
Prefix = args[which(args == "-prefix")+1]
TRE.path = args[which(args=="-TRE.path")+1]
Gene.path = args[which(args=="-gene.path")+1]
File.twoBit = args[which(args=="-2bit.path")+1]


#assign optional arguments

#optional system parameters
if ("-bigWig.path" %in% args ) BigWig.path = args[which(args=="-bigWig.path")+1]
if ("-ncores" %in% args ) Ncores = as.numeric(args[which(args=="-ncores")+1])
if ("-deseq" %in% args ) deseq.only = T
if ("-rtfbsdb" %in% args ) rtfbsdb.only = T
if ("-tfs.path" %in% args ) Tfs.path = args[which(args=="-tfs.path")+1]

#optional rtfbsdb parameters
if ("-mTH" %in% args ) MTH = as.numeric(args[which(args=="-mTH")+1])
if ("-cycles" %in% args ) Run.repeats = as.numeric(args[which(args=="-cycles")+1])
if ("-pval.up" %in% args ) Pval.cutoff.up = as.numeric(args[which(args=="-pval.up")+1])
if ("-pval.down" %in% args ) Pval.cutoff.down = as.numeric(args[which(args=="-pval.down")+1])
if ("-fdr.cutoff" %in% args ) Fdr.cutoff = as.numeric(args[which(args=="-fdr.cutoff")+1])


#optional TRE to gene parameters
if ("-dist" %in% args ) Distance.cutoff = as.numeric(args[which(args=="-dist")+1])

if ("-closest.N" %in% args ) {
	if (args[which(args=="-closest.N")+1]=="off") {
		Closest.N = NULL } else {
		Closest.N = as.numeric(args[which(args=="-closest.N")+1])
	}
}	

if ("-pval.gene" %in% args ) {
	if (args[which(args=="-pval.gene")+1]=="off") {
		Gene.pval.cutoff = NULL } else {
		Gene.pval.cutoff = as.numeric(args[which(args=="-pval.gene")+1])
	}
}


load(Tfs.path)



print("input file info")
cat("-bigWig.path:", BigWig.path, "\n")
cat("-query, plus files=", Plus.files.query, "\n")
cat("-query, minus files=", Minus.files.query, "\n")
cat("-control, plus files=", Plus.files.control, "\n")
cat("-control, minus files=", Minus.files.control, "\n")
cat("-TRE.path=", TRE.path, "\n")
cat("-gene.path=", Gene.path, "\n")


cat("\n")
print("motif enrichment parameters");
cat("-tfs.path=", Tfs.path, "\n");
cat("-mTH=", MTH, "\n");
cat("-cycles=", Run.repeats, "\n");
cat("-fdr.cutoff=", Fdr.cutoff, "\n");

cat("\n")
print("TRE-gene parameters");

cat("-dist=", Distance.cutoff, "\n")

if(is.null(Closest.N)) {
	cat("-closet.N= off \n")
} else {
	cat("-closet.N=", Closest.N, "\n")
}

if(is.null(Gene.pval.cutoff)) {
	cat("-pval.gene= off \n")
} else {
	cat("-pval.gene=", Gene.pval.cutoff, "\n")
}

cat("\n")
cat("running using ", Ncores, " cores", "\n")


#get deseq tables
print("running DESeq2")
tfTar <- diffTXN (TRE.path, Gene.path, BigWig.path, Plus.files.query, Plus.files.control, Minus.files.query, Minus.files.control, Prefix, Ncores)

if(!deseq.only){

	#run motif enrichment tests
	print("running rtfbsdb")
	tfTar <- searchTFBS(tfTar, tfs, File.twoBit, Pval.cutoff.up, Pval.cutoff.down, Half.size, MTH, Min.size, Run.repeats, ncores=Ncores)


	#filter and cluster motifs (link motif to TRE)
	print( "filtering motifs" )
	tfTar <- filter.rtfbsdb(tfTar, Fdr.cutoff, Sites.num.cutoff, Exp.cutoff, ncores=Ncores)

	print("plotting motifs")
	#plot motif clustering and enrichment
	plot( tfTar, Prefix)

	if(!rtfbsdb.only){
		#link TF, TRE and Gene
		print("associating TFs to TREs and genes")
		tfTar <- mapTF( tfTar, Prefix, Distance.cutoff, Closest.N, Gene.pval.cutoff);
	}
}

