get.gene.bed <- function(file.gene, file.tres, step=1000){

    refGene <- read.table(file.gene)
    if(ncol(refGene) < 6) stop("Genomic bed data should be in BED6 format")
    
    refGene <- refGene[,1:6]
    colnames(refGene) <- c("gene.chrom", "gene.chromStart", "gene.chromEnd", "transcript.id","gene.name","gene.strand")
    refGene <- refGene[grep("random|Un|hap", refGene$gene.chrom, invert=TRUE),]

    #filter chromosome without bigwigs using TREs
    refGene     <-refGene[refGene[,1] %in% unique(read.table(file.tres)[,1]),]

    refGene.ex  <- refGene[(refGene$gene.chromEnd-refGene$gene.chromStart)<=step,]
    bodies      <- refGene[(refGene$gene.chromEnd-refGene$gene.chromStart)>step,]
    
    bodies$gene.chromStart[bodies$gene.strand == "+"] <- bodies$gene.chromStart[bodies$gene.strand == "+"]+500
    bodies$gene.chromEnd[bodies$gene.strand == "-"]   <- bodies$gene.chromEnd[bodies$gene.strand == "-"]-500
    bodies <-unique(bodies)    
    
    list(bodies, refGene.ex)
}   

get.tres.bed <- function(file.tres){

    TREs    <- read.table(file.tres)[,1:3]
    TREs    <- center.bed(TREs,250,250)
    TREs[,2]<- sapply(TREs[,2],function(x)max(x,0))
    
    colnames(TREs) <- c("tre.chrom", "tre.chromStart", "tre.chromEnd")
    TREs
}

fun.counts <- function(file.plus, file.minus, intervals) {
        
    pl  <- load.bigWig(file.plus)
    mn  <- load.bigWig(file.minus)
    
    on.exit(unload.bigWig(pl))
    on.exit(unload.bigWig(mn))
        
    if (dim(intervals)[2] == 3){
       cts.pl    <- bed.region.bpQuery.bigWig(pl, intervals, abs.value = TRUE)
       cts.mn    <- bed.region.bpQuery.bigWig(mn, intervals, abs.value = TRUE)
       counts    <- cts.pl+ cts.mn
    }else if(dim(intervals)[2] == 6){
       counts    <- bed6.region.bpQuery.bigWig(pl, mn, intervals, abs.value = TRUE)
    }else{
       stop( "*: Genomic intervals should be in BED3 or BED6 format!") 
    }
        
    counts
}

get.raw.counts <- function(region.bed, files.query.plus, files.query.minus, files.control.plus, files.control.minus, ncores =3 ){

    files.plus    <- c(files.query.plus, files.control.plus)
    files.minus   <- c(files.query.minus, files.control.minus)
    
    #count reads between each region of TRE
    raw.counts <- do.call(cbind,  mclapply(1:length(files.plus), function(i) 
                                    fun.counts( file.plus   = files.plus[i], 
                                                file.minus  = files.minus[i],
                                                intervals   = region.bed),
                                                mc.cores    = ncores))
    
    exp.names  <- c(sapply(1:(length(files.query.plus)), 
                                function(i){paste("query.rep",i, sep="") } ),
                    sapply(1:(length(files.control.plus)), 
                                function(i){paste("control.rep",i, sep="") } ))

    colnames(raw.counts) <- exp.names

    raw.counts
}

dif.deseq2 <- function(files.query.plus, files.query.minus, files.control.plus, files.control.minus, region.bed, out.file, ncores=1){
    
    ## the data frame for deseq analysis
    raw.counts <- get.raw.counts(region.bed, files.query.plus, files.query.minus, files.control.plus, files.control.minus, ncores= ncores)
    
    colData    <- data.frame(Condition  = c(rep("query",  length(files.query.plus)), 
                                            rep("control",length(files.control.plus))),
                             row.names  = colnames(raw.counts))    

    ## do deseq analysis
    dds <- DESeqDataSetFromMatrix( countData = raw.counts, 
                                   colData   = colData, 
                                   design    = ~ Condition)
    dds$Condition <- relevel( dds$Condition, ref="control")
    dds      <- DESeq(dds)
    res      <- results(dds)
    
    rst.dif  <- cbind.data.frame(region.bed , res)
  
    rst.dif
}


dif.deseq2.tres <- function(files.query.plus, files.query.minus, files.control.plus, files.control.minus, TRE.file, out.file, ncores=1){
    
    ## searching regions
    tres.region     <- get.tres.bed(TRE.file)

    ## deseq analysis
    tres.tab <- dif.deseq2(  files.query.plus, files.query.minus, files.control.plus, files.control.minus, tres.region, out.tres, ncores = ncores)

    write.table( tres.tab, file= out.file, col.names=T, row.names=F, quote=F, sep="\t")

    tres.tab
}

dif.deseq2.gene <- function(files.query.plus, files.query.minus, files.control.plus, files.control.minus, gene.file, TRE.file, out.file, ncores=1){
    
    ## searching regions
    gene.region <- get.gene.bed(gene.file, TRE.file, step=600)
    gene.bodies <- gene.region[[1]]
    gene.ex     <- gene.region[[2]]

    ## deseq analysis
    gene.EL <- dif.deseq2(files.query.plus, files.query.minus, files.control.plus, files.control.minus, gene.bodies, out.file, ncores=ncores)
    
    gene.ex[,colnames(gene.EL)[-c(1:6)]]<-NA
    gene.tab <-rbind.data.frame(gene.EL, gene.ex)

    ## output
    gene.deseq.bed <-tempfile(tmpdir=".", fileext=".bed")
    on.exit(unlink(gene.deseq.bed))
    write.table(gene.tab, file= gene.deseq.bed, col.name=F, row.name=F, sep="\t", quote=F)

    gene.tab <- try( read.table(pipe(paste( "LC_COLLATE=C sort -k1,1 -k2,2n", gene.deseq.bed))), silent=TRUE)
    if (class(gene.tab)=="try-error")
       return(NULL)
   
    colnames(gene.tab) <- c("gene.chrom", "gene.chromStart", "gene.chromEnd",
                            "transcript.id","gene.name","gene.strand",
                            "gene.baseMean","gene.log2FoldChange", 
                            "gene.lfcSE", "gene.stat", "gene.pvalue", "gene.padj"); 
    ## save result
    write.table(gene.tab, file=out.file,
                          col.names=T,row.names=F,quote=F,sep="\t")   

    gene.tab
}

diffTXN <- function( files.query.plus, files.query.minus, files.control.plus, files.control.minus, file.gene, file.tres, out.prefix=NULL, ncores=3){

    ## output file 
    if(is.null(out.prefix)) out.prefix <- Sys.Date( )
    out.tres     <- paste(out.prefix, ".TRE.deseq.txt", sep="")
    out.gene     <- paste(out.prefix, ".gene.deseq.txt", sep="")

    ## deseq analysis
    tab.dif.gene <- dif.deseq2.gene(  files.query.plus, files.query.minus, files.control.plus, files.control.minus, file.gene, file.tres, out.gene, ncores = ncores)
    tab.dif.tres <- dif.deseq2.tres(  files.query.plus, files.query.minus, files.control.plus, files.control.minus, file.tres, out.tres, ncores = ncores)
    
    tftar.obj <- list(  tab.dif.tres        = tab.dif.tres,
                        tab.dif.gene        = tab.dif.gene,
                        files.query.plus    = files.query.plus, 
                        files.query.minus   =files.query.minus, 
                        files.control.plus  =files.control.plus, 
                        files.control.minus =files.control.minus,
                        TRE.file            = file.tres, 
                        gene.file           = file.gene,               
                        ncores= ncores)
    
    class(tftar.obj)<-"tfTarget"

    tftar.obj
}
