getCounts.TRE <- function(plus, minus, intervals) {
  pl <- load.bigWig(plus);
  mn <- load.bigWig(minus);
  counts.plus <- bed.region.bpQuery.bigWig(pl, intervals, abs.value = TRUE);
  counts.minus <- bed.region.bpQuery.bigWig(mn, intervals, abs.value = TRUE);
  counts<-(counts.plus+ counts.minus);
  unload.bigWig( pl );
  unload.bigWig( mn );
  
  return(counts);
}

get.deseq.TRE.tab<-function(TRE.path, bigWig.path, plus.files.query, plus.files.control, minus.files.query, minus.files.control, ncores){

  #prepare regions of TRE
  TREs <- read.table(TRE.path)[,1:3]
  TREs <- center.bed(TREs,250,250)
  TREs[,2]<-sapply(TREs[,2],function(x)max(x,0))


  #prepare file names
  plus.files.query.full   <- paste(bigWig.path, plus.files.query, sep="/")
  plus.files.control.full <- paste(bigWig.path, plus.files.control, sep="/")
  minus.files.query.full  <- paste(bigWig.path, minus.files.query, sep="/")
  minus.files.control.full<- paste(bigWig.path, minus.files.control, sep="/")

  plus.files  <- c(plus.files.query.full, plus.files.control.full)
  minus.files <- c(minus.files.query.full, minus.files.control.full)

  file.names <- gsub(".bw","", gsub(".bigWig","", c(plus.files.query, plus.files.control)))
  file.names <- gsub(".pl","", gsub(".plus","", file.names))
  file.names <- gsub("_pl","", gsub("_plus","", file.names))

  stopifnot(length(plus.files)==length(minus.files))

  #get reads count and do deseq
  raw_counts <- do.call(cbind,mclapply(1:length(plus.files),
       function(i) getCounts.TRE(plus.files[i], minus.files[i],intervals= TREs),
     mc.cores= ncores))

  colnames(raw_counts) <- file.names

  colData <- data.frame(Condition= c(rep("query",length(plus.files.query)), rep("control",length(plus.files.control))),
                        row.names=colnames(raw_counts))

  dds <- DESeqDataSetFromMatrix( countData= raw_counts, colData= colData, design= ~ Condition)

  ## Set the reference condition as the normal brain.
  dds$Condition <- relevel( dds$Condition, ref="control")

  dds <- DESeq(dds)

  res <- results(dds)

  EL <- cbind.data.frame(TREs, res)

  return(as.data.frame(EL))
}

getCounts.gene <- function(plus, minus, intervals) {
  pl <- load.bigWig(plus)
  mn <- load.bigWig(minus)
  counts <- bed6.region.bpQuery.bigWig(pl, mn, intervals, abs.value = TRUE)
  unload.bigWig( pl );
  unload.bigWig( mn );
  
  return(counts);
}



get.deseq.gene.tab <-function(gene.path, bigWig.path, plus.files.query, plus.files.control, minus.files.query, minus.files.control, ncores){

  #prepare regions of gene body
  refGene <- read.table(gene.path)[,1:6]
  refGene <- refGene[grep("random|Un|hap", refGene$V1, invert=TRUE),]
  refGene.excluded <- refGene[(refGene$V3-refGene$V2)<=1000,]
  refGene <- refGene[(refGene$V3-refGene$V2)>1000,]

  bodies <- refGene
  bodies$V2[bodies$V6 == "+"] <- bodies$V2[bodies$V6 == "+"]+500
  bodies$V3[bodies$V6 == "-"] <- bodies$V3[bodies$V6 == "-"]-500
  bodies <-unique(bodies)

  #prepare file names
  plus.files.query.full   <- paste(bigWig.path, plus.files.query, sep="/")
  plus.files.control.full <- paste(bigWig.path, plus.files.control, sep="/")
  minus.files.query.full  <- paste(bigWig.path, minus.files.query, sep="/")
  minus.files.control.full<- paste(bigWig.path, minus.files.control, sep="/")

  plus.files  <- c(plus.files.query.full, plus.files.control.full)
  minus.files <- c(minus.files.query.full, minus.files.control.full)

  file.names <- gsub(".bw","", gsub(".bigWig","", c(plus.files.query, plus.files.control)))
  file.names <- gsub(".pl","", gsub(".plus","", file.names))
  file.names <- gsub("_pl","", gsub("_plus","", file.names))

  stopifnot( length(plus.files) == length(minus.files) )

  #get reads count and do deseq
  raw_counts <- do.call(cbind,mclapply(1:length(plus.files),
      function(i) getCounts.gene(plus.files[i], minus.files[i],intervals= bodies),
    mc.cores= ncores))

  colnames(raw_counts)<-file.names
  colData <- data.frame(Condition= c(rep("query",length(plus.files.query)),
                                     rep("control",length(plus.files.control))),
                       row.names=colnames(raw_counts))

  dds <- DESeqDataSetFromMatrix( countData= raw_counts, colData= colData, design= ~ Condition)

  dds$Condition <- relevel( dds$Condition, ref="control") ## Set the reference condition as the normal brain.

  dds <- DESeq(dds)

  res <- results(dds)

  EL <- cbind.data.frame(refGene, res)

  refGene.excluded[,colnames(EL)[-c(1:6)]]<-NA
  gene.tab<-rbind.data.frame(EL, refGene.excluded)

  options(scipen =99)
  gene.deseq.bed <-tempfile()
  write.table(gene.tab,file= gene.deseq.bed, col.name=F, row.name=F, sep="\t", quote=F)
  gene.tab<-read.table(pipe(paste( "LC_COLLATE=C sort -k1,1 -k2,2n", gene.deseq.bed)))
  unlink(gene.deseq.bed)

  return(gene.tab)
}

diffTXN <- function( TRE.path, gene.path, bigWig.path, plus.files.query, plus.files.control, minus.files.query, minus.files.control, out.prefix=NULL, ncores=3){

  options(scipen =99);

  deseq.table.TRE  <- get.deseq.TRE.tab(TRE.path, bigWig.path, plus.files.query,
      plus.files.control, minus.files.query, minus.files.control, ncores)

  if(!is.null(out.prefix))
  {
    TRE.tab.filename <- paste( out.prefix,".TRE.deseq.bed", sep="" )
    write.table( deseq.table.TRE, file= TRE.tab.filename, col.names=F, row.names=F, quote=F, sep="\t")
  }

  deseq.table.gene<-NULL
  if(!is.null(gene.path)) {
  #if do deseq on genes
    deseq.table.gene<-get.deseq.gene.tab(gene.path, bigWig.path, plus.files.query, plus.files.control, minus.files.query, minus.files.control, ncores)
    if(!is.null( out.prefix ))
    {
      gene.tab.filename <- paste( out.prefix, ".gene.deseq.bed", sep="")
      write.table(deseq.table.gene, file=gene.tab.filename,col.names=F,row.names=F,quote=F,sep="\t")
    }
  }

  res.obj<-list(deseq.table.TRE =deseq.table.TRE,
        deseq.table.TRE = deseq.table.TRE,
        deseq.table.gene = deseq.table.gene,
        gene.path = gene.path,
        bigWig.path = bigWig.path,
        plus.files.query = plus.files.query,
        plus.files.control = plus.files.control,
        minus.files.query = minus.files.query,
        minus.files.control = minus.files.control,
        ncores= ncores)

  class(res.obj)<-"tfTarget";
  return(res.obj);
}
