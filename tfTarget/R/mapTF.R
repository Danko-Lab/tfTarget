get.proximal.genes <- function( enh.bed, gene.deseq.bed.sorted, distance.cutoff, closest.N, gene.pval.cutoff, up.down){

  #remove the lfcSE and stat column to avoid being interprated as bed12 format by bedtools
  
  enh.bed<-enh.bed[,-c(6,7)]
  gene.deseq.bed.sorted<-gene.deseq.bed.sorted[,-c(9,10)]
  
  options(scipen =99)
  enh.bed.file<-tempfile()
  write.table(enh.bed,file = enh.bed.file, col.name=F, row.name=F, sep="\t", quote=F)
  
  gene.bed.file<-tempfile()
  write.table(gene.deseq.bed.sorted,file = gene.bed.file, col.name=F, row.name=F, sep="\t", quote=F)

  if(!is.null(closest.N)){
  	#find closest nth gene within the distance cutoff

    get.closet.nth.command<-paste("cat" ,  gene.bed.file, "| awk 'BEGIN{OFS=\"\t\"} {print $1,$6==\"+\"?$2:$3-1,$6==\"+\"?$2+1:$3,$4,$5,$6,$7,$8,$9,$10}' | sort-bed - | bedtools closest -t \"first\" -d -k", closest.N  ,"-a", enh.bed.file, "-b stdin")

    closest.tab <- read.table(pipe(get.closet.nth.command))
    unlink(enh.bed.file)
    unlink(gene.bed.file)

    first.N <- rep(c(1: closest.N),nrow(closest.tab)/closest.N)
    closest.tab$closest.N <- first.N
    distance <- closest.tab[ ,(ncol(closest.tab)-1)]
    closest.tab <- closest.tab[distance <= distance.cutoff, ]
  	
  }
  
  else{
    #only use distance cutoff (default)
    
    enh.ext.bed.file<-tempfile()
    enh.ext.bed<-enh.bed[,1:3]
    enh.ext.bed[,2]<-sapply(enh.ext.bed[,2],function(x)max(0,x- distance.cutoff))
    enh.ext.bed[,3]<-enh.ext.bed[,3]+ distance.cutoff
    
    enh.ext.bed<-cbind.data.frame(enh.ext.bed,enh.bed)
    
    write.table(enh.ext.bed, file = enh.ext.bed.file, col.name=F, row.name=F, sep="\t", quote=F)

    get.proximal.command<-paste("cat", gene.bed.file, "| awk 'BEGIN{OFS=\"\t\"} {print $1,$6==\"+\"?$2:$3-1,$6==\"+\"?$2+1:$3,$4,$5,$6,$7,$8,$9,$10}' | sort-bed - | bedtools closest -t \"all\" -d -a", enh.ext.bed.file, "-b stdin | awk 'BEGIN{OFS=\"\t\"} $NF==0 {$NF=\"\"; print $0}' | bedtools overlap -i stdin -cols 5,6,12,13")

    closest.tab <- read.table(pipe(get.proximal.command))
    closest.tab<-closest.tab[,-c(1:3)]
    
    closest.tab[,ncol(closest.tab)] <- sapply(closest.tab[,ncol(closest.tab)],function(overlap) -(min(0, overlap)) )
    
    unlink(enh.bed.file)
    unlink(enh.ext.bed.file)
    unlink(gene.bed.file)
  
  }
    
  if(!is.null(gene.pval.cutoff)) {
  	closest.tab<-closest.tab[!is.na(closest.tab$V17) & closest.tab$V17<gene.pval.cutoff,]
  	if(up.down=="up") closest.tab<-closest.tab[closest.tab$V15>0,]
  	if(up.down=="down") closest.tab<-closest.tab[closest.tab$V15<0,]
  }
  
  return(closest.tab)
  
}


mapTF<-function(tfTar, out.prefix=NULL, distance.cutoff, closest.N, gene.pval.cutoff){

  if(class(tfTar)!="tfTarget")
     stop("The first parameter is not a 'tfTarget' object!");

  if(class(tfTar$tfs)!="tfbs")
     stop("No 'tfbs' object in tfTar!");

  options(scipen =99);

  #link TRE to gene
  up.gene.tab   <- get.proximal.genes( enh.bed = tfTar$enh.up.bed,
                       gene.deseq.bed.sorted = tfTar$deseq.table.gene,
                       distance.cutoff = distance.cutoff,
                       closest.N = closest.N, 
                       gene.pval.cutoff= gene.pval.cutoff,
                       up.down="up");

  down.gene.tab <- get.proximal.genes( enh.bed = tfTar$enh.down.bed,
                       gene.deseq.bed.sorted = tfTar$deseq.table.gene,
                       distance.cutoff = distance.cutoff,
                       closest.N = closest.N, 
                       gene.pval.cutoff= gene.pval.cutoff,
                       up.down="down");

  TRE.genes.tab <- rbind.data.frame(up.gene.tab, down.gene.tab);

  #link TF, TRE and Gene
  TF.TRE.tab <- tfTar$TF.TRE.tab
  TF.TRE.tab$TRE <- paste( TF.TRE.tab$tre.chrom, TF.TRE.tab$tre.chromStart, TF.TRE.tab$tre.chromEnd)
  TRE.genes.tab$TRE <- paste(TRE.genes.tab[ ,1], TRE.genes.tab[ ,2], TRE.genes.tab[ ,3])

  TF.TRE.gene.tab <- merge(TF.TRE.tab,TRE.genes.tab,by="TRE")
  TF.TRE.gene.tab.short <- TF.TRE.gene.tab[,-c(1,13:15)]
  header.vec<-c("tre.chrom", "tre.chromStart", "tre.chromEnd",
           "tf.chrom", "tf.chromStart", "tf.chromEnd", "score", "strand",
           "motif.name", "motif.id", "motif.idx",
           "TRE.baseMean", "TRE.log2FoldChange", "TRE.pvalue", "TRE.padj",
           "gene.TSS.chr", "gene.TSS.start", "gene.TSS.end", "transcript.id",
           "gene.name", "gene.strand", "gene.baseMean",
           "gene.log2FoldChange",
           "gene.pvalue", "gene.padj", "distance")
  if(!is.null(closest.N)) header.vec<-c(header.vec, "closest.N")         
           
  colnames(TF.TRE.gene.tab.short) <- header.vec       
          
  tfTar$TF.TRE.gene.tab <- TF.TRE.gene.tab.short

  if(!is.null(out.prefix))
  {
  	if(is.null(closest.N)) closest.N<-"off"
  	if(is.null(gene.pval.cutoff)) gene.pval.cutoff <-"off"
  	  	
  	txt.name <- paste( out.prefix,
                    "dist", distance.cutoff,
                    "closest.N", closest.N,
                    "gene.pval", gene.pval.cutoff,
                     sep="_");
    TF.TRE.gene.tab.filename <- paste(txt.name, ".TF.TRE.gene.txt", sep="")
    write.table( TF.TRE.gene.tab.short, file = TF.TRE.gene.tab.filename, col.names = T, row.names = F, quote = F,sep="\t" )
  }

  return(tfTar)
}
