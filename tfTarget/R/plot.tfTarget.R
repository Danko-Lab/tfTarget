## Write.seqfile.
## Writes a fasta file containing all sequences in a bed region.
read.seqfile.from.bed <- function(bed, twoBitPath, rm.dup = TRUE, tmpdir = getwd())
{
  # create temporary filenames
  tmp.seq = tempfile(tmpdir=tmpdir)
  tmp.fa = tempfile(tmpdir=tmpdir)

  # write sequence list
  seqList = paste(bed[,1],":", as.integer(bed[,2]), "-", as.integer(bed[,3]), sep="")

  if(length(seqList) != length(unique(seqList)))
  {
    seqStat <- table(seqList);
    if (length(which(seqStat>1)))
    {
      dup.id <- names(seqStat)[which(seqStat>1)]
      if(rm.dup)
      {
        warning( paste( "Identical loci are removed from the BED data. (e.g. ", paste(head(dup.id), collapse =",", se=""), ")" ) );
        seqList <- unique(seqList);
      }
      else
        warning( paste( "Identical loci are found in the BED data. (e.g. ", paste(head(dup.id), collapse =",", se=""), ")" ) );
    }
  }

  writeLines(seqList, tmp.seq);

  # generate fasta file
  cmd = paste("twoBitToFa -seqList=", tmp.seq, " ", twoBitPath, " ", tmp.fa, sep="")
  err_code <- system(cmd, wait = TRUE);
  if( err_code != 0 )
  {
    cat("Failed to call twoBitToFa to get the sequence data.\n");
    unlink(c(tmp.seq, tmp.fa))
    return(NULL);
  }
  else
  {
    # read data
    ms_data <- read.ms(tmp.fa)
    # clean up
    unlink(c(tmp.seq, tmp.fa))
    return(ms_data)
  }
}


tfbs_to_bed <- function(sites)
{
  ## Alternative formulation fails if sites is empty.
  if (NROW(sites) == 0)
    return(NULL)

  spl <- strsplit(as.character(sites$seqname), ":|-")

  chroms <- sapply(spl, function(pair) pair[1]);
  starts <- as.integer(sapply(spl, function(pair) pair[2]));
  ends   <- as.integer(sapply(spl, function(pair) pair[3]));

  bed = data.frame(
          tre.chrom     = chroms,
          tre.chromStart= as.integer(starts),
          tre.chromEnd  = as.integer(ends),
          tf.chrom      = chroms,
          tf.chromStart = as.integer(starts + sites$start - 1), # add site offset
          tf.chromEnd   = as.integer(starts + sites$end),
          score      = sites$score,
          strand     = sites$strand )

  return(bed)
}


locate.TF<-function(positive.bed, negative.bed, motif.id, half.size, mTH, ncores, tfbs, file.twoBit,return.type="list", background.order=2){

  stopifnot(class(tfbs) == "tfbs");

  # read sequences
  positive.ms = read.seqfile.from.bed(positive.bed, file.twoBit);
  negative.ms =read.seqfile.from.bed(negative.bed, file.twoBit);

  motif.idx <- match(motif.id, tfbs@ mgisymbols);
 
  #compute background models
  both.ms = concat.ms(positive.ms, negative.ms)
  background.mm = build.mm(both.ms, background.order)

  # iterate over TF set
  r.comp <- mclapply(motif.idx, function(i, ...) {
    # get PWM information
    pwm = tfbs@pwm[[i]];

    motif.id<-tfbs@mgisymbols[i];
    motif.name<-as.character(tfbs@tf_info$TF_Name[i]);

    pos.sites = score.ms(positive.ms, pwm, background.mm, threshold = mTH);

    result.bed = tfbs_to_bed(pos.sites);
    result.bed$motif.name <- motif.name;
    result.bed$motif.id <- motif.id;
    result.bed$motif.idx<-i;
    return(result.bed);
  }, mc.cores = ncores)

  if(return.type=="list")
  	return(r.comp)
  else
    return(do.call(rbind.data.frame, r.comp))

}


bedTools.merge<-function(bed)
{
  #create temp files
  a.file=tempfile()
  options(scipen =99) # not to use scientific notation when writing out

  #write bed formatted dataframes to tempfile
  write.table(bed,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)

  # create the command string and call the command using system()
  command <- paste("LC_COLLATE=C sort -k1,1 -k2,2n",a.file,"| mergeBed -s -c 6 -o distinct -i stdin ",sep=" ")

  res    <- read.table(pipe(command),header=F);
  res$V5 <- NA;
  res$V6 <- res[,4];
  res$V4 <- NA;

  unlink(a.file);

  return(res);
}

#function that calls bedtools and operate on two bed dataframes
bedTools.2in<-function(bed1,bed2)
{
  #create temp files
  a.file=tempfile();
  b.file=tempfile();
  options(scipen =99); # not to use scientific notation when writing out

  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F);
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F);

  # create the command string and call the command using system()
  command <- paste("bedtools intersect -s -c -a",a.file,"-b",b.file,sep=" ");

  #res=as.numeric(read.table(pipe(command),header=F)[,1])>0
  res <- read.table(pipe(command),header=F)
  unlink(a.file);
  unlink(b.file);

  #return(as.numeric(res))
  return(as.numeric(res[,ncol(res)]>0))
}

get.pos.mat<-function(TF.tab.list, ncores){

  TF.tab.df <- do.call(rbind.data.frame, TF.tab.list)

  TF.tab.df <- TF.tab.df[,c("tf.chrom","tf.chromStart","tf.chromEnd","motif.name","motif.id","strand")];
  TF.tab.df.merged <- bedTools.merge(TF.tab.df);
  TF.tab.df.merged.small <- TF.tab.df.merged[TF.tab.df.merged$V3-TF.tab.df.merged$V2<=30,];

  pos.mat <- do.call(cbind.data.frame,
     mclapply(TF.tab.list,
       FUN=function(x) {
		 bedTools.2in(bed1 = TF.tab.df.merged.small, bed2 = x[,c("tf.chrom","tf.chromStart","tf.chromEnd","motif.name","motif.id","strand")])}
	    ,mc.cores= ncores));

  return(pos.mat)
}

cluster.motif.pos<-function(motif.df, changed.bed, unchanged.bed, half.size, mTH, ncores, tfs, file.twoBit, pdf.name=NULL){

  query.motifs <- motif.df[,1];
  if(length(query.motifs)<=1)
    return(NULL)

  TF.tab <- locate.TF( changed.bed, unchanged.bed, query.motifs, half.size, mTH, ncores, tfs, file.twoBit);

  if(is.null(pdf.name))
    return(do.call(rbind.data.frame,TF.tab))

  pos.mat <- get.pos.mat(TF.tab, ncores);
  colnames(pos.mat)<-paste(motif.df[,1], motif.df[,2],sep="_");

  cor.mat <- cor(pos.mat,method="spearman");

  pdf( paste(pdf.name, ".cor.heatmap.pdf",sep=""), pointsize=8, useDingbats=FALSE, paper="letter");

  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

  
  hc <- hclust(as.dist(1-cor.mat), method="ward.D2");
  hist.order <- heatmap.2(as.matrix(cor.mat),
       col = my_palette,
       density.info = "none",
       trace ="none",
       Rowv = as.dendrogram(hc),
       Colv = as.dendrogram(hc),
       dendrogram = c("both"),
       symm = T,
       symkey = F,
       symbreaks = T,
       scale = "none",
       margins=c(12,12),
       key.xlab="correlation coefficient",
       key.title="",
       cexRow = 0.2 + 1.2/log10(nrow(cor.mat)+1),
       cexCol = 0.2 + 1.2/log10(ncol(cor.mat)+1));

  dev.off();
  motif.df.ordered<-motif.df[rev(hist.order$rowInd),];

  pdf.name<-paste(pdf.name,".motif.pdf",sep="");
  plot_motif_table(motif.df.ordered, tfs=tfs, pdf.name= pdf.name);

  return(NULL);

}


plot.tfTarget<-function( tfTar, out.prefix,  merge.motif=T){

  stopifnot(!is.null(tfTar$df.motif.up.all) & !is.null(tfTar$df.motif.down.all));
  df.motif.up <- tfTar$df.motif.up.all;
  df.motif.down <- tfTar$df.motif.down.all;

  if(merge.motif) {
    df.motif.up <- merge.motif.df.by.name(df.motif.up);
    df.motif.down <- merge.motif.df.by.name(df.motif.down);
  }

  pdf.name <- paste( out.prefix,
                    "pval.upBd", tfTar$pval.cutoff.up,
                    "pval.lowBd", tfTar$pval.cutoff.down,
                    "mTH", tfTar$mTH,
                    "rep", tfTar$run.repeats,
                    "fdr", tfTar$fdr.cutoff,
                    "up",sep="_");

  cluster.motif.pos(df.motif.up,
                    tfTar$enh.up.bed,
                    tfTar$enh.unc.bed,
                    tfTar$half.size,
                    tfTar$mTH,
                    tfTar$ncores,
                    tfTar$tfs,
                    tfTar$file.twoBit,
                    pdf.name);

  pdf.name <- paste( out.prefix,
                    "pval.upBd", tfTar$pval.cutoff.up,
                    "pval.lowBd", tfTar$pval.cutoff.down,
                    "mTH", tfTar$mTH,
                    "rep", tfTar$run.repeats,
                    "fdr", tfTar$fdr.cutoff,
                    "down",sep="_");

  cluster.motif.pos( df.motif.down,
                     tfTar$enh.down.bed,
                     tfTar$enh.unc.bed,
                     tfTar$half.size,
                     tfTar$mTH,
                     tfTar$ncores,
                     tfTar$tfs,
                     tfTar$file.twoBit,
                     pdf.name);

  invisible()

}
