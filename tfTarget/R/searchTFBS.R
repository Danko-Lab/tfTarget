background.check<-function( gc.pos, gc.neg, gc.correction, file.pdf.vioplot=NA, gc.min.sample = 500, verbose=TRUE )
{
  print("new function detected")
  pdf.output <- FALSE;
  gc.test <- wilcox.test(gc.pos, gc.neg, conf.int=TRUE, conf.level=0.9 );

  # Actually 0.01 is not good for wilcox.test,
  # We use the customizable pvalue in this function.
  if( verbose)
  {
    cat("! Difference between GC content in negative and positive TREs (use 'gc.correction.pdf' to see pdf figure):\n");
    cat("  p-value (Wilcoxon-Mann-Whitney test) =", gc.test$p.value, "\n");
    cat("  median/sample size of GC positive =", median(gc.pos), "/", length(gc.pos), "\n");
    cat("  median/sample size of GC negative =", median(gc.neg), "/", length(gc.neg), "\n");

    if( !is.null(file.pdf.vioplot) && !is.na( file.pdf.vioplot)  )
    {
      r.try <- try ( pdf(file.pdf.vioplot) );
      if( class(r.try) != "try-error" )
      {
        vioplot( gc.pos, gc.neg, names=c("Positive", "Negative") );
        abline(h=median(gc.pos), lty="dotted");
        pdf.output <- TRUE;
      }
      else
        cat("  Failed to output the vioplot figure for the results.\n");
    }
  }

  if( gc.correction==FALSE )
  {
    # PDF has not been finished until dev.off();
    if( pdf.output )
    {
      dev.off();
      cat("* Please check the vioplot figure to make sure, the vioplot figure: ", file.pdf.vioplot, "\n" );
    }
    return(NULL);
  }

  if( length(gc.neg)< 5000 )
  {
    if(verbose) cat("! Failed to make GC correction due to small bed data(size<5000).\n");
    return(NULL);
  }

  resample <- function(ref, orig, n=10000, nbins=10) {
    ## (1) Break gcContent down into 10 equally sized bins.
    breaks <- seq(0, 1, length.out=nbins) #seq(min(ref), max(ref), length.out=nbins)

    ## (2) Get the empirical frequencies of each bin.
    empir <- sapply(1:(NROW(breaks)-1), function(x) {sum(breaks[x]<= ref & ref < breaks[x+1])/NROW(ref)})

    ## (3) Re-sample TREs in the BG set w/ probability proportional to the bin.
    resamp_prob <- rep(1, NROW(orig))
    for(i in 1:NROW(empir)) {
      incl <- breaks[i] <= orig & orig < breaks[i+1]
      resamp_prob[incl] <- empir[i]/ sum(incl)
    }

    sample(1:NROW(resamp_prob), n, prob= resamp_prob, replace=FALSE)
  }

  try.sample<-function(n.sample)
  {
    indx.bgnew <- resample( gc.pos, gc.neg, n=n.sample)
    gc.testx <- wilcox.test(gc.pos, gc.neg[indx.bgnew], conf.int=TRUE, conf.level=0.9 );
    return(gc.testx$p.value);
  }


  ns.sample <- c ( round(length(gc.neg)/c( 2:10,15,20,25)), 1000, 800, 500 );
  ns.sample <- ns.sample[ ns.sample >= gc.min.sample  ];

  ns.pvalue <- unlist(lapply(ns.sample, try.sample));
  n.sample  <- ns.sample[ which.max(ns.pvalue) ];

  indx.bgnew <- resample( gc.pos, gc.neg, n=n.sample)
  gc.test2 <- wilcox.test(gc.pos, gc.neg[indx.bgnew], conf.int=TRUE, conf.level=0.9 );

  if(verbose)
  {
    cat("* After the resampling from negative TREs:\n");
    cat("  p-value (Wilcoxon-Mann-Whitney test) =", gc.test2$p.value, "\n");
    cat("  median/sample size of GC negative =", median(gc.neg[indx.bgnew]), "/", length(indx.bgnew), "\n");
  }

  if( pdf.output )
  {
    vioplot(gc.pos, gc.neg, gc.neg[indx.bgnew], names=c("Positive", "Negative", "Negative.resample"));
    abline(h=median(gc.pos), lty="dotted")
    dev.off();
    cat("* Please check the vioplot figure to make sure, the vioplot figure: ", file.pdf.vioplot, "\n" );
  }

  ## return sampling background.
  return( indx.bgnew );
}


tfbs_enrichmentTest<-function( tfbs, file.genome,
          positive.bed,
          negative.bed = NA,
          file.prefix = NA,
          use.cluster = FALSE,
          ncores = 1,
          gc.correction = TRUE,
          gc.correction.pdf = NA,
          gc.min.sample = 500,
          gc.robust.rep = NA,
          threshold = 6,
          threshold.type = c("score", "fdr"),
          gc.groups = 1,
          background.order = 2,
          background.length = 100000,
          pv.adj = p.adjust.methods)
{
  stopifnot(class(tfbs) == "tfbs")

  if( !check_bed( positive.bed ) )
    stop("Wrong format in the parameter of 'positive.bed', at least three columns including chromosome, strat, stop.");

  if( !missing(negative.bed) && !check_bed( negative.bed ) )
    stop("Wrong format in the parameter of 'negative.bed', at least three columns including chromosome, strat, stop.");

  if( missing( threshold.type ) ) threshold.type <- "score";
  if( threshold.type == "score" && missing( threshold ) ) threshold <- 6;
  if( threshold.type == "fdr" && missing( threshold ) ) threshold <- 0.1;

  if( !missing(pv.adj)) pv.adj <- match.arg(pv.adj)
  if( missing(pv.adj) ) pv.adj <- "bonferroni";
  if( pv.adj == "fdr" ) pv.adj <- "BH";

  if( missing(gc.robust.rep) || is.na(gc.robust.rep)) gc.robust.rep <- 1;
  if( gc.robust.rep >1 && gc.robust.rep<3 ) gc.robust.rep <- 3;

  if( missing( gc.min.sample ) ) gc.min.sample<- 500;
  if( missing( ncores) ) ncores <- 1;
  if( missing( gc.groups) ) gc.groups <- 1;
  if( missing( gc.correction) ) gc.correction <- FALSE;
  if( missing( background.order ) ) background.order <- 2;
  if( missing( background.length ) ) background.length <- 100000;
  if( missing( use.cluster) ) use.cluster <- FALSE;
  if( use.cluster && NROW(tfbs@cluster)==0 )
    stop("No cluster information in the tfbs object");
  if( use.cluster && NCOL(tfbs@cluster)!=3 )
    stop("No selected motif for each cluster in the tfbs object");

  file.twoBit = file.genome;
  if( tolower( file_ext( file.genome ) ) != "2bit" )
  {
    file.twoBit = tempfile(fileext=".2bit")

    # generate fasta file
    err_code <- system(paste("faToTwoBit ", file.genome, " ", file.twoBit), wait = TRUE);
    if( err_code != 0 || !file.exists (file.twoBit) )
      stop("Failed to call faToTwoBit to convert FASTFA file.");
  }

  if( use.cluster )
  {
    r.mat <- range( tfbs@cluster[,1] );
    if( r.mat[1]<1 || r.mat[2]>tfbs@ntfs )
      stop("The first column of 'tfbs@cluster' exceeds the range of motif data set.");

    cluster.mat <- tfbs@cluster[ tfbs@cluster[,3]==1, c(1,2), drop=F];
  }
  else
    cluster.mat <- cbind( 1:tfbs@ntfs, 1);

  if( missing(negative.bed) )
  {
    negative.bed <- background.generate( positive.bed );
    cat("*", NROW(negative.bed),  "GC negative loci are randomly generated.\n");
  }

  # read sequences
  positive.ms = read.seqfile.from.bed(positive.bed, file.twoBit)
  negative.ms = read.seqfile.from.bed(negative.bed, file.twoBit)

  r.comp <- comparative_scanDb_rtfbs( tfbs,
              file.twoBit,
              positive.bed,
              negative.bed,
              positive.ms,
              negative.ms,
              file.prefix,
              cluster.mat = cluster.mat,
              ncores = ncores,
              gc.correction = gc.correction,
              gc.correction.pdf = gc.correction.pdf,
              gc.min.sample = gc.min.sample,
              gc.correction.verbose = TRUE,
              threshold = threshold ,
              threshold.type = threshold.type,
              gc.groups = gc.groups,
              background.order = background.order,
              background.length = background.length,
              pv.adj=pv.adj );

  if( is.null(r.comp) ) return(NULL);
  ret <- r.comp$ret;

  if( !is.null(r.comp$bg.sample)  && gc.robust.rep > 1 )
  {
    ret.list <- list();
    ret.list[[1]] <- ret;

    for(i in 2:gc.robust.rep)
    {
      cat("* GC robust replication for background resampling, loop=", i, "\n");
      r.comp <- comparative_scanDb_rtfbs( tfbs,
                file.twoBit,
                positive.bed,
                negative.bed,
                positive.ms,
                negative.ms,
                file.prefix,
                cluster.mat = cluster.mat,
                ncores = ncores,
                gc.correction = gc.correction,
                gc.correction.pdf = gc.correction.pdf,
                gc.min.sample = gc.min.sample,
                gc.correction.verbose = FALSE,
                threshold = threshold ,
                threshold.type = threshold.type,
                gc.groups = gc.groups,
                background.order = background.order,
                background.length = background.length,
                pv.adj=pv.adj );
      ret.list[[i]] <- r.comp$ret;
    }

    for(i in 1:NROW(ret))
    {
      Npos.list <- c();
      Nneg.list <- c();

      for(k in 1:gc.robust.rep)
      {
        Npos.list <- c( Npos.list, ret.list[[k]][i,'Npos'] );
        Nneg.list <- c( Nneg.list, ret.list[[k]][i,'Nneg'] );
      }

      Npos <- median( Npos.list );
      Nneg <- median( Nneg.list );
      n.gc.pos <- ret.list[[1]][i,'gc.pos'];
      n.gc.neg <- ret.list[[1]][i,'gc.neg'];
      Nneg.expected <- Nneg * n.gc.pos/ n.gc.neg ;

      Npos <- ifelse( Npos==0, 1, Npos );
      Nneg <- ifelse( Nneg==0, 1, Nneg );

      fe.ratio <- ( Npos / n.gc.pos )/( Nneg / n.gc.neg );

      pos.left <-  n.gc.pos - Npos;
      if( pos.left<0 ) pos.left  <- 0;
      neg.left <- n.gc.neg - Nneg;
      if( neg.left<0 ) neg.left  <- 0;

      tbl  = rbind( c( Npos, Nneg ), c( pos.left, neg.left ) )
      pval = fisher.test(tbl)$p.value;

      ret[i, 'Npos']     <- Npos;
      ret[i, 'Nneg']     <- Nneg;
      ret[i, 'expected'] <- round( Nneg.expected, 1);
      ret[i, 'fe.ratio'] <- fe.ratio;
      ret[i, 'pvalue']   <- pval;
    }

    ret$pv.adj <- adjust.pvale( ret$pvalue, cluster.mat, pv.adj );
  }

  ret$Nneg   <- NULL;
  ret$gc.pos <- NULL;
  ret$gc.neg <- NULL;

  if(!is.null(ret))
  {
    if( NROW(tfbs@tf_info) > 0 )
    {
      tf.motifid   <- unlist(strsplit(as.character( ret$tf.name ), "@"))[seq(1, length(ret$tf.name)*2-1, 2)];
      tf.idx       <- match( tf.motifid, tfbs@tf_info$Motif_ID );
      ret$tf.name  <- tfbs@tf_info$TF_Name[ tf.idx ];
      ret$motif.id <- tfbs@tf_info$Motif_ID[ tf.idx ];
    }
  }

  r.parm <- list( file.genome   = file.genome,
        file.prefix       = file.prefix,
        use.cluster       = use.cluster,
        cluster.mat       = cluster.mat,
        ncores            = ncores,
        gc.robust.rep     = gc.robust.rep,
        threshold.type    = threshold.type,
        threshold         = threshold,
        pv.adj            = pv.adj,
        background.order  = background.order,
        background.length = background.length,
        gc.min.sample     = gc.min.sample,
        gc.correction     = gc.correction );

  r <- list( result = ret, parm = r.parm);
  class(r) <- c("tfbs.enrichment");

  return(r);
}

# enh.up.bed: upregulated bed file, in dataframe
# enh.unc.bed: control bed file, in dataframe
# half.size: bp from the center, default= 150
# mTH: threshold of motif score, default=10
# min.size: min of sampling size, default=2500
# run.repeats: number of runs for robust test, defult=50


#do not use cluster
tfbs.enrichmentTest.multiple<-function( tfbs, file.twoBit, enh.up.bed, enh.unc.bed, mTH, min.size, run.repeats, ncores ){

  print("repeating how many times")
  print(run.repeats)

  print(paste("TREs changed/background: ", NROW(enh.up.bed), NROW(enh.unc.bed)))

  motifs_list<-list()

  for (i in 1: run.repeats){
    print(i)
    set.seed(i)
    t.comp <- tfbs.enrichmentTest(
          tfbs,
          file.twoBit,
          enh.up.bed,
          negative.bed= enh.unc.bed,
          gc.correction=TRUE,
          gc.min.sample= min.size,
          threshold = mTH,
          pv.adj = "fdr",
          ncores = ncores,
          use.cluster=FALSE);
      res<-t.comp$result
      motifs_list[[i]]<-res
    }

  return(motifs_list)

}

searchTFBS <- function(tfTar, tfs, file.twoBit, pval.cutoff.up=0.01, pval.cutoff.down=0.1, half.size=150, mTH=7, min.size=150, run.repeats=2, ncores=1 ){
  if(class(tfTar)!="tfTarget")
     stop("The first parameter is not a 'tfTarget' object!");

  if(class(tfs)!="tfbs")
     stop("The second parameter is not a 'tfbs' object!");

  if(!all(file.exists(file.twoBit)))
     stop("The 2bit file is not found!");

  options(scipen =99);

  deseq.table.TRE <-tfTar$deseq.table.TRE;
  #ncores <- tfTar$ncores;

  deseq.table.sig <- center.bed(deseq.table.TRE[!is.na(deseq.table.TRE$TRE.padj) & deseq.table.TRE$TRE.padj <pval.cutoff.up,], half.size, half.size)
  enh.unc.bed     <- center.bed(deseq.table.TRE[!is.na(deseq.table.TRE$TRE.padj) & deseq.table.TRE$TRE.padj> pval.cutoff.down,], half.size, half.size)
  enh.up.bed      <- deseq.table.sig[deseq.table.sig$TRE.log2FoldChange>0,]
  enh.down.bed    <- deseq.table.sig[deseq.table.sig$TRE.log2FoldChange<0,]

  motif.list.up <- tfbs.enrichmentTest.multiple(tfs, file.twoBit,  enh.up.bed, enh.unc.bed, mTH, min.size, run.repeats, ncores);
  motif.list.down <- tfbs.enrichmentTest.multiple(tfs, file.twoBit, enh.down.bed, enh.unc.bed, mTH, min.size, run.repeats, ncores);

  tfTar$pval.cutoff.up <- pval.cutoff.up;
  tfTar$pval.cutoff.down <- pval.cutoff.down;
  tfTar$half.size <- half.size;
  tfTar$mTH <- mTH;
  tfTar$min.size <- min.size;
  tfTar$run.repeats <- run.repeats;
  tfTar$tfs <- tfs;
  tfTar$file.twoBit <- file.twoBit;

  tfTar$enh.up.bed <- enh.up.bed;
  tfTar$enh.down.bed <- enh.down.bed;
  tfTar$enh.unc.bed <- enh.unc.bed;
  tfTar$motif.list.up <- motif.list.up;
  tfTar$motif.list.down <- motif.list.down;

  return(tfTar);

}