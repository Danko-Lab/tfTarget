summary.tfTarget <- function( tfTar ){

  if(!is.null( tfTar$enh.up.bed ) & !is.null( tfTar$enh.down.bed ) & !is.null( tfTar$enh.unc.bed )){
    TRE.up.num <- nrow( tfTar$enh.up.bed )
    TRE.down.num <- nrow( tfTar$enh.down.bed )
    TRE.unc.num <- nrow( tfTar$enh.unc.bed )
  }
  else
    return(NULL)

  if(!is.null( tfTar$df.motif.up.all ) & !is.null( tfTar$df.motif.down.all )){
    up.tf.num <- length( unique(tfTar$df.motif.up.all$motif.names) )
    down.tf.num <- length( unique(tfTar$df.motif.down.all$motif.names) )
    res.df <- data.frame( count=c(TRE.up.num, TRE.down.num, TRE.unc.num, up.tf.num, down.tf.num))
    rownames(res.df) <- c( "TRE.up.num", "TRE.down.num", "TRE.unc.num", "up.tf.num", "down.tf.num")
    return(res.df)
  }

  else{
    res.df <- data.frame( count=c(TRE.up.num, TRE.down.num, TRE.unc.num))
    rownames(res.df) <- c("TRE.up.num", "TRE.down.num", "TRE.unc.num")
    return(res.df)
  }
}

print.tfTarget<-function(tfTar){

  if(!is.null( tfTar$enh.up.bed ) & !is.null( tfTar$enh.down.bed ) & !is.null( tfTar$enh.unc.bed )){
    TRE.up.num <- nrow(tfTar$enh.up.bed)
    TRE.down.num <- nrow(tfTar$enh.down.bed)
    TRE.unc.num <- nrow(tfTar$enh.unc.bed)

    cat("number of up-regulated TREs=", TRE.up.num,"\n")
    cat("number of down-regulated TREs=", TRE.down.num,"\n")
    cat("number of unc-regulated TREs=", TRE.unc.num,"\n")
  }
  else
    if(!is.null( tfTar$df.motif.up.all ) & !is.null( tfTar$df.motif.down.all )){
      up.tf.num <- length( unique(tfTar$df.motif.up.all$motif.names) )
      down.tf.num <- length( unique(tfTar$df.motif.down.all$motif.names) )

      cat( "number of TFs enriched in up-regulated TREs=", up.tf.num,"\n")
      print( as.character(unique(tfTar$df.motif.up.all$motif.names)))
      cat( "number of TFs enriched in down-regulated TREs=", down.tf.num,"\n")
      print( as.character(unique(tfTar$df.motif.down.all$motif.names)))
  }

  invisible()
}