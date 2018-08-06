
get_robust_motif.median.pval.sites <-function(motif.list, fdr.cutoff=0.02, sites.num.cutoff =10){

  fe.ratios<-do.call(cbind,lapply(motif.list,FUN=function(x) x$fe.ratio))
  fe.ratios<-apply(fe.ratios ,FUN=median,MARGIN=1,na.rm =T)

  pval.adj<-do.call(cbind,lapply(motif.list,FUN=function(x) x$pvalue))
   pval.adj <-apply(pval.adj ,FUN=median,MARGIN=1,na.rm =T)

   sites.num<-do.call(cbind,lapply(motif.list,FUN=function(x) x$Npos))
   sites.num <-apply(sites.num ,FUN=median,MARGIN=1,na.rm =T)

  return(which(pval.adj <fdr.cutoff & fe.ratios>1 & sites.num> sites.num.cutoff))
}


merge.motif.list<-function(motif.list, motifs.to.plot){
  fe.ratios<-do.call(cbind,lapply(motif.list,FUN=function(x) x[motifs.to.plot,5]))
  fe.ratios<-apply(fe.ratios ,FUN=median,MARGIN=1,na.rm =T)
  pval.adj<-do.call(cbind,lapply(motif.list,FUN=function(x) x[motifs.to.plot,6]))
  pval.adj <-apply(pval.adj ,FUN=median,MARGIN=1,na.rm =T)

  motif.ids<-motif.list[[1]][motifs.to.plot,1]
  motif.names<-motif.list[[1]][motifs.to.plot,2]
  motif.df<- cbind.data.frame(motif.ids, motif.names, fe.pval=paste(fe.ratios, pval.adj,sep=","), motif.ids)

}


merge.motif.df.by.name<-function(motif.df){

  get.fe.pval.prod<-function(x){unlist(lapply(x, FUN=function(i){
    fe.pval<- as.numeric(unlist(strsplit(i,split=",")))
    log2.fe <- log2(fe.pval[1])
    if(log2.fe<0) log2.fe<-0
    pval<-fe.pval[2]
    if(pval==0) pval<-0.001
    miunus.log10.pval<-  -log10(pval)
    return(log2.fe*miunus.log10.pval)
  }))}

  merge.by.name<-function(motif.name, motif.df){
    motif.df.tf<-motif.df[motif.df$motif.names== motif.name,]
    if(nrow(motif.df.tf)==1) return(motif.df.tf)

    motif.df.dat<- as.character(motif.df.tf[,3])
    motif.df.dat.prod<-get.fe.pval.prod(motif.df.dat)
    return(motif.df.tf[which.max(motif.df.dat.prod),])
  }

  all.tf.names<-unique(as.character(motif.df $ motif.names))
  do.call(rbind.data.frame, lapply(all.tf.names, FUN= merge.by.name, motif.df= motif.df))
}


get.motif.df<-function(motif.list, fdr.cutoff, sites.num.cutoff, exp.cutoff, plus.files.query, minus.files.query, plus.files.control, minus.files.control, bigWig.path, tf.gene.path){

  motifs.to.plot <-get_robust_motif.median.pval.sites(motif.list, fdr.cutoff, sites.num.cutoff)
  df.motif <-merge.motif.list(motif.list, motifs.to.plot)

  tf.exp.mat<-get.raw.count(bigWig.path= bigWig.path, plus.files=c(plus.files.query, plus.files.control), minus.files=c(minus.files.query, minus.files.control), tf.names= df.motif[,2], tf.gene.path)
  
  if.exp.query<-apply(tf.exp.mat[,c(1:length(plus.files.query))]> exp.cutoff, 1, prod)
  if.exp.control<-apply(tf.exp.mat[,-c(1:length(plus.files.query))]> exp.cutoff, 1, prod)
  
  if.exp<- if.exp.query+if.exp.control

  df.motif.exp<-df.motif[!is.na(if.exp) & if.exp>0,]

  return(df.motif.exp)
}


filter.rtfbsdb<-function(tfTar, fdr.cutoff=0.05, sites.num.cutoff=10, exp.cutoff=2, ncores=1){

  if(class(tfTar)!="tfTarget")
     stop("The first parameter is not a 'tfTarget' object!");

  if(class(tfTar$tfs)!="tfbs")
     stop("No 'tfbs' object in tfTar!");

  options(scipen =99);

  df.motif.up.all <- get.motif.df( motif.list = tfTar$motif.list.up,
                                   fdr.cutoff = fdr.cutoff,
                                   sites.num.cutoff = sites.num.cutoff,
                                   exp.cutoff = exp.cutoff,
                                   plus.files.query = tfTar$plus.files.query,
                                   minus.files.query = tfTar$minus.files.query,
                                   plus.files.control = tfTar$plus.files.control,
                                   minus.files.control = tfTar$minus.files.control,
                                   bigWig.path = tfTar$bigWig.path,
                                   tf.gene.path = tfTar$gene.path);

  TF.TRE.tab.up <- cluster.motif.pos( df.motif.up.all,
                                    tfTar$enh.up.bed,
                                    tfTar$enh.unc.bed,
                                    tfTar$half.size,
                                    tfTar$mTH,
                                    ncores,
                                    tfTar$tfs,
                                    tfTar$file.twoBit,NULL);

  df.motif.down.all <- get.motif.df(motif.list = tfTar$motif.list.down,
                                    fdr.cutoff = fdr.cutoff,
                                    sites.num.cutoff = sites.num.cutoff,
                                    exp.cutoff = exp.cutoff,
                                    plus.files.query = tfTar$plus.files.query,
                                    minus.files.query = tfTar$minus.files.query,
                                    plus.files.control = tfTar$plus.files.control,
                                    minus.files.control = tfTar$minus.files.control,                                    
                                    bigWig.path = tfTar$bigWig.path,
                                    tf.gene.path = tfTar$gene.path)

  TF.TRE.tab.down<-cluster.motif.pos( df.motif.down.all,
                                     tfTar$enh.down.bed,
                                     tfTar$enh.unc.bed,
                                     tfTar$half.size,
                                     tfTar$mTH,
                                     ncores,
                                     tfTar$tfs,
                                     tfTar$file.twoBit,NULL)

  TF.TRE.tab <- rbind.data.frame( TF.TRE.tab.up, TF.TRE.tab.down );

  tfTar$TF.TRE.tab <- TF.TRE.tab;
  tfTar$df.motif.up.all <- df.motif.up.all;
  tfTar$df.motif.down.all <- df.motif.down.all;
  tfTar$TF.TRE.tab <- TF.TRE.tab;

  tfTar$fdr.cutoff <- fdr.cutoff;
  tfTar$sites.num.cutoff <- sites.num.cutoff;
  tfTar$exp.cutoff <- exp.cutoff;

  return(tfTar);
}
