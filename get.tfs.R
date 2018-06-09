
library("rtfbsdb")

args<-commandArgs(TRUE)

db <- CisBP.extdata(args);
tfs <- tfbs.createFromCisBP(db);

nonredundant.idx <-!duplicated(tfs@tf_info[,c("Motif_ID","TF_Name")])

tfs@ntfs<-as.integer(table(nonredundant.idx)["TRUE"])
tfs@mgisymbols<-tfs@mgisymbols[nonredundant.idx]
tfs@filename<-tfs@filename[nonredundant.idx]
tfs@pwm<-tfs@pwm[nonredundant.idx]
tfs@tf_info<-tfs@tf_info[nonredundant.idx,]


tfs@cluster<-cbind(tfs@cluster,1)
attr(tfs@cluster,"dimnames")[[2]]<-c("subset", "clusters" ,"selected")

save(file=paste(args,".tfs.rdata",sep=""))










