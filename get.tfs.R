
library("rtfbsdb")

args<-commandArgs(TRUE)

db <- CisBP.extdata(args);
tfs <- tfbs.createFromCisBP(db);


#merge redundant motif IDs

tf_info.ori<-tfs@tf_info
motif.ids.uniq<-unique(as.character(tfs@tf_info[,c("Motif_ID")]))

merge.id<-function(motif.id){
	ret.df<-tfs@tf_info[tfs@tf_info$Motif_ID== motif.id,,drop=F]
	unlist(apply(ret.df,2,function(input.col) paste(unique(input.col),collapse="/")) )
}

tf_info.uniq<-unique(do.call(rbind.data.frame,lapply(motif.ids.uniq, merge.id)))
colnames(tf_info.uniq)<-colnames(tfs@tf_info)

tfs@tf_info<-tf_info.uniq

tfs@ntfs<-nrow(tf_info.uniq)
tfs@mgisymbols<-as.character(tf_info.uniq$Motif_ID)
tfs@filename<-tfs@filename[match(tf_info.uniq$Motif_ID,tf_info.ori$Motif_ID)]
tfs@pwm<-tfs@pwm[match(tf_info.uniq$Motif_ID,tf_info.ori$Motif_ID)]


save(tfs, file=paste(args,".tfs.rdata",sep=""))







