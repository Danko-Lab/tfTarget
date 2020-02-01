library("rtfbsdb")

## format R --args species cisbp.zip.file

args<-commandArgs(TRUE)

if(args[1] %in% c("Homo_sapiens", "Mus_musculus")) {
   db <- CisBP.extdata(args)} else 
{
   if (NROW(args)==1)
     db <- CisBP.download(args[1])
   else
   {
      if (file.exists(args[2]))
         db <- CisBP.zipload(args[2], species=args[1])
      else
         stop(paste(args[2], "is not existing for", args[1], ".\n"))   
   }      
}

tfs <- tfbs.createFromCisBP(db)

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


save(tfs, file=paste(args[1],".tfs.rdata",sep=""))
cat("TFS data is stored into", paste(args[1],".tfs.rdata",sep=""), ".\n");






