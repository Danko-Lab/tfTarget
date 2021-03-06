
getCounts.gene.rpkm <- function(plus, minus, intervals= bodies) {
  pl <- load.bigWig(plus)
  mn <- load.bigWig(minus)
  counts <- bed6.region.bpQuery.bigWig(pl, mn, intervals, abs.value = TRUE) ## Get counts
  counts<- counts/ ((pl$basesCovered*pl$mean+abs(mn$basesCovered*mn$mean))/10^6)
  counts<- counts/((intervals[,3]-intervals[,2])/10^3)
  unload.bigWig(pl);
  unload.bigWig(mn);
}

getCounts.gene <- function(plus, minus, intervals= bodies) {
  pl <- load.bigWig(plus)
  mn <- load.bigWig(minus)
  counts <- bed6.region.bpQuery.bigWig(pl, mn, intervals, abs.value = TRUE) ## Get counts
  unload.bigWig(pl);
  unload.bigWig(mn);
}

get.raw.count <- function(bigWig.path, plus.files, minus.files, tf.names, tf.gene.path){

  stopifnot(length(plus.files)==length(minus.files));

  #739 unique TF names, 726 found in gencode.v25lift37.annotation.gene.gtf.sorted.bed #724 with length >1000
  #TF does not pass length filter: FERD3L=639 and SRY=844
  #TF gene length distribution
  #Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 639    5814   18544   60438   59376 1216867

  #multiple position TFs: (chrY filtered)
  # chrX 585579 620146  N  N  + SHOX protein_coding
  # chrY 535579 570146  N  N  + SHOX protein_coding

  # chrX 2404455 2418508  N  N  - ZBED1 protein_coding
  # chrY 2354455 2368508  N  N  - ZBED1 protein_coding

  refGene <- read.table(tf.gene.path);
  refGene <- refGene[grep("random|Un|hap", refGene$V1, invert=TRUE),];
  refGene <-refGene[refGene[,1] %in% unique(read.table(TRE.path)[,1]),]
  refGene <- refGene[(refGene$V3-refGene$V2)>600,];

  bodies <- refGene;
  bodies$V2[bodies$V6 == "+"] <-bodies$V2[bodies$V6 == "+"]+500;
  bodies$V3[bodies$V6 == "-"] <- bodies$V3[bodies$V6 == "-"]-500;
  bodies$V4<-"N";
  bodies <-unique(bodies);

  file.names <- gsub(".bw","", gsub(".bigWig","", plus.files));
  file.names <- gsub(".pl","", gsub(".plus","", file.names));
  file.names <- gsub("_pl","", gsub("_plus","", file.names));

  plus.files <- paste(bigWig.path, plus.files, sep="/");
  minus.files<- paste(bigWig.path, minus.files, sep="/");

  raw_counts <- do.call(cbind,mclapply(1:length(plus.files),
      function(i) getCounts.gene(plus.files[i], minus.files[i],intervals= bodies),
    mc.cores=length(plus.files)) );
  raw_counts<-as.matrix(raw_counts,drop=F);

  colnames(raw_counts) <-file.names;
  rownames(raw_counts) <- bodies$V5;
  raw_counts <-raw_counts[match(tf.names ,rownames(raw_counts) ), ,drop=F];

  return(raw_counts);
}