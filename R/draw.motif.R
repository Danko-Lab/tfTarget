
plot_motif_info <-function(tfbs, motif.df, file.pdf, report.size, report.title, report.style, report.note="")
{
  scheme2 <- data.frame( min.col =c( 1,0,0 ),hinge.col=c(1,0.5,0), max.col = c(0,1,0) );
  scheme1 <- data.frame( min.col =c( 0,0.5,1 ),hinge.col=c(1,1,0), max.col = c(1,0,0) );

  get_rgbcol_hinge<-function(pvalue, p.min=0.001,p.hinge=0.05, p.max=1, log10=T, scheme=2)
  {
    if(is.na(pvalue)==T) return(c(0,0,0))
    if(pvalue > p.max) pvalue <- p.max;
    if(pvalue < p.min) pvalue <- p.min;

    if(log10)
    {
      pvalue <- log10(pvalue)
      p.min  <- log10(p.min)
      p.hinge  <- log10(p.hinge)
      p.max  <- log10(p.max)
    }


    if(scheme==1) col.scheme <- scheme1 else col.scheme <- scheme2;

    if(pvalue<=p.hinge) {
      col.scheme <- col.scheme[,1:2]
      c1 <- (col.scheme[1,1] - col.scheme[1,2]) / (p.min-p.hinge)*(pvalue-p.hinge) + col.scheme[1,2];
      c2 <- (col.scheme[2,1] - col.scheme[2,2]) / (p.min-p.hinge)*(pvalue-p.hinge) + col.scheme[2,2];
      c3 <- (col.scheme[3,1] - col.scheme[3,2]) / (p.min-p.hinge)*(pvalue-p.hinge) + col.scheme[3,2];
    }

    else{
      col.scheme <- col.scheme[,2:3]
      c1 <- (col.scheme[1,1] - col.scheme[1,2]) / (p.hinge-p.max)*(pvalue-p.max) + col.scheme[1,2];
      c2 <- (col.scheme[2,1] - col.scheme[2,2]) / (p.hinge-p.max)*(pvalue-p.max) + col.scheme[2,2];
      c3 <- (col.scheme[3,1] - col.scheme[3,2]) / (p.hinge-p.max)*(pvalue-p.max) + col.scheme[3,2];
    }

    c <- rgb(c1,c2,c3)
    return(c);
  }


  drawlegend_hinge<-function( x0, y0, width, height, title, pv.min=1e-6,pv.hinge=0.05,pv.max=1, pv.log10=T, scheme=1)
  {
    bar.len <- 50;
    if(pv.log10)
      pv <- 10^seq(log10(pv.min), log10(pv.max), length.out=bar.len)
    else
      pv <- seq(pv.min, pv.max, length.out=bar.len);

    bar.width <- width*0.6;

    col.list<-c()

    for( i in 1:length(pv) )
    {
      col <- get_rgbcol_hinge( pv[i], pv.min, pv.hinge, pv.max, log10=pv.log10, scheme=scheme);
      col.list<-c(col.list,col)
      grid.rect(  x0 + 0.3*width + i*bar.width/bar.len,
            y0,
            width = bar.width/bar.len,
            height = height,
            gp = gpar(fill=col, col=NA),
            just = c("right", "centre") );
    }


    grid.text(  title,
          x0 + 0.25*width-0.02,
          y0,
          gp = gpar(cex=0.6),
          just = c("right",  "centre"));
    grid.text( sprintf("%1g", pv.min),
          x0 + 0.28*width,
          y0 - height,
          gp = gpar(cex=0.6),
          just = c("centre", "centre"));
    grid.text( sprintf("%1g", pv.max),
          x0 + 0.9*width,
          y0 - height,
          gp = gpar(cex=0.6),
          just = c("centre", "centre"));
    grid.text( sprintf("%1g", pv.hinge),
          x0 + 0.6*width,
          y0 - height,
          gp = gpar(cex=0.6),
          just = c("centre", "centre"));
  }



  drawlegend_circle<-function( x0, y0, width, height)
  {

    pval.vec<-c(1:4,4)

    base.r<-0.001
    r.adj<-pval.vec/300
    plot.r<- base.r+ r.adj

    plot.r <- plot.r *0.7/0.04/10*8

    x0.vec<-x0+c(0,0.03,0.07,0.12,0.19)

    for(i in 1:(length(x0.vec)-1))grid.circle(x = x0.vec[i],y = y0,r= plot.r[i],gp = gpar(lwd=0.5,lty="solid"))

    grid.circle(x = x0.vec[i+1],y = y0,r= plot.r[i+1],gp = gpar(lwd=1.5,lty="solid"))

    grid.text(  "-log10 p value",
          x0-0.03,
          y0,
          gp = gpar(cex=0.6),
          just = c("right",  "centre"));



    grid.text( c("1","2","3","4",">10"),
          x0.vec,
          y0 - height,
          gp = gpar(cex=0.6),
          just = c("centre", "centre"));


  }

  get_short_value<-function(val, log=T)
  {
    if(is.na(val)==T) return("NA")

    if(log)
    {
      n.log <- ceiling(log10(val));
      if(n.log > -2 || val>=0.05 )
        str <- round(val, digits=2)
      else
        str <- sprintf("%.1g", val);

      return(str);
    }
    else
      return(round(val, digits=1));

  }

  draw_bottom <- function(){
    grid.text(report.note, x=1, y=0.95,  gp=gpar(cex=0.4), just=c("right", "top"));

    x0 <- 0.65;

    k=3
    p.min  = as.numeric(as.character(report.style$extra1[k]));
    p.hinge= as.numeric(as.character(report.style$extra2[k]));
    p.max  = as.numeric(as.character(report.style$extra3[k]));
    log10  = as.numeric(as.character(report.style$extra4[k]));
    scheme = as.numeric(as.character(report.style$extra5[k]));


    drawlegend_hinge( x0, 0.5, width=0.2, height=0.4, title=report.style$header[k], p.min, p.hinge, p.max, log10, scheme );
    x0 <- x0 - 0.4;
    drawlegend_circle(x0,0.5, width=0.2, height=0.4)
  }

  draw_item<-function(r.comp.sort, idx.start, idx.stop)
  {
    M40 <- view.lines;

    y0 <- 0;
    for(i in idx.start:idx.stop)

    {
      y0 <- 1 - (i - idx.start +1 )/M40;
      for(k in 1:NCOL(motif.df))
      {
        if(report.style$style[k]=="text")
        {
          if(as.character(report.style$hjust[k]) == "left")
            grid.text(  motif.df[i, k],
                  x = report.style$position[k]+0,
                  y = y0+0.005,
                  rot = 0,
                  gp = gpar(cex=0.5),
                  check.overlap = F,
                  just = c("left", "centre") );
          if(as.character(report.style$hjust[k]) == "centre")
            grid.text( motif.df[i, k],
                  x = report.style$position[k]+report.style$width[k]/2,
                  y = y0+0.005,
                  rot = 0,
                  gp = gpar(cex=0.5),
                  check.overlap = F,
                  just=c("centre", "centre") );
          if(as.character(report.style$hjust[k])=="right")
            grid.text( motif.df[i, k],
                  x = report.style$position[k]+report.style$width[k],
                  y = y0+0.005,
                  rot = 0,
                  gp = gpar(cex=0.5),
                  check.overlap = F,
                  just = c("right", "centre") );
        }


        if(report.style$style[k]=="circo.2D")
        {
          info<-unlist(strsplit (as.character(motif.df[i,k]),","))
          sig<-as.numeric(info[2]);    #correcting pval=0 case

          pval.min<-0.0001
          highlight.circle<-F
          if(sig <1E-10) highlight.circle<-T
          if(sig <pval.min) sig <- pval.min;

          sig <- -log10(sig)

          fe.ratio<-as.numeric(info[1])

          col.fill <- get_rgbcol_hinge( fe.ratio,
          p.min = as.numeric(as.character(report.style$extra1[k])),
          p.hinge = as.numeric(as.character(report.style$extra2[k])),
          p.max = as.numeric(as.character(report.style$extra3[k])),
          log10 = as.numeric(as.character(report.style$extra4[k])),
          scheme = as.numeric(as.character(report.style$extra5[k]))) ;

          base.r <- 0.001
          r.adj  <- sig/300

          grid.circle(x = report.style$position[k] + report.style$width[k]/2.0,
                      y = y0+0.005,
                      r= base.r+ r.adj,
                      gp = gpar(lwd=0,
                      fill = col.fill,lty="blank"));

          if (highlight.circle)
             grid.circle(x = report.style$position[k] + report.style$width[k]/2.0,
                         y = y0+0.005,
                         #r = report.style$width[k]-0.02,
                         r= base.r+ r.adj,
                         gp = gpar(lwd=1,
                         lty="solid"));
         }

        if(report.style$style[k]=="logo" && !is.null(tfbs@tf_info))
        {
           pushViewport( viewport( y = y0,
                     x = report.style$position[k],
                     width = report.style$width[k],
                     height = 1/M40,
                     just = c("left","bottom")));

          idx <-which( as.character(motif.df[i,k]) == as.character(tfbs@tf_info$Motif_ID) );

          if(length(idx)>1)
          {
             warning("Multiple matrice in the tfbs object for Motif ID")
             #, first matrix (index:", idx[1], ") is used to draw logo.\n"));
             idx <- idx[1];
          }

          seqLogo( exp(t(tfbs@pwm[[idx]])), xaxis = FALSE, yaxis = FALSE);
          popViewport();
        }
      }
    }

    grid.lines( x = c(0, 1), y=c( y0-1/M40, y0-1/M40 ), gp=gpar(col="black", lwd=1) );
    return(1 - ( y0 - 1/M40) );
  }

  view.height <- 1;
  view.lines  <- 28;
  if(!is.na(file.pdf))
  {
    if (report.size=="letter")
    {
      view.height <- 10/6;
      view.lines  <- as.integer(28*10/6);
      pdf( file.pdf, width=8, height=10,pointsize=15,useDingbats=FALSE );
    }
    else{
      pdf( file.pdf, width=8, height=60,pointsize=15,useDingbats=FALSE );
      view.height <- 10/6;
      view.lines  <- 200;
      }
  }

  for(i in 1:ceiling(NROW(motif.df)/view.lines))
  {
    grid.newpage();
    pushViewport( viewport(x=0, y=0,
            width=1, height=1,
            just=c("left","bottom")) );

    pushViewport( viewport(x=0.15, y=0.90,
            width=0.7, height=0.82,
            just=c("left","top"),
            xscale = c(0, 1),
            yscale = c(0, 1) ) )

    y0 <- draw_item( motif.df, (i-1)*view.lines+1, min( i*view.lines, NROW(motif.df) ) );
    popViewport();
    pushViewport( viewport(x=0.15, y=0.9-0.82*y0,
            width=0.7, height=0.04,
            just=c("left","top"),
            xscale = c(0, 1),
            yscale = c(0, 1) ) )
        draw_bottom();
        popViewport();
  }

  popViewport();
  if(!is.na(file.pdf)) dev.off();
}

plot_motif_table<-function(motif.df, tfs, pdf.name){

  cols.to.draw <- ncol(motif.df)-3

  hinge.min <- as.character(0.2)
  hinge.med <- as.character (1)
  hinge.max <- as.character (6)

  df.style <- data.frame(
      position = numeric(0),
      width    = numeric(0),
      header   = character(0),
      hjust    = character(0),
      style    = character(0),
      extra1   = character(0),
      extra2   = character(0),
      extra3   = character(0),
      extra4   = character(0));

  df.style <- rbind(df.style,
      data.frame( position = 0.1,
            width    = 0.05,
            header   = "motif ID",
            hjust    = "left",
            style    = "text",
            extra1   = "0",
            extra2   = "0",
            extra3   = "0",
            extra4   = "0",
            extra5   = NA));

  df.style <- rbind(df.style,
      data.frame( position = 0.27,
            width    = 0.05,
            header   = "motif name",
            hjust    = "left",
            style    = "text",
            extra1   = "0",
            extra2   = "0",
            extra3   = "0",
            extra4   = "0",
            extra5   = NA));

  for (nth.col in 1:cols.to.draw){
    df.style <- rbind(df.style,
      data.frame( position = 0.4+0.03* (nth.col-1),
            width    = 0.03,
            header   = "fold enrichment",
            hjust    = "centre",
            style    = "circo.2D",
            extra1   = hinge.min,
            extra2   = hinge.med,
            extra3   = hinge.max,
            extra4   = "1",
            extra5   = "1"));
  }

  df.style <- rbind(df.style,
      data.frame( position = 0.5 + 0.03* cols.to.draw,
            width    = 0.4,
            header   = "Motif Logo",
            hjust    = "centre",
            style    = "logo",
            extra1   = "0",
            extra2   = "0",
            extra3   = "0",
            extra4   = "0",
            extra5   = "0"))

  plot_motif_info(tfbs=tfs, motif.df =motif.df, file.pdf= pdf.name, report.size="letter",
                  report.title="", report.style= df.style, report.note="")

  return(NULL)
}
