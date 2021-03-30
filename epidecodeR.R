if(!require("BiocManager")) {
  install.packages("BiocManager")
  #library(BiocManager)
}

if(!require("EnvStats")){
  install.packages("EnvStats")
  #library(EnvStats)
}

if(!require("ggplot2")){
  install.packages("ggplot2")
  #library(ggplot2)
}

if(!require("rtracklayer")){
  BiocManager::install("rtracklayer")
  #library(rtracklayer)
}

if(!require("GenomicRanges")){
  BiocManager::install("GenomicRanges")
  #library(GenomicRanges)
}

if (!require("rstatix")){
  BiocManager::install("rstatix")
  #library(rstatix)
}

if (!require("ggpubr")){
  BiocManager::install("ggpubr")
  #library(ggpubr)
}

`%notin%` <- Negate(`%in%`)



epidecodeR<-function (events, deg, gtf_file, id_type, boundaries, pval, param, ints) {
  
  ## Check inputs and throw errors
  if (missing(events) || is.null(events)) {
    message ("Error: events file is mandatory. Please check users manual for accepted formats!")
    ## Throw error and exit
  }
  if (missing(deg) || is.null(deg)) {
    message ("Error: DEG file is mandatory. Please check users manual for accepted formats!")
    ## Throw error and exit
  }
  if (missing(gtf_file)) {
    gtf_file<-NULL
  }
  if (missing(id_type) || is.null(id_type)) {
    id_type="gene_name"
  }
  if (missing(boundaries) || is.null(boundaries)) {
    boundaries = 0
  }
  if (missing(pval) || is.null(pval)) {
    pval=0.05
  }
  if (missing(param) || is.null(param)) {
    param = 2
    ints = c(1, 1)
  } else if (param == 1 | param == 2 | param == 3) {
    if (missing(ints) || is.null(ints) || length(ints) != 2) {
      if (param==2) {
        ints = c(1, 4)
      } else if (param == 3) {
        ints = c(2, 4)
      }
    }
  } else {
    message ("Invalid group selected! Allowed values are 
             param = 1; 
             param = 2; 
             param = 3")
    ## Throw error
  }
  
  ## Import events file 
  file<-read.table(events, header = FALSE, row.names = NULL, stringsAsFactors = F, sep = "\t", fill = TRUE)
  if (length(colnames(file))==2) { ## The file type is list of genes and event counts
    if (! is.numeric(file[1,2]) | is.null(file[1,2]) | is.na(file[1,2])) { ## Check if header is present and if present remove first row
      file<-file[-1,]
    }
    eventcounts<-as.numeric(file[,2]) ## store counts in the eventcounts variable
    names(eventcounts)<-file[,1] ## store gene names for counts in the eventcounts variable
  } else if ((length(colnames(file))>=3)==T) { 
    if (is.null(gtf_file)) { ## Input is bed and GTF not required; continue with event count per gene on column 4
      eventcounts<-unlist(lapply (unique (file[,4]), function (x) length(which (file[,4]==x))))
      names(eventcounts)<-unique(file[,4])
    } else {
      bed<-with(file, GRanges(seqnames = file$V1, IRanges(start = file$V2, end = file$V3)))
      gtf<-import(gtf_file)
      #gtf<-gtf[gtf$type=="gene",]
      #rownames(gtf)<-c(1:length(gtf[,1]))
      gtf<-gtf[,c("type", "gene_id", "gene_name")]
      gtf<-gtf[elementMetadata(gtf)[,"type"]=="gene"]
      #gtf<-makeGRangesFromDataFrame(gtf, keep.extra.columns=T)
      gtf<-gtf+boundaries
      overlap<-suppressWarnings(findOverlaps(bed, gtf, select = "arbitrary"))
      names(overlap)<-c(1:length(overlap))
      overlap<-na.omit(overlap)
      q<-as.data.frame(bed[as.numeric(names(overlap))])
      s<-as.data.frame(gtf[overlap])
      if (id_type == "merge") {
        df<-unique(data.frame(chr=q$seqnames, start=q$start, end=q$end, id=paste(s[,"gene_id"],s[,"gene_name"], sep="|"), width=q$end-q$start, strand=s$strand))
      } else {
        df<-unique(data.frame(chr=q$seqnames, start=q$start, end=q$end, id=s[,id_type], width=q$end-q$start, strand=s$strand))        
      }

      eventcounts<-unlist(lapply (unique (df$id), function (x) length(which (df$id==x))))
      names(eventcounts)<-unique(df$id)
    }
  } else {
    message ("Error in input file! Input should be tab separated gene list and event counts or BED file of events. Exiting...")
  }
  
  ## Import DEG list
  del<-read.table(deg, header = T, row.names = NULL, stringsAsFactors = F, sep = "\t", fill = T)
  #del<-del[,1:3]
  #colnames(del)<-c("ID", "Log2FC", "Pval")
  if (! is.numeric(del[,2]) && ! is.numeric(del[,3])) {
    message ("DEG list not in proper format!
             Input should be three columns separated by tab with column 1 as id, column 2 as log2FC and column 3 as P value")
  }
  del<-na.omit(del)
  del<-del[del[,3]<=pval,]

  grp1<-del[del[,1] %notin% names(eventcounts),]
  cdfdf<-data.frame(cdfPlot(distribution = "norm", param.list = list(mean=mean(grp1[,2]), sd=sd(grp1[,2])), plot.it = F), grp="0")
  ecdfdf<-data.frame(ecdfPlot(grp1[,2], plot.it = F), grp="0")
  
  grptables<-list("0"=grp1)
  
  if (param == 1) {
    grp2<-del[del[,1] %in% names(eventcounts[eventcounts>0]),]
    if (dim(grp2)[1]>1) {
      cdfdf<-rbind(cdfdf, data.frame(cdfPlot(distribution = "norm", param.list = list(mean=mean(grp2[,2]), sd=sd(grp2[,2])), plot.it = F), grp="1+"))
      ecdfdf<-rbind(ecdfdf, data.frame(ecdfPlot(grp2[,2], plot.it = F), grp="1+"))
    } else if (dim(grp2)[1]==1) {
      ecdfdf<-rbind(ecdfdf, data.frame(ecdfPlot(grp2[,2], plot.it = F), grp="1+"))
    } else {
      message("No genes map to \"1+\" group!")
    }
    grptables[["1+"]]<-grp2
    grpcounts=c("0"=length(grp1[,1]), "1+"=length(grp2[,1]))
  }
  
  
  if (param == 2) {
    grp2name=paste0("1to", ints[2], sep = "")
    grp2<-del[del[,1] %in% names(eventcounts[eventcounts>=1 & eventcounts<=ints[2]]),]
    if (dim (grp2)[1]>1) {
      cdfdf<-rbind(cdfdf, data.frame(cdfPlot(distribution = "norm", param.list = list(mean=mean(grp2[,2]), sd=sd(grp2[,2])), plot.it = F), grp=grp2name))
      ecdfdf<-rbind(ecdfdf, data.frame(ecdfPlot(grp2[,2], plot.it = F), grp=grp2name))
    } else if (dim(grp2)[1]==1) {
      ecdfdf<-rbind(ecdfdf, data.frame(ecdfPlot(grp2[,2], plot.it = F), grp=grp2name))
    } else {
      message(paste ("No genes map to \"",grp2name,"\" group", sep = ""))
    }
    
    grp3name=paste0(ints[2]+1, "+", sep = "")
    grp3<-del[del[,1] %in% names(eventcounts[eventcounts>ints[2]]),]
    if (dim (grp3)[1]>1) {
      cdfdf<-rbind(cdfdf, data.frame(cdfPlot(distribution = "norm", param.list = list(mean=mean(grp3[,2]), sd=sd(grp3[,2])), plot.it = F), grp=grp3name))
      ecdfdf<-rbind(ecdfdf, data.frame(ecdfPlot(grp3[,2], plot.it = F), grp=grp3name))
    } else if (dim(grp3)[1]==1) {
      ecdfdf<-rbind(ecdfdf, data.frame(ecdfPlot(grp3[,2], plot.it = F), grp=grp3name))
    } else {
      message(paste ("No genes map to \"",grp3name,"\" group", sep=""))
    }
    #grptables[["1"]]<-grp2
    grptables[[grp2name]]<-grp2
    #grptables[["2+"]]<-grp3
    grptables[[grp3name]]<-grp3
    #grpcounts=c("0"=length(grp1[,1]), "1"=length(grp2[,1]), "2+"=length(grp3[,1]))
    grpcounts=c(length(grp1[,1]), length(grp2[,1]), length(grp3[,1]))
    names(grpcounts)<-names(grptables)
  }
  
  
  if (param == 3) {
    grp2name="1"
    grp2<-del[del[,1] %in% names(eventcounts[eventcounts==1]),]
    if (dim (grp2)[1]>1) {
      cdfdf<-rbind(cdfdf, data.frame(cdfPlot(distribution = "norm", param.list = list(mean=mean(grp2[,2]), sd=sd(grp2[,2])), plot.it = F), grp=grp2name))
      ecdfdf<-rbind(ecdfdf, data.frame(ecdfPlot(grp2[,2], plot.it = F), grp=grp2name))
    } else if (dim (grp2)[1]==1) {
      ecdfdf<-rbind(ecdfdf, data.frame(ecdfPlot(grp2[,2], plot.it = F), grp=grp2name))
    } else {
      message("No genes map to \"1\" group")
    }
    
    grp3name=paste0("2to", ints[2], sep="")
    grp3<-del[del[,1] %in% names(eventcounts[eventcounts>1 & eventcounts<=ints[2]]),]
    if (dim (grp3)[1]>1) {
      cdfdf<-rbind(cdfdf, data.frame(cdfPlot(distribution = "norm", param.list = list(mean=mean(grp3[,2]), sd=sd(grp3[,2])), plot.it = F), grp=grp3name))
      ecdfdf<-rbind(ecdfdf, data.frame(ecdfPlot(grp3[,2], plot.it = F), grp=grp3name))
    } else if (dim (grp3)[1]==1) {
      ecdfdf<-rbind(ecdfdf, data.frame(ecdfPlot(grp3[,2], plot.it = F), grp=grp3name))
    } else {
      message (paste("No genes map to \"", grp3name, "\" group"))
    }
    
    grp4name=paste0(ints[2]+1, "+", sep="")
    grp4<-del[del[,1] %in% names(eventcounts[eventcounts>ints[2]]),]
    if (dim (grp4)[1]>1) {
      cdfdf<-rbind(cdfdf, data.frame(cdfPlot(distribution = "norm", param.list = list(mean=mean(grp4[,2]), sd=sd(grp4[,2])), plot.it = F), grp=grp4name))
      ecdfdf<-rbind(ecdfdf, data.frame(ecdfPlot(grp4[,2], plot.it = F), grp=grp4name))
    } else if (dim(grp4)[1]==1) {
      ecdfdf<-rbind(ecdfdf, data.frame(ecdfPlot(grp4[,2], plot.it = F), grp=grp4name))
    } else {
      message(paste("No genes map to \"", grp4name, "\" group"))
    }
    grptables[["1"]]<-grp2
    grptables[[grp3name]]<-grp3
    grptables[[grp4name]]<-grp4
    grpcounts=c(length(grp1[,1]), length(grp2[,1]), length(grp3[,1]), length(grp4[,1]))
    names(grpcounts)<-names(grptables)
  }
  
  test<-data.frame(aov(Order.Statistics~grp, data = ecdfdf) %>% tukey_hsd())
  
  setClass("epidecodeR", slots=list(t="data.frame", e="data.frame", eventcounts="numeric", grptables="list", grpcounts="integer",sign.test="data.frame"))
  ro<-new("epidecodeR", t=cdfdf, e=ecdfdf, eventcounts=eventcounts, grptables=grptables, grpcounts=grpcounts, sign.test=test)
  return(ro)
}


plottingfunc<-function (objdf, title, lim, xlab, ylab) {
  p<-ggplot(objdf) +
    ggtitle(title) +
    geom_vline(xintercept = 0, color="grey", alpha=0.5, linetype="dashed") +
    geom_hline(yintercept = 0.5, color="grey", alpha=0.5, linetype="dashed") +
    xlab(xlab) +
    ylab(ylab) +
    xlim(lim) +
    ylim(0,1) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

makeplot<-function (anobj, type, lim, title, xlab, ylab) {
  if (missing(anobj)) {
    message ("No epidecodeR object provided!")
  }
  if (missing(type)) {
    type="both"
  }
  if (missing(lim)) {
    lim=c(-1,1)
  }
  if (missing(title)) {
    title=""
  }
  if (missing(xlab)) {
    xlab="Log2FC"
  }
  if (missing(ylab)) {
    ylab="Cumulative Probabilities"
  }
  pt<-plottingfunc(objdf=anobj@t, title=title, lim=lim, xlab=xlab, ylab=ylab)
  pe<-plottingfunc(objdf=anobj@e, title=title, lim=lim, xlab=xlab, ylab=ylab)
  pb<-plottingfunc(objdf = NULL, title=title, lim=lim, xlab=xlab, ylab=ylab)

  for (i in 1:length(anobj@grpcounts)) {
    if (anobj@grpcounts[i]==0) {
      pt<-pt+geom_blank(data=data.frame(cdfPlot(distribution = "norm", param.list = list(mean=0, sd=1), plot.it = F), 
                    grp=names(anobj@grpcounts[i])), aes(Quantiles, Cumulative.Probabilities, group=grp, color=grp))
      pe<-pe+geom_blank(data=data.frame(cdfPlot(distribution = "norm", param.list = list(mean=0, sd=1), plot.it = F), 
                                        grp=names(anobj@grpcounts[i])), aes(Quantiles, Cumulative.Probabilities, group=grp, color=grp))
      pb<-pb+geom_blank(data=data.frame(cdfPlot(distribution = "norm", param.list = list(mean=0, sd=1), plot.it = F), 
                                        grp=names(anobj@grpcounts[i])), aes(Quantiles, Cumulative.Probabilities, group=grp, color=grp))
    } else if (anobj@grpcounts[i]==1) {
      message(paste("Theoretical plot for ", names(anobj@grpcounts[i]), " group not included since only one gene found to belong to the group!", sep = ""))
    }
  }
  pt<-pt+geom_line(aes(Quantiles, Cumulative.Probabilities, group=grp, color=grp))
  pt<-pt+scale_color_discrete(name="Groups (genes)", breaks=names(anobj@grpcounts), 
                          labels=paste(names(anobj@grpcounts)," (",anobj@grpcounts,")", sep = "" ), drop=FALSE)
  pe<-pe+geom_point(aes(Order.Statistics, Cumulative.Probabilities, group=grp, color=grp), size=1) +
         geom_line(aes(Order.Statistics, Cumulative.Probabilities, group=grp, color=grp))
  pe<-pe+scale_color_discrete(name="Groups (genes)", breaks=names(anobj@grpcounts), 
                              labels=paste(names(anobj@grpcounts)," (",anobj@grpcounts,")", sep = "" ), drop=FALSE)
  pb<-pb+geom_line(data=anobj@t, aes(Quantiles, Cumulative.Probabilities, group=grp, color=grp))
  pb<-pb+geom_point(data=anobj@e, aes(Order.Statistics, Cumulative.Probabilities, group=grp, color=grp),size=1) +
         geom_line(data=anobj@e, aes(Order.Statistics, Cumulative.Probabilities, group=grp, color=grp), linetype="dashed")
  pb<-pb+scale_color_discrete(name="Groups (genes)", breaks=names(anobj@grpcounts), 
                              labels=paste(names(anobj@grpcounts)," (",anobj@grpcounts,")", sep = "" ), drop=FALSE)
  if (type == "t") {
    p<-pt
  } else if (type == "e") {
    p<-pe
  } else if (type == "both") {
    p<-pb
  }
  return(p)
}

plot.test<-function (obj, title, ylab) {
    
    if (missing(title) || is.null(title)) {
      title = ""
    }
    if (missing(ylab) || is.null(ylab)) {
      ylab="Log2FoldChange"
    }
    bplot<-ggplot()+
    ggtitle(title) +
    xlab("Groups") +
    ylab(ylab) +
    theme_classic()
  
  over<-NULL  
    
  for (i in 1:length(obj@grpcounts)) {
    if (obj@grpcounts[i]==0) {
      bplot<-bplot+geom_blank(data=data.frame(x=rnorm(100), grp=names(obj@grpcounts[i])), aes(grp, x, group=grp, color=grp))
    } else {
      bplot<-bplot+geom_boxplot(data=(subset(obj@e, obj@e$grp %in% names(obj@grpcounts[i]))), aes(grp, Order.Statistics, group=grp, color=grp),
                                width=0.4)
    }
  }
  
  #for (i in 1:dim(obj@sign.test)[1]) {
  #  over[i]=max(obj@e$Order.Statistics)+2*i
  #}
  
  if (length (unique(obj@sign.test$p.adj.signif))==1 && unique(obj@sign.test$p.adj.signif)=="ns") {
    bplot<-bplot + 
      scale_color_discrete(name="Groups (genes)", breaks=names(obj@grpcounts), 
                           labels=paste(names(obj@grpcounts)," (",obj@grpcounts,")", sep = "" ), drop=FALSE) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
  } else {
    bplot<-bplot + 
      scale_color_discrete(name="Groups (genes)", breaks=names(obj@grpcounts), 
                           labels=paste(names(obj@grpcounts)," (",obj@grpcounts,")", sep = "" ), drop=FALSE) +
      stat_pvalue_manual(obj@sign.test, hide.ns = T, label = "p.adj", y.position = max(obj@e$Order.Statistics)*1.2, step.increase = 0.2) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
  }
  
  
  return(bplot)
}

