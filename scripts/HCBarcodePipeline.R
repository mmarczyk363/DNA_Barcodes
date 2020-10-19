## ===================================================
## High Complexity Barcode Library Analysis Pipeline
## 6.9.2017 CH, Yale U
## ===================================================

library(Biobase)
library(Biostrings)
library(ShortRead)

#  ----------------------------------------------.------------------
# |     BC1 (18)   | CGAA |    BC2 (18)   | CGAA |    JUNK (32)     |   
#  ----------------------------------------------.------------------

## ---------------------------------------------------------------------------------------------------
preprocReads <- 
function(fqseq, spacer="CGAA", maxMis=0, minSpacerCt=2, minQual=25, nLowQBases=2) {
# Extract barcodes after quality filtering
# 6.9.2017 CH

  stopifnot(inherits(fqseq, "ShortReadQ"))
	
	# 1. Keep reads that have both spacers
	cat("- Scanning for spacer sequences ...\n")
	matches <- vmatchPattern(spacer, sread(fqseq), max.mismatch= maxMis)
	res <- elementNROWS(matches) >= minSpacerCt
	cat("Reads with both spacers:", round(100*sum(res)/length(res),digits=2),"%\n")
	
	# 2. Trim sequence after second spacer ("junk")
	cat("- Trimming reads ...\n")
	matches <- matches[res]
	ends <- as(end(matches), "matrix")[, 1:2]
	fqseq <- narrow(fqseq[res], 1, ends[, 2])
	rm(res)
	
	# 3. Keep reads with Phred scores >= minQual with nLowQBases exceptions
	cat("- Filtering for base quality ...\n")
	qmat <- as(quality(fqseq), "matrix") <= minQual
	qcount <- rowSums(qmat, na.rm=TRUE)
	
	qpass <- qcount <= nLowQBases
	fqseq <- fqseq[qpass]
	matches <- matches[qpass]
	qmat <- qmat[qpass, ]
	ends <- ends[qpass, ]
	rm(qcount)
	cat("Reads that pass QC:", round(100*sum(qpass)/length(qpass),digits=2),"%\n")
	
	starts <- as(start(matches), "matrix")
		
	# 4. Replace reads with Phred scores < minQual with N
	cat("- Replacing low quality reads with N ...\n")
	rseq <- sread(fqseq)
	rm(fqseq)
	
	rangs <- apply(qmat, 1, function(s) {
		sm <- (1:length(s))[s]
		sm <- sm[!is.na(sm)]
		IRanges(start=sm, width=rep(1, length(sm)))
	})
	rm(qmat)
	
	rangs <- do.call("IRangesList", rangs)
	rseq <- replaceAt(rseq, rangs, "N")
	rm(rangs)
	
	# 5. Split reads to extract barcodes
	cat("- Extracting barcodes ...\n")
	res <- DNAStringSetList(
		BC1 = narrow(rseq, 1, starts[, 1] - 1),
		BC2 = narrow(rseq, ends[, 1] + 1, starts[, 2] - 1))
	res
}
## ---------------------------------------------------------------------------------------------------
assembleBCs <- 
function(bcseq, bclib=bc18dict, bcWidth=c(16,20), maxMM=2) {
# Match barcodes to dictionary
# 6.9.2017 CH
# 10.4.2018 MM
  
	stopifnot(inherits(bcseq, "DNAStringSetList"))
	
  # 1. Remove too short or too long barcodes
  lenBC1 <- nchar(bcseq$BC1); lenBC2 <- nchar(bcseq$BC2)
  ind <- lenBC1 >= bcWidth[1] & lenBC1 <= bcWidth[2] & lenBC2 >= bcWidth[1] & lenBC2 <= bcWidth[2]
  cat("Reads with both barcodes from 16 to 20 bp:", round(100*sum(ind)/length(ind),digits=2),"%\n")
  bcseq <- DNAStringSetList(lapply(bcseq, "[", ind))
  
  # 2. Match to library dictionary (increase MM sequentially)
  bcseq_part <- bcseq
  n_seq <- length(bcseq$BC1)
  libmatch <- list(BC1=integer(n_seq),BC2=integer(n_seq))
  ind_list <- 1:n_seq
  for (b in 0:maxMM){
    match_tmp <- lapply(bcseq_part, function(s) 
      vwhichPDict(PDict(ct.bcodes18.rev, max.mismatch = b), s, max.mismatch=b, with.indels=F))
    ind_tmp <- do.call("&", lapply(match_tmp, function(z) unlist(lapply(z, function(s) length(s) == 1))))
    # cat("Reads with single barcode matched:", round(100*sum(ind_tmp)/length(ind_tmp),digits=2),"% when MM =",b,"\n")
    
    #store results and remove single matches to dictionary from pool of sequences remained to match
    libmatch$BC1[ind_list[ind_tmp]] <- unlist(match_tmp$BC1[ind_tmp])
    libmatch$BC2[ind_list[ind_tmp]] <- unlist(match_tmp$BC2[ind_tmp])
    bcseq_part <- DNAStringSetList(lapply(bcseq_part, "[", !ind_tmp))
    ind_list <- ind_list[!ind_tmp]
  }
  cat("Reads matched to dictionary:", round(100*sum(libmatch$BC1 > 0)/n_seq,digits=2),"% \n")
  # cat("Reads matched:", sum(libmatch$BC1 > 0),"\n")
  
  # 3. keep only reads matching library
  res <- IntegerList(BC1=unlist(libmatch$BC1[libmatch$BC1 > 0]), BC2=unlist(libmatch$BC2[libmatch$BC1 > 0]))
  res
}

## ---------------------------------------------------------------------------------------------------
getBCAbundance <-
function(bcidx) {
# Calculate abundance of each unique barcode
# 6.9.2017 CH

  stopifnot(inherits(bcidx, "IntegerList"))
	
	# 1. concatenate indices
	bcID <- do.call("paste", bcidx)
	
	abundance <- sort(table(bcID), decreasing=TRUE)
	distr <- as.data.frame(table(abundance))
	colnames(distr) <- c("nOccurences", "nReads")
	res <- list(nUniqueBC=length(abundance), bcAbundance= abundance, distribution=distr)
	res
} 

## ---------------------------------------------------------------------------------------------------
runBCPipeline <-
function(fqfln) {
# Run all pre-processing steps and get bcLib class
# 6.9.2017 CH

  if (inherits(fqfln,"ShortReadQ")){
  res <- fqfln
  } else {
    cat(paste("Reading FASTQ file:", fqfln, "\n"))
    res <- readFastq(fqfln)
  }
  nfq <- length(res)
  cat(paste("> Total reads in file:", nfq, "\n"))
  
	cat("Preprocessing...\n")
	res <- preprocReads(res)
	bcseq <- res
	npre <- as.integer(elementNROWS(res)[1])
	cat(paste("> Reads after preprocessing:", npre, "\n"))
	
	cat("Assembling...\n")
	res <- assembleBCs(res)
	bcids <- res
	nassm <- as.integer(elementNROWS(res)[1])
	cat(paste("> Reads after assembly:", nassm, "\n"))
	
	cat("Computing distribution...\n")
	res <- getBCAbundance(res)
	nubc <- res$nUniqueBC
	cat("> Unique barcodes detected:", nubc, "\n")
	
	rr <- list(n=c(reads=nfq,npre=npre,nas=nassm,nubc=nubc),bcseq=bcseq,bcid=bcids,bc=res)
	class(rr) <- "bcLib"
	rr
}

## ---------------------------------------------------------------------------------------------------
combineBCLibs2 <- 
function(lb1, lb2) {
# combine 2 BC libraries
# 10.24.2017 CH
  stopifnot(all(inherits(lb1, "bcLib"),inherits(lb2, "bcLib")))
  ns <- lb1$n + lb2$n
  
  if ("bcseq" %in% names(lb1) && "bcseq" %in% names(lb2))
    bcseq <- DNAStringSetList(BC1=append(lb1[["bcseq"]]$BC1, lb2[["bcseq"]]$BC1),
                            BC2=append(lb1[["bcseq"]]$BC2, lb2[["bcseq"]]$BC2))
  else
    bcseq <- DNAStringSetList(BC1="", BC2="")
  
  bcid <- IntegerList(BC1=append(lb1[["bcid"]]$BC1, lb2[["bcid"]]$BC1),
                      BC2=append(lb1[["bcid"]]$BC2, lb2[["bcid"]]$BC2))
  
  abd <- getBCAbundance(bcid)
  ns["nubc"] <- abd$nUniqueBC
  
  res <- list(n=ns, bcseq=bcseq, bcid=bcid, bc=abd)
  class(res) <- "bcLib"
  res
}

## ---------------------------------------------------------------------------------------------------
combineBCLibs <- 
function(lb1, lb2, ...) {
# combine many BC libraries using recursive function
# 10.24.2017 CH
  if (length(list(...)) > 0){
    res <- combineBCLibs(lb1, combineBCLibs(lb2, combineBCLibs(...)))
  } else {
    res <- lb1
    if (!missing(lb2)){
      res <- combineBCLibs2(lb1, lb2)
    }
  }
}

## ---------------------------------------------------------------------------------------------------
print.bcLib <- 
function(lb) {
# print function for bclib object
# 10.24.2017 CH
  
  cat("\nLibrary size=\t",lb$n["reads"])
  cat("\nUnique barcodes=\t",lb$n["nubc"])
  cat("\n\nBarcode sequences:\n")
  print(lb$bcseq)
  cat("\nBarcode abundance:\n")
  print(head(lb$bc$bcAbundance, n=10))
  cat("\nBarcode abundance distribution:\n")
  print(head(lb$bc$distribution, n=10))
}

## ---------------------------------------------------------------------------------------------------
plot.bcLib <- 
function(lb, scale = "orig", fRange = NULL, ...) {
# plot barcode abundance distribution
# 10.24.2017 CH
  
  data_plot <- lb$bc$distribution
  data_plot[,1] <- as.integer(levels(data_plot[,1]))[data_plot[,1]]
  if (!is.null(fRange) | length(fRange) == 2){
    data_plot <- data_plot[data_plot[,1] >= fRange[1] & data_plot[,1] <= fRange[2],]
  }
  
  
  if (scale == "percent"){
    data_plot[,1] <- 100*data_plot[,1]/max(data_plot[,1])
    xlab <- "Abundance [%]"
  } else if (scale == "orig")
    xlab = "Abundance"
  
  par(mar=c(4, 4, 4, 3.5) + 0.1)
  
  # main abundance plot
  plot(data_plot[,1], data_plot[,2], type="h", xlab=xlab, ylab="Number of Barcodes", bty="n",...)
  
  ## Cumulative Frequency Lines 
  cum_freq <- cumsum(data_plot[,2]/sum(data_plot[,2]))
  max_y <- max(data_plot[,2], na.rm = T)
  box(col = "black")
  lines(data_plot[,1], cum_freq * max_y, type = "b", cex = 0.7, pch = 19, col="blue")
  
  ## Annotate Right Axis
  axis(side = 4, at = c(0,cum_freq[1]*max_y, .5 * max_y, max_y), 
       labels = paste(c(0,round(100*cum_freq[1], digits=2),50,100) ,"%",sep=""), 
       las = 1, col.axis = "grey62", col = "blue", cex.axis = 0.8, col.axis = "blue")

  invisible()
}

## ---------------------------------------------------------------------------------------------------
venn.bc <- 
function(x, snames=NULL, n=NULL, main=NULL, fRange=NULL, ifplot=T, scale="orig") {
# plot venn diagram 
# 10.24.2017 CH
 
  stopifnot(all(sapply(x, function(s) inherits(s, "bcLib")))) 
  require(gplots)
  
  bcs <- lapply(x, function(s){ 
    if (scale == "percent")
      100 * s$bc$bcAbundance/max(s$bc$bcAbundance)
    else
      s$bc$bcAbundance
    })
  
  if (!is.null(n)){
    bcs <- lapply(bcs, head, n)
  } else {
    if (!is.null(fRange)){
      if (scale == "orig"){
        if (is.matrix(fRange)){
          bcs <- lapply(1:length(bcs), function(s) {
            fRange_tmp <- eval(parse(text=paste(fRange[s,],collapse=":")))
            bcs[[s]][bcs[[s]] %in% fRange_tmp]
          })  
        } else {
          fRange <- eval(parse(text=paste(fRange,collapse=":")))
          bcs <- lapply(bcs, function(s) s[s %in% fRange])  
        }
      } else if (scale == "percent"){
        if (is.matrix(fRange)){
          bcs <- lapply(1:length(bcs), function(s) {
            bcs[[s]][bcs[[s]] >= fRange[s,1] & bcs[[s]] <= fRange[s,2]]
          })
        } else {
          bcs <- lapply(bcs, function(s) {
            s[s >= fRange[1] & s <= fRange[2]] 
          })  
        }
        
      }
    }
  }
  
  # get names of barcodes to compare
  bcs_all <- bcs
  bcs <- lapply(bcs, names)
  if (!is.null(snames)) names(bcs) <- snames
  
  # get set interection and union to calculate Jaccard index
  inter <- intersect2(bcs)
  uni <- unique(unlist(bcs))
  Jaccard <- length(inter)/length(uni)
  if (is.null(main)) main <- paste("Jaccard index:",round(Jaccard,digits=3))
  
  # calculate correlation betweeen abundances of common barcodes
  if (length(bcs_all) == 2){
    r <- cor(as.numeric(bcs_all[[1]][inter]),as.numeric(bcs_all[[2]][inter]))
    tau <- cor(as.numeric(bcs_all[[1]][inter]),as.numeric(bcs_all[[2]][inter]),method="kendall")
    
    if (ifplot) {
      par(mfrow=c(1,2))
      venn(bcs)
      title(main=paste0(main,", J=",round(Jaccard,digits=3)))
      rep1 <- as.numeric(bcs_all[[1]][inter])
      rep2 <- as.numeric(bcs_all[[2]][inter])
      plot(rep1,rep2,main=paste0(main,", R=",round(r,digits=3),", Tau=",round(tau,digits=3)))
      reg <- lm(rep2~rep1)
      abline(reg)
      
      
    }
  } else if (ifplot){
    venn(bcs)
    title(main=paste0(main,", J=",round(Jaccard,digits=3)))
    r <- NA
    tau <- NA
  }
  
  res <- list()
  res$Jaccard <- Jaccard
  res$r <- r
  res$tau <- tau
  res$inter <- length(inter)
  res$uni <- length(uni)
  return(res)
}

## ---------------------------------------------------------------------------------------------------
intersect2 <- 
function(x){
# Calcualte intersect of 2 or more sets from the list
# 8.29.2018 MM
  
  nsets <- length(x)

  if(nsets <= 1) {
    if(nsets == 1 && is.list(x[[1]])) {
      do.call("intersect2", x[[1]])
    } else {
      stop("Cannot evaluate intersection of fewer than 2 sets.")
    }
  } else if(nsets == 2) {
    intersect(x[[1]], x[[2]])
  } else {
    intersect(x[[1]], intersect2(x[-1]))
  }
}

## ---------------------------------------------------------------------------------------------------
filterBCLib <- 
function(lb, n=NULL, fRange=NULL, fnames = NULL) {
# filter BC libraries
# 8.29.2018 MM

  stopifnot(inherits(lb, "bcLib"))
  
  bcs <- lb$bc$bcAbundance
  nbcs <- length(bcs)
  ind <- rep(TRUE,nbcs) #index of barcodes to remain
  if (!is.null(n)){
    # take only first n most abundant barcodes
    ind[(n+1):nbcs] <- FALSE
  }
  if (!is.null(fRange)){
    # filter barcodes by frequency
    stopifnot(length(fRange) == 2)
    fRange <- eval(parse(text=paste(fRange,collapse=":")))
    
    ind <- ind & (bcs %in% fRange)
  }
  if (!is.null(fnames)){
    # filter barcodes by name
    ind <- ind & !(names(bcs) %in% fnames)
  }
  
  # save no. of unique barcodes after filtering
  lb$n[["nubc"]] <- sum(ind)
  
  # get names of remaining barcodes and filter bcid
  bc_keep <- names(bcs)[ind]
  bcID <- do.call("paste", lb$bcid)
  ind <- bcID %in% bc_keep
  lb$bcid <- IntegerList(BC1=lb$bcid$BC1[ind], BC2=lb$bcid$BC2[ind])
  
  # get new abundances
  lb$bc <- getBCAbundance(lb$bcid)
  lb$n[["reads_filt"]] <- sum(lb$bc$bcAbundance)
  
  # return filtered bcLib
  res <- lb
  class(res) <- "bcLib"
  return(res)
}

## ---------------------------------------------------------------------------------------------------
plot_distr_comm <- 
function(x, snames=NULL, main=NULL, fRange=NULL, comm = "all") {
# Plot distribution of abundance for common and exclusive barcodes
# 10.1.2018 MM
  require(ggplot2)
  
  stopifnot(all(sapply(x, function(s) inherits(s, "bcLib"))))
  
  bcs <- lapply(x, function(s) s$bc$bcAbundance )
  
  if (!is.null(fRange)) {
    if (is.matrix(fRange)) {
      bcs <- lapply(1:length(bcs), function(s) {
        fRange_tmp <- eval(parse(text = paste(fRange[s, ], collapse = ":")))
        bcs[[s]][bcs[[s]] %in% fRange_tmp]
      })
    } else {
      fRange <- eval(parse(text = paste(fRange, collapse = ":")))
      bcs <- lapply(bcs, function(s)
        s[s %in% fRange])
    }
  }
  
  # get names of barcodes to compare
  bcs_all <- bcs
  bcs <- lapply(bcs, names)
  if (!is.null(snames)) names(bcs) <- snames
  
  # get set interection and union to group barcodes
  if (comm == "all"){
    inter <- rep(list(intersect2(bcs)),length(bcs))
    diff_all <- sapply(bcs,function(s) setdiff(s,inter))
  }
  else if (comm == "single"){
    diff_all <- lapply(1:length(bcs),function(s) setdiff(bcs[[s]],unique(unlist(bcs[names(bcs)[-s]]))))
    names(diff_all) <- names(bcs)
    inter <- lapply(1:length(bcs), function(s) setdiff(bcs[[s]],diff_all[[s]]))
    names(inter) <- names(bcs)
  }
  
  # format data to plot
  data_plot <- lapply(names(bcs_all),function(s){
    group <- character(length(bcs_all[[s]]))
    names(group) <- rownames(bcs_all[[s]])
    group[inter[[s]]] <- "Common"
    group[diff_all[[s]]] <- "Exclusive"
    data.frame(Abundance = as.numeric(bcs_all[[s]]), Group = group, stringsAsFactors = F)
  })
  names(data_plot) <- names(bcs_all)
  data_plot <- reshape::melt(data_plot)
  data_plot$variable <- NULL
  colnames(data_plot)[2:3] <- c("Abundance","Sample")
  
  p <- ggplot(data_plot, aes(x=Group, y=Abundance, color = Group)) + geom_violin() + 
    geom_boxplot(width=0.05) + facet_grid(. ~ Sample) + scale_y_continuous(trans='log10') + 
    scale_color_brewer(type = "qual", palette = "Set1") + theme_bw()
  return(p)
  
  # par(mfrow=c(1,length(bcs_all)))
  # lapply(names(bcs_all),function(s){
  #     group <- character(length(bcs_all[[s]]))
  #     names(group) <- rownames(bcs_all[[s]])
  #     group[inter] <- "Common"
  #     group[diff_all[[s]]] <- "Exclusive"
  #     boxplot(log10(as.numeric(bcs_all[[s]])) ~ group,ylab="Log10 abundance",main=s,xlab="Barcodes")
  #     return(1)
  # })
  # invisible()
}