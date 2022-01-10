library (methylKit)
library (graphics)
args <- commandArgs(trailingOnly = T)

wkdir <- args[1]
interval <- as.numeric(args[2])
dir <- args[3]
diff <- as.numeric(args[4])
qvalue <- as.numeric(args[5])
in1 <- args[6]
in2 <- args[7]
output <- args[8]
id1 <- args[9]
id2 <- args[10]
locount <- as.numeric(args[11])

setwd(wkdir)
genome <- args[12]
refgene <- args[13]
refcpg <- args[14]

fn <- function(input_a1,input_b1,out,length){
    pdf(paste(dir,"/QC/QC_", out,"_",length,".pdf", sep=""))
    myobj <- read(list(input_a1, input_b1), sample.id=list(id1, id2), assembly=genome, treatment=c(0,1), context="CpG")
    gc();gc()   
    filtered.myobj <- filterByCoverage(myobj, lo.count = locount, lo.perc = NULL, hi.count = NULL, hi.perc = NULL)
    getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
    getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
    getCoverageStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
    getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
    remove(myobj)
    gc();gc()   

    tiles  <-  tileMethylCounts(filtered.myobj,win.size=interval,step.size=interval)

    meth2  <-  unite(tiles, destrand=FALSE)
    remove(tiles)

    getCorrelation(meth2,plot=TRUE)

    myDiff  <-  calculateDiffMeth(meth2)
    myDiff25p <- getMethylDiff(myDiff,difference=as.numeric(diff),qvalue=as.numeric(qvalue))
    remove(myDiff)

    diffMethPerChr(myDiff25p,plot=TRUE)
    write.table (myDiff25p, paste0(dir,"/", out,".txt") ,sep="\t", quote=F)
}

fn1 <- function(length){
  try(fn( in1, in2, output, length))
}

try(fn1(interval))
