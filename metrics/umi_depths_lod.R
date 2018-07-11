# compute theoretical LOD for each locus
# Chang Xu, 12May2016, modified on 13Jun2016

# average prediction index for each barcode. set to 3.5, which is the PI for a 8-read pair barcode with 7 read pairs being the true allele. 
pi.per.barcode <- 3.5

# get input params
args <- commandArgs(TRUE)
args.meanMtDepth = as.numeric(args[1])
args.fileIn = args[2]
args.fileOut = args[3]
args.outPrefix = args[4]

# read in the MT depths
dat <- read.delim(args.fileIn, sep='|', header=F)
colnames(dat) <- c('chrom', 'locL', 'locR', 'MTs')

# compute the cutoff and the number of variant barcodes needed to call, based on 20 FP/Mb 
cutoff.20 <- 14.0 + 0.012 * args.meanMtDepth
barcode.needed.20 <- ceiling(cutoff.20 / pi.per.barcode)
print(paste("cutoff.20:", cutoff.20, "barcode.needed.20:", barcode.needed.20))

# function to find LOD = the smallest allele fraction such that P(# of variant barcodes >= barcode.needed.20) >= 0.95
#   requires at least 5 barcodes on the locus. otherwise findp is one point and uniroot will give error
findLOD <- function(barcode.depth){
   if(!is.na(barcode.depth) & is.numeric(barcode.depth) & barcode.depth >= 5) {
      findp <- function(p) {pbinom(barcode.needed.20-1, barcode.depth, p) - .05}
      lod.out <- try(uniroot(findp, interval=c(0,1))$root, silent=T)
      if('try-error' %in% class(lod.out)) {
         lod.out <- 1.0
      }
   }
   else{
     lod.out <- 1.0
   }
   return(round(lod.out,4))
}
# function to collapse same value(lod/coverage) columns
# and write a bedgraph file
output_bedgraph <- function(df,outfile,header,val_col="foo"){
#                  @param dataframe df      : The input dataframe to iterate over
#                  @param string    val_col : The column name of the value in the bedgraph
#                  @param string    outfile : The output file path
#                  @param string    header  : The header for the bedgraph file

   file_handle <- file(outfile,"w")
   cat(header,file=file_handle)
   prev_val <- NULL
   test <- c("lod","sumt")
   for (row in 1:nrow(df)) {
      val <- df[row,val_col]
      chr <- df[row, "chrom"]
      pos <- df[row, "locR"]
      if (is.null(prev_val)) {
	 prev_val <- val
	 prev_chr <- chr
	 prev_pos <- pos
	 init_pos <- pos
	 next # skip first iteration of loop
      }
      else {
	 if (prev_chr != chr) {
            out <- paste(prev_chr,"\t",as.integer(init_pos-1),"\t",as.integer(prev_pos),"\t",round(prev_val,5),"\n",sep="")
            cat(out,file=file_handle)
	    init_pos <- pos
	 }
	 else if (prev_val != val) {
            out <- paste(prev_chr,"\t",as.integer(init_pos-1),"\t",as.integer(prev_pos),"\t",round(prev_val,5),"\n",sep="")
            cat(out,file=file_handle)
	    init_pos <- pos
	 }
	 prev_val <- val
	 prev_chr <- chr
	 prev_pos <- pos
      }
   }
   # finish out last line of the file
   out <- paste(prev_chr,"\t",as.integer(init_pos-1),"\t",as.integer(prev_pos),"\t",round(prev_val,3),"\n",sep="")
   cat(out,file=file_handle)
   close(file_handle)
}


# compute LOD for each target locus 
dat$lod <- sapply(dat$MTs, findLOD)

# drop MTs column
dat$MTs <- NULL

# output bedgraph file 
header <- sprintf("track type=bedGraph name='%s.variant-calling-lod'\n",args.outPrefix)
output_bedgraph(dat,args.fileOut,header,"lod")

# output quantiles of the LOD
lod.quantiles <- quantile(dat$lod,probs=c(0.01,0.05,0.10,0.50,0.90,0.95,0.99))
write.table(lod.quantiles, paste(args.fileOut,".quantiles.txt",sep=""), sep='|', row.names=T, col.names=F, quote=F)
