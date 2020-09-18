# TODO: Add comment
# 
# Author: fscabo,jnunez,ctorroja
###############################################################################
library(optparse,quietly=T, warn.conflicts=F, verbose=F)

##############################
## DEFINE OPTION PARAMETERS ##
##############################

opls <- list()

opls[[1]] <- make_option(c('--species','-s'), action="store", type="character", dest="species", default=NULL, help="Species used as Transcriptome Reference (mouse, human, zebrafish, sheep)")
opls[[2]] <- make_option(c('--biomart','-b'), action="store", type="character", dest="bioMArchive", default="www", help="Ensembl BioMart database to look for annotations (i.e: dec2011.archive) [default: %default]")
opls[[3]] <- make_option(c('--analysis','-a'), action="store", type="character", dest="analysis.type", default=NULL, help="type of DE analysis (genes, isoforms, miRNAs)")
opls[[4]] <- make_option(c('--pipeline','-e'), action="store", type="character", dest="pipeline.type", default="edgeR.method", help="type of Pipeline for DE analysis (edgeR.limma, edgeR.method, DESeq, edgeR.NoRep) [default: %default]")
opls[[5]] <- make_option(c('--target_file','-t'), action="store", type="character", dest="target.file", default="Targets.file.txt", help="File containing sample information (with sample,cond,[random,f1,f2,...]) [default: %default")
opls[[6]] <- make_option(c('--design_file','-d'), action="store", type="character", dest="setComp.file", default="Design.file.txt", help="File containing design information (with cond1,cond2,outputPrefix,MultiplicationFactor) [default: %default")
opls[[7]] <- make_option(c('--paired','-p'), action="store_true", type="logical", dest="paired", default=FALSE, help="Whether the analysis should be paired [default: %default]")
opls[[8]] <- make_option(c('--batch'), action="store_true", type="logical", dest="batch", default=FALSE, help="Whether the counts should be corrected for batch effect [default: %default]")
opls[[9]] <- make_option(c('--factors','-f'), action="store", type="integer", dest="factors", default=1, help="Number of factors to be used in comparisons [default: %default]")
opls[[10]] <- make_option(c('--fdr'), action="store", type="numeric", dest="fdr", default=0.05, help="False Discovery Rate filter [default: %default]")
opls[[11]] <- make_option(c('--min_cpms','-c'), action="store", type="integer", dest="min.cpm", default=1, help="Min number of COunts Per Million reads to consider the element for the analysis [default: %default]")
opls[[12]] <- make_option(c('--min_samples','-l'), action="store", type="integer", dest="min.ns", default=NULL, help="Minimum number of samples that express the element to be analyzed (required)")
opls[[13]] <- make_option(c('--filterby'), action="store", type="character", dest="filterby", default=NULL, help="Filter the elements grouping by this column name in target file [default: %default]")
opls[[14]] <- make_option(c('--workdir','-w'), action="store", type="character", default=getwd(), dest="path", help="Working directory to output results [default: %default]")
opls[[15]] <- make_option(c('--input_dir','-i'), action="store", type="character", default="RSEM.", dest="dir.data", help="Directory containing the input files with expression values")
opls[[16]] <- make_option(c('--matrix_file','-m'), action="store", type="character", dest="matrix.file", default=NULL, help="File name containing the matrix of expression values for all samples")
opls[[17]] <- make_option(c('--suffix','-x'), action="store", type="character", dest="suffix", default=NULL, help="Suffix of files with counts [default: %default]")
opls[[18]] <- make_option(c('--counts_col'), action="store", type="integer", dest="counts.column", default=5, help="Column containing the counts [default: %default]")
opls[[19]] <- make_option(c('--matchID_col'), action="store", type="integer", dest="matchID.column", default=2, help="Column containing the matching IDs (transcriptIDs or geneIDs) to the ones being analyzed  [default: %default]")
opls[[20]] <- make_option(c('--id_col'), action="store", type="integer", dest="id.column", default=1, help="Column containing the IDs (transcriptIDs, geneIDs or miRNAIDs) [default: %default]")
opls[[21]] <- make_option(c('--header'), action="store_true", type="logical", dest="header", default=FALSE, help="Whether Count files contains header  [default: %default]")
opls[[22]] <- make_option(c('--sep'), action="store", type="character", dest="sep", default="\t", help="Separation between columns [default: %default]")
opls[[23]] <- make_option(c('--outprefix','-o'), action="store", type="character", dest="analysis.prefix", default="DiffAll", help="Prefix for the output files [default: %default]")
opls[[24]] <- make_option(c('--no_mit'), action="store_true", type="logical", dest="no.mit", default=FALSE, help="Remove counts from Mitochondrial genes  [default: %default]")
opls[[25]] <- make_option(c('--noreplicates'), action="store_true", type="logical", dest="noreplicates", default=FALSE, help="Whether biological replicates are absent [default: %default]")
opls[[26]] <- make_option(c('--commondisp'), action="store", type="numeric", dest="commondisp", default=0.01, help="Common dispersion [default: %default;  guidelines: human data, 0.16; genetically identical model organisms, 0.01; technical replicates, < 0.0001]")
opls[[27]] <- make_option(c('--norm'), action="store", type="character", dest="normMethod", default="TMM", help="Norm method to apply (TMM,upperquartile,RLE,none) [default: %default]")
opls[[28]] <- make_option(c('--logratioTrim'), action="store", type="numeric", dest="logratioTrim", default=0.3, help="amount of trim to use on log-ratios (\"M\" values) for method=\"TMM\" [default: %default]")
opls[[29]] <- make_option(c('--sumTrim'), action="store", type="numeric", dest="sumTrim", default=0.05, help="amount of trim to use on the combined absolute levels (\"A\" values) for method=\"TMM\" [default: %default]")

opts <- OptionParser(usage = "usage: %prog [options]",
		option_list = opls, add_help_option = TRUE,
		prog = "DiffExpRNASeq.R", description = "", epilogue = "")

args <- parse_args(opts,
		args = commandArgs(trailingOnly = TRUE),
		print_help_and_exit = TRUE,
		positional_arguments = FALSE)

#######################
## ANALYSIS SETTINGS ##
#######################
species <- args$species #### default mouse; USER INPUT
targets.file <- args$target.file #### default Targets.file.txt; USER INPUT
setComp.file <- args$setComp.file
analysis.type <- args$analysis.type #### default genes; USER INPUT
analysis.pipeline <- args$pipeline.type #### default RSEM; USER INPUT
paired <- args$paired #### default F; USER INPUT
batch <- args$batch #### default F; USER INPUT
factors <- args$factors #### default 1; USER INPUT
min.ns <- args$min.ns #### USER INPUT
min.cpm <- args$min.cpm #### default 1; USER INPUT
fdr <- args$fdr #### default 0.05; USER INPUT
filterby <- args$filterby
no.mit <- args$no.mit
normMethod <- args$normMethod
logratioTrim <- args$logratioTrim
sumTrim <- args$sumTrim
####################
## INPUT SETTINGS ##
####################
path <- args$path #### USER INPUT
dir.data <- args$dir.data #### USER INPUT
bioMArchive <- args$bioMArchive #### USER INPUT
matrix.file <- args$matrix.file ####paste("RSEM.",analysis.type,".results.matrix",sep="") #### USER INPUT
suffix <- args$suffix #### USER INPUT
noreplicates <- args$noreplicates #### default F; USER INPUT
commondisp <- args$commondisp #### default NULL; USER INPUT

#####################
## OUTPUT SETTINGS ##
#####################
analysis.prefix <- args$analysis.prefix ##### USER INPUT

###################
## DEbug Options ##
###################
#args$path <- "/home/ctorroja/Documents/Analysis/JLPompa/MetaJag1KOvsMib1KO/DE-Analysis"
#args$species <- "mouse" #### default mouse; USER INPUT
#args$setComp.file <- "Design.file.txt"
#args$target.file <- "Targets.file.txt" #### default Targets.file.txt; USER INPUT
#args$analysis.type <- "genes" #### default genes; USER INPUT
#args$pipeline.type <- "edgeR.method"
#args$paired <- TRUE #### default F; USER INPUT
#args$factors <- 1 #### default 1; USER INPUT
#args$min.ns <- 3 #### USER INPUT
#args$min.cpm <- 2 #### default 1; USER INPUT
#args$fdr <- 0.05
#args$suffix <- ".genes.results"
#args$dir.data <- "RSEM.genes"
#args$bioMArchive <- "jan2013.archive"
#args$analysis.prefix <- "Diff.Test"
#args$matrix.file <- NULL
#args$id.column <- 1
#args$matchID.column <- 2
#args$counts.column <- 5
#args$header <- T
#args$batch <- F
#args$filterby <- random #### default NULL; USER INPUT
#args$no.mit <- F
#args$noreplicates <- F #### default F; USER INPUT
#args$commondisp <- NULL #### default NULL; USER INPUT
#args$normMethod <- "TMM"
#args$logratioTrim <- 0.3
#args$sumTrim <- 0.05
########################

expath <- strsplit(commandArgs(trailingOnly=F)[4],"=",fixed=T)[[1]][2]

library(biomaRt,quietly=T, warn.conflicts=F, verbose=F)
if (analysis.pipeline == "edgeR.limma") {
  source(paste(dirname(expath),"EdgeR.limma.New.R",sep="/"))
  source(paste(dirname(expath),"ComBat.original.New.R",sep="/"))
}  else if (analysis.pipeline == "edgeR.method") {
  source(paste(dirname(expath),"EdgeR.method.New.R",sep="/"))
  source(paste(dirname(expath),"ComBat.original.New.R",sep="/"))
}


sink(stderr())
cat("\n\nSESSION INFO:\n")
sessionInfo()
sink()
message("\n\n")
message(paste("PATH:",expath),"\n\n")

sink(stderr())
cat(as.character("RUNNING ARGUMENTS:\n"))
print(t(as.data.frame(args,row.names=c("Values"))))
sink()
message("\n\n")

cat("\n\nRUNNING ARGUMENTS:\n")
print(t(as.data.frame(args,row.names=c("Values"))))
cat("\n")
cat("\n")


setwd(path)

counts.file <- paste("Counts",analysis.type,analysis.prefix,"txt",sep=".")
dir.results <- paste(path,paste(analysis.prefix,analysis.type,analysis.pipeline,normMethod,paste("n",min.ns,sep=""),paste("cpm",min.cpm,sep=""),sep="."),sep="/")
if (!is.null(filterby)) {
  dir.results <- paste(dir.results,paste("filterBy",filterby,sep=""),sep=".")
}
if (paired) { 
	dir.results <- paste(dir.results,"paired",sep=".")
}
if (batch) {
  dir.results <- paste(dir.results,"batch",sep=".")
}

# Set Contrasts: Cond1,Cond2,Name_of_Contrast,MultiplicationFactor. Cond1/Cond2, Condition/Control, ConditionVsControl
message(paste("Load Design File:",setComp.file))
setComps <- read.delim(setComp.file,header=F,sep="\t",strip.white = T,comment.char="#")

sink(stderr())
cat("DESIGN:\n")
print(setComps)
sink()
message("\n\n")

species.list <- as.data.frame(rbind(
				c("human","Homo sapiens","Homo_sapiens","hsapiens","9606","hgnc_symbol"),
				c("mouse","Mus musculus","Mus_musculus","mmusculus","10090","mgi_symbol"),
				c("zebrafish","Danio rerio","Danio_rerio","drerio","7955","zfin_symbol"),
				c("sheep","Ovis aries","Ovis_aries","oaries","9940","external_gene_id"),
				c("pig","Sus scrofa","Sus_scrofa","sscrofa","9823","external_gene_id")
		))

names(species.list) <- c("common","preferred","preferred2","reduced","taxid","gene.name")
cat("ORGANISM:\n")
print(as.data.frame(species.list[match(species,species.list$common),]),row.names=FALSE)
cat("\n")
cat("\n")

########################
## CONNECT TO BIOMART ##
########################
message(paste("Connect To BioMart at",paste(bioMArchive,"ensembl.org",sep=".")))
#Connect to Version 60 (nov2010) of ensembl archive
sink(stderr())
mart <- useMart(host=paste(bioMArchive,"ensembl.org",sep="."), path="/biomart/martservice", biomart = "ENSEMBL_MART_ENSEMBL", port=80,  verbose = TRUE)
datasets <- listDatasets(mart)

mart<-useDataset(paste(species.list[which(species.list$common == species),"reduced"],"gene","ensembl",sep="_"),mart)
sink()
message("\n\n")
filters = listFilters(mart)
attributes = listAttributes(mart)

################
## GET COUNTS ##
################

if (is.null(matrix.file)==FALSE) {
	
	cat(paste("Load Expression Matrix from file",matrix.file,sep=" "))
	message(paste("Load Expression Matrix from file",matrix.file))
	counts = read.delim(matrix.file,header=T,sep="\t",row.names=1,comment.char="#")
	# get all IDs
	message("get all IDs")
	element.id=rownames(counts)
	matchIDs = counts[element.id,length(names(counts))]
	message(paste(length(element.id),analysis.type))
	message("Counts in first record:")
	sink(stderr())
	print(counts[1,])
	sink()
	message("\n\n")
	
} else {
	
	cat("Build Expression Matrix from sample files:\n")
	message("Build Expression Matrix from sample files:")
	file.names=dir(dir.data)
	print(file.names,row.names=FALSE,quote=FALSE)
	sink(stderr())
	print(file.names,row.names=FALSE,quote=FALSE)
	sink()
	# get all IDs
	message(paste("get all IDs from column",args$id.column))
	ns=length(file.names)
	element.id=vector()
	for (i in 1:ns){
		data=read.delim(paste(dir.data,"/",file.names[i],sep=""),header=args$header,sep="\t",comment.char="#")
		element.id=union(element.id,data[,args$id.column])
	}
	
	matchIDs = data[match(element.id,data[,args$id.column]),args$matchID.column]
	
	message(paste(length(element.id),analysis.type))
	
	#Build Expression Matrix Files
	counts=matrix(0,length(element.id),ns)
	for (i in 1:ns){
		data=read.delim(paste(dir.data,"/",file.names[i],sep=""),header=args$header,sep="\t",comment.char="#")
		
		dummy=match(data[,args$id.column],element.id)
		counts[,i]=data[dummy,args$counts.column]
	}
	counts=as.data.frame(counts)
	names(counts)=sub(suffix,"",file.names)
	message("Counts in first record:")
	sink(stderr())
	print(counts[1,])
	sink()
	message("\n\n")
}
cat("\n")
cat("\n")

######################
# Load Targets file ##
######################
if (file.exists(targets.file)==FALSE) {
	
	message("Build Targets file from counts files")
	cnames <- names(counts)
	library.names <- sub("[.][[:alnum:][:punct:]]+$","",cnames[1:length(cnames)-1])
	sample.names <- sub("[[:alnum:]_]+[.]","",cnames[1:length(cnames)-1])
	cond.names <- sub("[[:alnum:]_]+[.]","",sample.names)
	ref.names <- sub("[[:alnum:]_]+[.]","",cond.names)
	sample.names <- sub("[.][[:alnum:][:punct:]]+$","",sample.names)
	cond.names <- sub("[.][[:alnum:][:punct:]]+$","",cond.names)
	ref.names <- sub("[.][[:alnum:][:punct:]]+$","",ref.names)
	
	targets <- as.data.frame(cbind("sample.file"=cnames,"sample"=sample.names,"cond"=cond.names))
	targets <- targets[sort(as.character(paste(targets$cond,targets$sample,sep="_")),index.return=T)[[2]],]
	write.table(targets,targets.file,row.names=F, sep="\t", col.names=T, quote=F)
	rownames(targets) = targets$sample.file
	targets$sample.file = NULL
	
} else {
	message(paste("Load Targets file:",targets.file))
	targets=read.delim(targets.file,header=T,sep="\t",row.names=1,comment.char="#")
}
sink(stderr())
cat("TARGETS:\n")
print(targets,quote=FALSE)
sink()
message("\n\n")

if (any(is.element(rownames(targets),colnames(counts))==FALSE)) {
	sink(stderr())
	print("ERROR -> Different filenames in counts directory from filenames in targets file.")
	print(paste("Sample file.names from files:",colnames(counts)))
	print(paste("Sample file.names in target file:",rownames(targets)))
	sink()
	quit(save="no")
}


# Reduce Samples to the ones present in the Targets.file
counts=counts[,rownames(targets)]
names(counts)=targets$sample

if (analysis.type == "isoforms") {
	name.id <- "IsoformID"
	name.to.id <- "GeneID"
} else if (analysis.type == "genes") {
	name.id <- "GeneID"
	name.to.id <- "TranscriptID"
} else if (analysis.type == "miRNAs") {
	name.id <- "miRNAID"
	name.to.id <- "pmiRNAID"
}
counts=cbind(element.id,counts,matchIDs)
names(counts)=c(name.id,as.character(targets[,"sample"]),name.to.id)

# Write Expression Matrix Files
message(paste("Write",analysis.type,"results matrix",counts.file),"\n\n")

write.table(counts, counts.file, row.names=F, sep="\t", col.names=T)

##############################################################
##  Run Analysis Pipeline
##############################################################

dir.create(dir.results,showWarnings = FALSE)

 spa <- as.data.frame(species.list[match(species,species.list$common),])
 
 if (analysis.pipeline == "edgeR.limma") {
	message(paste("Analysing",analysis.type,"differential expression with",analysis.pipeline,"pipeline"))
	EdgeR.limma(dir.results,counts.file,targets.file,min.ns,min.cpm,setComps,batch,paired,factors,analysis.type,spa,fdr,filterby,no.mit)
} else if (analysis.pipeline == "edgeR.method") {
	message(paste("Analysing",analysis.type,"differential expression with",analysis.pipeline,"pipeline"))
	EdgeR.method(dir.results,counts.file,targets.file,min.ns,min.cpm,setComps,batch,paired,factors,analysis.type,spa,fdr,filterby,no.mit,noreplicates,commondisp,normMethod,logratioTrim, sumTrim)
} else if (analysis.pipeline == "edgeR.NoRep") {
	cat("NOT IMPLEMENTED YET\n")
#	source(paste(dirname(expath),"EdgeR.NoRep.New.R",sep="/"))
#	source(paste(dirname(expath),"ComBat.original.New.R",sep="/"))
#	message(paste("Analysing",analysis.type,"differential expression with",analysis.pipeline,"pipeline"))
#	EdgeR.method(dir.results,counts.file,targets.file,min.ns,min.cpm,setComps,batch,paired,factors,analysis.type,spa,fdr,filterby,no.mit)
} else if (analysis.pipeline == "edgeR.miRNA") {
	cat("NOT IMPLEMENTED YET\n")
# 	source(paste(dirname(expath),"EdgeR.miRNA.New.R",sep="/"))
# 	source(paste(dirname(expath),"ComBat.original.New.R",sep="/"))
# 	message(paste("Analysing",analysis.type,"differential expression with",analysis.pipeline,"pipeline"))
# 	EdgeR.miRNA(dir.results,counts.file,targets.file,min.ns,min.cpm,setComps,batch,paired,factors,analysis.type,spa,fdr,filterby,no.mit)
} else if (analysis.pipeline == "DESeq") {
	cat("NOT IMPLEMENTED YET\n")
#	source(paste(dirname(expath),"DESeq.methd.R",sep="/"))
#	source(paste(dirname(expath),"ComBat.original.New.R",sep="/"))
#	message(paste("Analysing",analysis.type,"differential expression with",analysis.pipeline,"pipeline"))
#	DESeq.method(dir.results,counts.file,targets.file,min.ns,min.cpm,setComps,batch,paired,factors,analysis.type,spa,fdr,filterby,no.mit)
} else if (analysis.pipeline == "EBSeq") {
	cat("NOT IMPLEMENTED YET\n")
#	source(paste(dirname(expath),"EBSeq.methd.R",sep="/"))
#	source(paste(dirname(expath),"ComBat.original.New.R",sep="/"))
#	message(paste("Analysing",analysis.type,"differential expression with",analysis.pipeline,"pipeline"))
#	EBSeq.method(dir.results,counts.file,targets.file,min.ns,min.cpm,setComps,batch,paired,factors,analysis.type,spa,fdr,filterby,no.mit)
} else {
	message(paste("Incorrect Method:",analysis.pipeline))
}
sink(stderr())
warnings()
sink()

cat("DONE\n")
