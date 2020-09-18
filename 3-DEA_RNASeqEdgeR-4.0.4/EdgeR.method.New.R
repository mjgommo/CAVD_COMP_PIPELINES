# TODO: Add comment
# 
# Author: ctorroja
###############################################################################
library(WriteXLS,quietly=T, warn.conflicts=F, verbose=F)
library(edgeR,quietly=T, warn.conflicts=F, verbose=F)
library(gplots,quietly=T, warn.conflicts=F, verbose=F)
library(RColorBrewer,quietly=Tq, warn.conflicts=F, verbose=F)
library(biomaRt,quietly=T, warn.conflicts=F, verbose=F)

sink(stderr())
cat("\n\n")
sessionInfo()
sink()
message("\n\n")

EdgeR.method=function(dir.results,counts.file,targets.file,min.ns,min.cpm,setComps,batch,paired,factors,analysis.type,species,fdr,filterby,no.mit,noreplicates,commondisp,normMethod,logratioTrim,sumTrim) { 
	
	sink(stderr())
	cat("\n\n")
	print(species)
	sink()
	message("\n\n")
	
	dir.diagnostics=paste(dir.results,"/Diagnostics",sep="")
	dir.create(dir.diagnostics,showWarnings = FALSE)
	cat(paste("RESULTS DIR: ",dir.results,"\n"))
	cat("\n")
	cat("\n")
	
	#### Get DATA ###############################################
	message(paste("Load data from",counts.file))
	data=read.delim(counts.file,header=T,sep="\t",comment.char="#",row.names=1)
	counts=round(data[,1:(dim(data)[2]-1)],0)
	rownames(counts) <- rownames(data)
	
	element.id=rownames(data)
	celement.id=as.data.frame(data[,ncol(data)])
	rownames(celement.id) <- rownames(data)
	#############################################################
	
	#### Get BioMart Annotations #######################
	if (analysis.type == "miRNAs") {
		element.type <- "miRNA"
		celement.type <- "miRNA"
	} else if (analysis.type == "isoforms" ) {
		element.type <- "transcript"
		celement.type <- "gene"
	} else {
		element.type <- "gene"
		celement.type <- "transcript"
	}
	
	if (analysis.type != "miRNAs") {
		message("Get Annotations from BioMart\n\n")
		ann.all=getBM(attributes = c("chromosome_name","start_position","end_position","strand",as.character(species$gene.name),"description",as.character(paste(element.type,"biotype",sep="_")),as.character(paste("ensembl",element.type,"id",sep="_"))), 
				filters = paste("ensembl",element.type,"id",sep="_"),values = element.id,uniqueRows=T, mart = mart )
		ann.all <- ann.all[match(unique(ann.all[,paste("ensembl",element.type,"id",sep="_")]),ann.all[,paste("ensembl",element.type,"id",sep="_")]),]
		
		rownames(ann.all) <- ann.all[,paste("ensembl",element.type,"id",sep="_")]
	}
	####################################################
	
	#### Get Targets ############################################
	targets=read.delim(targets.file,header=T,sep="\t",row.names=2,comment.char="#")

	message("Targets:")
	sink(stderr())
	print(targets)
	sink()
	message("\n\n")
	#############################################################
		
	#### Define conditions ######################################
	conds <- NULL
	samples.by.cond <- NULL
	
	random=as.factor(targets$random)
	
	if (factors == 1) {
		conds=as.factor(targets$cond)
		while (!is.na(grep("f\\d",colnames(targets),perl=T)[1])) {
			targets[grep("f\\d",colnames(targets),perl=T)[1]] <- NULL
		}
	} else if (factors > 1) {
		for (r in 1:length(rownames(targets))) {
			fcond <- ""
			for (f in grep("^f\\d+",colnames(targets),perl=T,value=T) ) { #4:(4 + factors - 1)
				if (fcond == "") {
					fcond <- targets[r,f]
				} else {
					fcond <- paste(fcond,targets[r,f],sep=".")
				}
			}
			conds=rbind(conds,fcond)
		}
		conds <- as.factor(conds)
	}

	samples.by.cond <- cbind(as.character(conds),rownames(targets))
	colnames(samples.by.cond) <- c("cond","sample")
	
	message("Conditions:")
	sink(stderr())
	print(conds)
	print(samples.by.cond)
	sink()	
	message("\n\n")  
	#############################################################
		
	#### Define Conditions Colors #####################
	cond.names=levels(conds)
	color.list <- c(brewer.pal(12, "Paired"),brewer.pal(12, "Set3"))
	cond.color <- as.data.frame(cbind("cond"=cond.names,"color"=color.list[1:length(cond.names)]))
	rownames(cond.color) <- cond.names
	###################################################
	
	#########Filter out MT genes ##################
	if (analysis.type != "miRNAs") {
		mt.genes <- ann.all[which(ann.all$chromosome_name == "MT"),]
		mt.counts <- counts[match(mt.genes[,paste("ensembl",element.type,"id",sep="_")],rownames(counts)),]
		if (no.mit) {
			cat("Remove Mitochondrial genes from analysis\n")
			element.id <- element.id[-c(match(mt.genes[,paste("ensembl",element.type,"id",sep="_")],element.id))]
			celement.id <- as.data.frame(celement.id[element.id,], row.names=element.id)
			colnames(celement.id) <- celement.type
			counts <- counts[element.id,]
			ann.all <- ann.all[element.id,]
		}
		
		cat(paste("TOTAL SAMPLES: ",length(colnames(counts)),"\n",sep=""))
		cat(paste("TOTAL ",toupper(analysis.type),": ",length(element.id),"\n",sep=""))
		
		tot.nucgenes <- colSums(counts)
		tot.mitgenes <- colSums(mt.counts)
		ppmit <- tot.mitgenes/(tot.nucgenes+tot.mitgenes)
	}
	################################################

	#### Fitler Low Expressed Elements #########################
	if (!is.null(filterby)) {
		message(paste("Filter low expressed",analysis.type,"by",filterby),"\n\n")
		groups <- unique(targets[,filterby])
		isexpr <- rowSums(cpm(counts)) >= 0
		for (g in 1:length(groups)) {
			isexpr[rowSums(cpm(counts[,rownames(targets)[which(targets[,filterby] == groups[g])]])>min.cpm) < min.ns] <- F
		}
	} else {
		message(paste("Filter low expressed",analysis.type),"\n\n")
		isexpr <- rowSums(cpm(counts)>min.cpm) >= min.ns
		# To include a filter by average expression !!!
		# isexpr <- rowSums(cpm(counts) >= min.cpm) >= min.ns & rowMeans(cpm(counts)) >= 10		
	}

	counts <- counts[isexpr,]
	element.id=element.id[isexpr]
	celement.id <- as.data.frame(celement.id[element.id,], row.names=element.id)
	colnames(celement.id) <- celement.type
	
	message(paste("SELECTED",toupper(analysis.type),":",table(isexpr)[2]),"\n\n")
	cat(paste("SELECTED ",toupper(analysis.type),": ",table(isexpr)[2],"\n",sep=""))
	cat("\n")
	cat("\n")
	
	sink(stderr())
	print(counts[1:5,])
	sink()
	message("\n\n")
	############################################################
	

	#### Create DGEList Object and Normalize ##############
	y <- DGEList(counts=counts,group=conds, genes=rownames(counts))
  print(paste("Normalization Method:",normMethod))
	y <- calcNormFactors(y,method=normMethod,logratioTrim=logratioTrim,sumTrim=sumTrim)
# 	y <- calcNormFactors(y)
	#### Report Selected Samples and Constrasts ################
	cat("SAMPLES:\n")
	df <- cbind(Sample=rownames(targets),targets,
			Lib.Size=y$samples$lib.size,Norm.Factors=y$samples$norm.factors,
			CondColor=cond.color[y$samples$group,"color"])
	lns <- sum(nchar(as.matrix(df)))
	max(apply(nchar(as.matrix(df)),1,sum))
	options(width=(max(apply(nchar(as.matrix(df)),1,sum))+dim(df)[2]*4))
	print(df,row.names=FALSE)
	write.table(df,paste(dir.results,"/sample_stats.txt",sep=""),quote=F,sep="\t",row.names=F)
	options(width=80)
	cat("\n")
	cat("\n")
	
	cat("CONTRASTS DESIGN:\n")
	colnames(setComps) <- c("cond1","cond2","Name","Factor1","Factor2")
	print(setComps,quote=FALSE,row.names=FALSE)
	cat("\n")
	cat("\n")
	###########################################################
	
	#### BATCH EFFECT #########################################
	if (batch) {
		print("Correct Batch Effects")
		bfactors <- data.frame()
		bfactors <- targets[,which(regexpr("^b",colnames(targets)) == 1)]
		batch.cpms <- log2(1e06*t(t(y$counts+0.01)/y$samples$lib.size))
		batch.logCPM=ComBat(batch.cpms,cbind(Sample.Name=rownames(targets),Batch=targets$Batch,bfactors),plot.to=dir.diagnostics,write=F,par.prior=T)
		batch.counts = t(apply(batch.logCPM,1,function (x) {round((2^x)*y$samples$lib.size/10^6)}))
		y.old <- y
		y <- DGEList(counts=batch.counts,group=conds, genes=rownames(counts))
		y <- calcNormFactors(y,method=normMethod,logratioTrim=logratioTrim,sumTrim=sumTrim)
# 		y <- calcNormFactors(y)
	}
	###########################################################
	
	#### Build Design Matrix ##################################
	if (paired == T) {
		design=model.matrix(~0+conds+random,data=y$samples)
		colnames(design) <- sub("conds","",colnames(design))
	} else {
		design=model.matrix(~0+conds,data=y$samples)
		colnames(design) <- levels(y$samples$group)	
	}
	
	sink(stderr())
	cat("DESIGN MATRIX:\n")
	print(design)
	sink()
	message("\n\n")
	
	# Disperssion #############################################
	
	###########################################################
	# error handling function
	###########################################################
	tryCatch.W.E <- function(expr){
		W <- NULL
		w.handler <- function(w){ # warning handler
			W <<- w
			invokeRestart("muffleWarning")
		}
		list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),warning = w.handler),warning = W)
	}
	
	if (noreplicates) {
		print("No biological repicates exist")
		print("Using default or user provided common dispersion")	# Recommended common dispersion values are 0.01 for genetically identical model organisms and 0.16 for human data,
		y$common.dispersion <- commondisp							# according to the EdgeR Users Guide and the Chipster Manual (but see below).
																	# Technical replicates were considered to have a value of around 0.0001 or less (but see below).
		y$prior.df <- 20											# Adds to the DGEList object the parameter $prior.df, with the default value that estimateGLMTagwiseDisp would add.
																	# The prior.df value is required for glmQLFTest(), further below.
	} else {

		er.msg=tryCatch.W.E(estimateGLMCommonDisp(y, design))
		er.msg2=unlist(strsplit(strsplit(as.character(er.msg$value),":")[[1]][3]," "))
		coefs2remove=sub("\n","",er.msg2[grep("random",er.msg2)])
		if (length(coefs2remove) > 0) {
			message(paste("Removing:",coefs2remove))
			design=design[,-match(coefs2remove,colnames(design))]
		}
		
		y <- estimateGLMCommonDisp(y, design)
		y <- estimateGLMTrendedDisp(y, design)
		y <- estimateGLMTagwiseDisp(y, design,prior.df=10)
	}
	###################################################	
	
	#### Create Norm Counts Table #####################
	norm.counts <- cpm(y, normalized.lib.sizes=T)
	if (analysis.type != "miRNAs") { namesids <- ann.all[rownames(norm.counts),as.character(species$gene.name)] } else { namesids <- rownames(norm.counts) }
	write.table(cbind(ID=rownames(norm.counts),Name=namesids,norm.counts),file = paste(dir.results,"/NormCounts.txt",sep=""),row.names=F,quote=F,sep="\t")
	###################################################
	
	#### PLOT ######################################
	message("Plot TagWise Dispersion")
	png(paste(dir.diagnostics,"/ExpDispersion.png",sep=""))
	plotBCV(y,main="Experiment Dispersion")
	dev.off()
	################################################
	
	if (analysis.type != "miRNAs") {
		#### PLOT Mit Content ##########################
		message("Plot Mitochondrial Content") 
		png(paste(dir.diagnostics,"/MTContent.png",sep=""), width = 80+80*length(rownames(targets)), height = 220+220*(length(rownames(targets))/6), units = "px", pointsize = 12)
		barplot(ppmit,ylim=c(0,1),ylab="Proportion of total counts",main="Mitochondrial Counts Content", col=as.character(cond.color[targets$cond,"color"]))
		dev.off()
		################################################
	}
	
	#### PLOT ######################################
	message("Plot counts distributions")   
	ylimit <- 0
	for (f in rownames(targets)) {
		ylimit <- ylimit + quantile(counts[,f])[4]
	}
	ylimit <- ylimit/length(rownames(targets))
	png(paste(dir.diagnostics,"/GeneExpDistrib.png",sep=""), width = 80+80*length(rownames(targets)), height = 220+220*(length(rownames(targets))/6), units = "px", pointsize = 12)
	boxplot(counts,ylim=c(0,ylimit*3), las=2, yaxt="n", ylab="counts",main="Gene Expression Distribution Across Samples", col=as.character(cond.color[targets$cond,"color"]),pars=list(cex.axis=9/max(nchar(rownames(targets)))))#
  axis(2)
	dev.off()
	################################################
	
	#### PLOT ######################################
	message("Plot counts distributions") 
	png(paste(dir.diagnostics,"/GeneExpDistribFullRange.png",sep=""), width = 80+80*length(rownames(targets)), height = 220+220*(length(rownames(targets))/6), units = "px", pointsize = 12)
	boxplot(counts, las=2, yaxt="n", ylab="counts",main="Gene Expression Distribution Across Samples", col=as.character(cond.color[targets$cond,"color"]),pars=list(cex.axis=9/max(nchar(rownames(targets)))))#
  axis(2)
	dev.off()
	################################################
  
	#### PLOT ######################################
	message("Plot Normalized Gene Expression Distributions")
	ylimit <- 0
	for (f in rownames(targets)) {
		ylimit <- ylimit + quantile(norm.counts[,f])[4]
	}
	ylimit <- ylimit*2/length(rownames(targets))
	png(paste(dir.diagnostics,"/GeneExpDistribNorm.png",sep=""), width = 80*length(rownames(targets)), height = 220+220*(length(rownames(targets))/6), units = "px", pointsize = 12)
	boxplot(norm.counts,ylim=c(0,ylimit), las=2, yaxt="n", ylab="counts",main="Gene Expression Distribution Across Samples After Normalisation",col=as.character(cond.color[targets$cond,"color"]),pars=list(cex.axis=9/max(nchar(rownames(targets)))))#
  axis(2)
	dev.off()
	################################################
	
	#### PLOT ######################################
	message("Plot Normalized Gene Expression Distributions")
	png(paste(dir.diagnostics,"/GeneExpDistribNormFullRange.png",sep=""), width = 80*length(rownames(targets)), height = 220+220*(length(rownames(targets))/6), units = "px", pointsize = 12)
	boxplot(norm.counts, las=2, yaxt="n", ylab="counts",main="Gene Expression Distribution Across Samples After Normalisation",col=as.character(cond.color[targets$cond,"color"]),pars=list(cex.axis=9/max(nchar(rownames(targets)))))#
  axis(2)
	dev.off()
	################################################
  
	#### PLOT ######################################
	message("Plot Heatmap of Sample Distances")
	dist.mat=dist(t(log2(as.matrix(norm.counts)+0.01)))
	png(paste(dir.diagnostics,"/Heatmap_LogNormCPM.png",sep=""))
	heatmap(as.matrix(dist.mat),margins = c(0.5*max(nchar(rownames(targets))),0.5*max(nchar(rownames(targets)))),symm=T,main=paste("LogNormCPM With",length(rownames(norm.counts)),"expressed",analysis.type))
	dev.off()
	################################################
	
	#### PLOT ######################################
	message("Plot PC Clustering Samples")
	pc <- prcomp(t(log2(norm.counts+0.1)))
	jpeg(filename = paste(dir.diagnostics,"/PC_Clustering_LogNormCPM.jpeg",sep=""))   
	plot(pc$x,type="n",xlab=paste("PC1 (",summary(pc)[[6]][2],")",sep=""),main="LogNormCPM PC Clustering", ylab=paste("PC2 (",summary(pc)[[6]][5],")",sep=""),col=as.character(cond.color[targets$cond,"color"]))
	text(pc$x,cex=12/max(nchar(rownames(targets))),xpd=NA,labels=as.factor(row.names(targets)),col=as.character(cond.color[targets$cond,"color"]))
	dev.off()
	################################################
	
	#### PLOT ######################################
	message("Plot PC Clustering Samples by EdgeR")
	jpeg(filename = paste(dir.diagnostics,"/PC_Clustering_MDSEdgeR.jpeg",sep=""))   
	plotMDS(y,cex=12/max(nchar(rownames(targets))),main="MDS distance",labels=as.factor(row.names(targets)),col=as.character(cond.color[targets$cond,"color"]))
	dev.off()
	################################################
	
	#### Create Norm Average Counts Matrix #########
	count.col.names <- colnames(norm.counts)
	colnames(norm.counts) <- c(count.col.names)
	
	ave.expr <- as.data.frame(y$genes)
	ave.expr=cbind(ave.expr,apply(cpm(y,normalized.lib.sizes=T),1,mean))
	colnames(ave.expr) <- c("ID","AvrExp")
	
	ave.conds <- as.data.frame(y$genes)
	names(ave.conds) <- c("ID")
	for(cond in levels(conds)) {
		ave.conds <- cbind(ave.conds,apply(as.data.frame(norm.counts[,which(conds==cond)]),1,mean))
		colnames(ave.conds)[ncol(ave.conds)] <- paste("ave_",cond,sep="")
	}
	################################################
	
	#### FIT Model #################################
	message("Fit model")
	fit<-glmFit(y,design)
	#png(paste(dir.diagnostics,"/QQplot.png",sep=""))
	#gof(fit, pcutoff=0.1, adjust="BH", plot=T, main="qq-plot of genewise goodness of fit")
	#dev.off()
	################################################
	
	#### Make Contrasts ############################
	message("Make contrasts")
	contrast.names <- paste(setComps[,3],paste(paste(setComps[,4],"*(",setComps[,1],")",sep=""),paste(setComps[,5],"*(",setComps[,2],")",sep=""),sep="-"),sep="=")
	cont.matrix<-makeContrasts(contrasts=contrast.names,levels=colnames(design))
	colnames(cont.matrix) <- setComps[,3]
	################################################
	
	#### ANOVA #####################################
	print(cont.matrix)
	alrt <- glmLRT(fit, contrast=cont.matrix)
	anova.LRT <- topTags(alrt,n=length(element.id))[which(alrt$table$PValue<=fdr),]
	################################################
	
	#### DE Tests ##################################
	message("Write Output Files from Analysis")
	file.names<-paste(dir.results,setComps[,3],sep="/")
	gid.deg <- c()
		
	cat("CONTRASTS:\n")

	#### Create Mosaic Plot for SmearPlot ##########
	np <- ceiling(sqrt(length(file.names)))
	npc <- np
	if (npc > 5) { npc <- 5 }
	npr <- ceiling(length(file.names)/npc)
	png(paste(dir.diagnostics,"/SmearPlots.png",sep=""))
	par(mfrow=c(npc,npr))
	################################################
	
	for (i in 1:length(file.names)){
		
		lrt <- glmLRT(fit, contrast=cont.matrix[,i])
		
		x<-topTags(lrt,n=length(element.id),adjust.method="BH",sort.by="p.value")$table
		rownames(x) <- x$genes
    
		if (analysis.type == "miRNAs") {
			xout <- data.frame(ID=x$genes, Name=x$genes, AvrExp=ave.expr[match(x$genes,ave.expr$ID),"AvrExp"])
			x.header <- c("ID","miRNA","AvrExp")
		} else {
			xout <- data.frame(ID=x$genes, Name=ann.all[x$genes,5], AvrExp=ave.expr[match(x$genes,ave.expr$ID),"AvrExp"])  
			x.header <- c("ID",colnames(ann.all)[5],"AvrExp")
		}
		
		for (cn in c(1,2)) {
			for (co in unlist(strsplit(as.character(setComps[i,cn]),"[+-]",perl=T))) {
				co <- gsub("\\w*\\(","",co,perl=T)
				co <- gsub("\\).*","",co,perl=T)
				xout <- data.frame(xout, ave.conds[rownames(x),paste("ave_",co,sep="")] )
				x.header <- c(x.header, paste("ave_",co,sep=""))
			}
		}
		
		fc <- ifelse(x[,"logFC"]>0,2^x[,"logFC"],-(2^abs(x[,"logFC"])))

		if (analysis.type == "miRNAs") {
			xout <- data.frame(xout,
					foldChange=fc,
					x[,c("logFC","PValue","FDR")])
			x.header <- c(x.header,"foldChange","logFC","P.Value","adj.P.Val")
			
		} else {
			xout <- data.frame(xout,
					foldChange=fc,
					x[,c("logFC","PValue","FDR")],
					ann.all[rownames(x),1:4],
					ann.all[rownames(x),c(7,6)],
					celement.id[rownames(x),])
			x.header <- c(x.header,"foldChange","logFC","P.Value","adj.P.Val","chr","start","end","strand","biotype","desc",paste("ensembl",celement.type,"id",sep="_"))
		}

		for (cn in c(1,2)) {
			for (co in unlist(strsplit(as.character(setComps[i,cn]),"[+-]",perl=T))) {
				co <- gsub("\\w*\\(","",co,perl=T)
				co <- gsub("\\).*","",co,perl=T)
				for (sa in unlist(samples.by.cond[which(samples.by.cond[,"cond"] == co),"sample"])) {
					xout <- data.frame(xout, norm.counts[rownames(x),sa] )
					xout <- data.frame(xout, counts[rownames(x),sa] )
					x.header <- c(x.header,paste("Norm",sa,sep="_"))
					x.header <- c(x.header,paste("Raw",sa,sep="_"))
				}
			}
		}
				
		colnames(xout) <- x.header
		
		message(paste("Write",file.names[i],"files"))
		write.table(xout,paste(file.names[i],".txt",sep=""),row.names=FALSE,quote=FALSE,sep="\t")
		WriteXLS("xout",ExcelFileName=paste(file.names[i],".xlsx",sep=""),SheetName=c(as.character(setComps[i,3])))
		
		x.deg <- xout[which(xout$adj.P.Val<=fdr),]
		x.deg <- x.deg[order(-abs(x.deg$logFC)),]
		contrastName <- as.character(setComps[i,3])
		setComps[i,"DEGenes"] <- dim(x.deg)[1]
		cat(paste(setComps[i,3],length(x.deg$ID),sep="\t"),"\n")
		
		
		gid.deg <- c(as.character(gid.deg),as.character(x.deg$ID))
		
#		png(paste(dir.diagnostics,"/LogFCvslogCPM_",setComps[i,3],".png",sep=""))
#		plotSmear(lrt,de.tags=x.deg$ID,main=setComps[i,3])

		plot(log2(xout$AvrExp),xout$logFC,type="n",main=paste(contrastName,"SmearPlot"),
				xlab="Log2(AvrExp)",ylab="log2FC")
		grid(col="cadetblue1")
		points(log2(xout$AvrExp),xout$logFC,pch = ".")
		points(log2(x.deg$AvrExp),x.deg$logFC,col="red",pch = ".")
		abline(h = c(-1, 1), col = "blue")
#		dev.off()
		assign(contrastName,x.deg)

	}  
	###############################################

  	dev.off()
	cat("\n")
	cat("\n")
	
	message("Write DEG Excel File")
	WriteXLS(as.character(setComps[,3]),ExcelFileName=paste(dir.results,"/DEG_AllContrasts.xlsx",sep=""))
	
	message("Write DEG_AllContrats.stats file")
	
	reps <- table(targets$cond)
	rownames(setComps) <- setComps[,3]
	
	de.stats <- data.frame(Contrasts=setComps[,"Name"],
						DEGenes=setComps[,"DEGenes"],
						Cond1=as.character(setComps[,"cond1"]),
						nCond1=reps[as.character(setComps[,"cond1"])],
						Cond2=as.character(setComps[,"cond2"]),
						nCond2=reps[as.character(setComps[,"cond2"])],
						Disp=rep(y$common.dispersion,dim(setComps)[1]))

	write.table(de.stats,paste(dir.results,"/DEG_AllContrasts.stats",sep=""),col.names=T, row.names=F, quote=F, sep="\t")

	
	#### Plot HeatMaps ############################
	if (length(gid.deg) >= 2) {
		message(paste("Plot HeatMap with Differentialy Expressed",analysis.type))
		a <- as.matrix(log2(norm.counts[unique(gid.deg),1:ncol(norm.counts)]+0.01))
		# Calculate zscore by row
		a <- t(scale(t(a), center = TRUE, scale = TRUE))
		# Hierarquical Clustering with Euclidean Distance
		colv<- as.dendrogram(hclust(as.dist(sqrt(2*dim(a)[1]*(1-cor(a)))), method="average"))
		rowv<- as.dendrogram(hclust(as.dist(sqrt(2*dim(a)[2]*(1-cor(t(a))))), method="average"))
		jpeg(filename = paste(dir.diagnostics,"/heatmapSamples_DE_LogNormCounts.jpeg",sep=""), width =1000, height = 1000,quality=100,units = "px", pointsize = 12, bg="white",res = NA)
		heatmap.2(a, Rowv=rowv, Colv=colv, col=redgreen(75), scale="none",ColSideColors = as.character(cond.color[targets$cond,"color"]),key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, margins=(c(0.7*max(nchar(rownames(targets))),6)))
		dev.off()
	}
	
	if (dim(anova.LRT$table)[1] >= 2) {
		message(paste("Plot HeatMap with LRT detected",analysis.type))
    
		a <- as.matrix(log2(norm.counts[anova.LRT$table$genes,1:ncol(norm.counts)]+0.01))
		a <- t(scale(t(a), center = TRUE, scale = TRUE))
		colv<- as.dendrogram(hclust(as.dist(sqrt(2*dim(a)[1]*(1-cor(a)))), method="average"))
		rowv<- as.dendrogram(hclust(as.dist(sqrt(2*dim(a)[2]*(1-cor(t(a))))), method="average"))
		jpeg(filename = paste(dir.diagnostics,"/heatmapSamples_LRT_LogNormCounts.jpeg",sep=""), width =1000, height = 1000,quality=100,units = "px", pointsize = 12, bg="white",res = NA)
		heatmap.2(a, Rowv=rowv, Colv=colv, col=redgreen(75), scale="none",ColSideColors = as.character(cond.color[targets$cond,"color"]),key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,margins=(c(0.7*max(nchar(rownames(targets))),6)))
		dev.off()
	}
	##########################################
}

