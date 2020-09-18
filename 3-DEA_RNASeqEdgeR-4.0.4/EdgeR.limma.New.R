# TODO: Add comment
# 
# Author: ctorroja
###############################################################################
library(WriteXLS,quietly=T, warn.conflicts=F, verbose=F)
library(edgeR,quietly=T, warn.conflicts=F, verbose=F)
library(gplots,quietly=T, warn.conflicts=F, verbose=F)
library(RColorBrewer,quietly=Tq, warn.conflicts=F, verbose=F)
library(biomaRt,quietly=T, warn.conflicts=F, verbose=F)

EdgeR.limma=function(dir.results,counts.file,targets.file,min.ns,min.cpm,setComps,batch,paired,factors,analysis.type,species,fdr,filterby,no.mit){
	
	sink(stderr())
	cat("\n\n")
	sessionInfo()
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
	
	if (paired == T) {   
		random=as.factor(targets$random)
		random.names=levels(as.factor(targets$random))
	}
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
	
	#### Fitler Out Low Expressed Elements #########################
	if (!is.null(filterby)) {
		message(paste("Filter low expressed",analysis.type,"by",filterby),"\n\n")
		groups <- unique(targets[,"random"])
		isexpr <- rowSums(cpm(counts)) >= 0
		for (g in 1:length(groups)) {
			isexpr[rowSums(cpm(counts[,rownames(targets)[which(targets$random == groups[g])]])>min.cpm) < min.ns] <- F
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
	
	#### Build Design Matrix ##################################
	design=model.matrix(~0+conds)
	colnames(design)=cond.names
	
	sink(stderr())
	cat("DESIGN MATRIX:\n")
	print(design)
	sink()
	message("\n\n")
	
	#####################################
	
	#### Normalize Data #################
	message("Normalizing Data")
	nf=calcNormFactors(counts)
	raw.lib.size <- apply(counts,2,sum)
	jpeg(filename = paste(dir.diagnostics,"/ExpDispersion.jpeg",sep=""))
	y=voom(counts,design,plot=T,lib.size=colSums(counts)*nf)
	dev.off()
	y$genes=element.id
	#####################################
	
	#### Report Selected Samples and Constrasts ################
	cat("SAMPLES:\n")
	df <- cbind(Sample=rownames(targets),
			targets,
			Lib.Size=raw.lib.size,
			Norm.Factors=nf,
			CondColor=cond.color[targets$cond,"color"])
	lns <- sum(nchar(as.matrix(df)))
	max(apply(nchar(as.matrix(df)),1,sum))
	options(width=(max(apply(nchar(as.matrix(df)),1,sum))+dim(df)[2]*4))
	print(df,row.names=FALSE)
	write.table(df,paste(dir.results,"/sample_stats.txt",sep=""),quote=F,sep="\t",row.names=F)
	options(width=80)
	cat("\n")
	cat("\n")
	
	cat("CONTRASTS DESIGN:\n")
	colnames(setComps) <- c("cond1","cond2","Name","Factor")
	print(setComps,quote=FALSE,row.names=FALSE)
	cat("\n")
	cat("\n")
	###########################################################
	
	#### BATCH EFFECT #########################################
	if (batch) {
		print("Correct Batch Effects")
		bfactors <- data.frame()
		bfactors <- targets[,which(regexpr("^b",colnames(targets)) == 1)]
		y.old <- y
		batch.logCPM <- ComBat(y$E[,rownames(targets)],cbind(Sample.Name=rownames(targets),Batch=targets$Batch,bfactors),plot.to=dir.diagnostics,write=F,par.prior=T)
		y$E <- batch.logCPM
		norm.counts.old <- 2^y.old$E[,rownames(targets)]
	}
	###########################################################
	
	#### Create Norm Counts Table #####################
	norm.counts <- 2^y$E[,rownames(targets)]
	if (analysis.type != "miRNAs") { namesids <- ann.all[rownames(norm.counts),as.character(species$gene.name)] } else { namesids <- rownames(norm.counts) }
	write.table(cbind(ID=rownames(norm.counts),Name=namesids,norm.counts),file = paste(dir.results,"/NormCounts.txt",sep=""),row.names=F,quote=F,sep="\t")
	
	#### PLOT ######################################
	message("Plot Gene Expression distributions")   
	ylimit <- 0
	for (f in rownames(targets)) {
		ylimit <- ylimit + quantile(counts[,f])[4]
	}
	ylimit <- ylimit/length(rownames(targets))
	png(paste(dir.diagnostics,"/GeneExpDistrib.png",sep=""), width = 80*length(rownames(targets)), height = 220*(length(rownames(targets))/6), units = "px", pointsize = 12)
	boxplot(counts,ylim=c(0,ylimit*3), las=2, ylab="counts",main="Gene Expression Distribution Across Samples", col=as.character(cond.color[targets$cond,"color"]), pars=list(cex.axis=9/max(nchar(rownames(targets)))))
	dev.off()
	################################################
	
	#### PLOT ######################################
	message("Plot Normalized Gene Expression Distributions")
	ylimit <- 0
	for (f in rownames(targets)) {
		ylimit <- ylimit + quantile(norm.counts[,f])[4]
	}
	ylimit <- ylimit*2/length(rownames(targets))
	png(paste(dir.diagnostics,"/GeneExpDistribNorm.png",sep=""), width = 80*length(rownames(targets)), height = 220*(length(rownames(targets))/6), units = "px", pointsize = 12)
	boxplot(norm.counts,ylim=c(0,ylimit), las=2, ylab="counts",main="Gene Expression Distribution Across Samples After Normalisation",col=as.character(cond.color[targets$cond,"color"]), pars=list(cex.axis=9/max(nchar(rownames(targets)))))
	dev.off()
	################################################
	
	#### PLOT Mit Content ##########################
	message("Plot Mitochondrial Content") 
	png(paste(dir.diagnostics,"/MTContent.png",sep=""), width = 80+80*length(rownames(targets)), height = 220+220*(length(rownames(targets))/6), units = "px", pointsize = 12)
	barplot(ppmit,ylim=c(0,1),ylab="Proportion of total counts",main="Mitochondrial Counts Content", col=as.character(cond.color[targets$cond,"color"]))
	dev.off()
	################################################
	
	#### PLOT ######################################
	message("Plot Heatmap of Sample Distances")
	png(paste(dir.diagnostics,"//Heatmap_LogNormCPM.png",sep=""))
	dist.mat=dist(t(y$E[,rownames(targets)]))
	heatmap(as.matrix(dist.mat),margins = c(0.5*max(nchar(rownames(targets))),0.5*max(nchar(rownames(targets)))),symm=T,main=paste("With",length(rownames(norm.counts)),"expressed",analysis.type))
	dev.off()
	################################################
	
	#### PLOT ######################################
	message("Plot PC Clustering Samples")
	pc <- prcomp(t(y$E[,rownames(targets)]))
	jpeg(filename = paste(dir.diagnostics,"/PC_Clustering_LogNormCPM.jpeg",sep=""))   
	plot(pc$x,type="n", main="LogNormCPM PC Clustering", xlab=paste("PC1 (",summary(pc)[[6]][2],")",sep=""), ylab=paste("PC2 (",summary(pc)[[6]][5],")",sep=""),col=as.character(cond.color[targets$cond,"color"]))
	text(pc$x,cex=12/max(nchar(rownames(targets))),xpd=NA,labels=as.factor(row.names(targets)),col=as.character(cond.color[targets$cond,"color"]))
	dev.off()
	################################################
	
	#### PLOT ######################################
	message("Plot PC Clustering Samples by EdgeR")
	jpeg(filename = paste(dir.diagnostics,"/PC_Clustering_BCVEdgeR.jpeg",sep=""))   
	plotMDS(y,cex=12/max(nchar(rownames(targets))),main="BCV Distance",labels=as.factor(row.names(targets)),col=as.character(cond.color[targets$cond,"color"]))
	dev.off()
	################################################
	
	#### Create Norm Average Counts Matrix #########
	ave.expr <- as.data.frame(y$genes,row.names=y$genes)
	ave.expr=cbind(ave.expr,apply(2^y$E,1,mean))
	colnames(ave.expr) <- c("ID","AvrExp")
	
	ave.conds <- as.data.frame(y$genes,row.names=y$genes)
	names(ave.conds) <- c("ID")
	for(cond in levels(conds)) {
		ave.conds <- cbind(ave.conds,apply(as.data.frame(2^y$E[,which(conds==cond)]),1,mean))
		colnames(ave.conds)[ncol(ave.conds)] <- paste("ave_",cond,sep="")
	}
	################################################
	
	#### FIT Model #################################
	message("Fit model")
	if (paired == T) {
		corfit<-duplicateCorrelation(y,design,block=random)
		print(paste("Correlación del Test Pareado:",corfit$consensus))
		warnings()
		fit<-lmFit(y,design,block=random,cor=corfit$consensus)
	} else {
		fit<-lmFit(y,design)
	}
	###############################################
	
	#### Make Contrasts ###########################
	message("Make contrasts")
#	contrast.names <- paste(setComps[,3],paste(paste(setComps[,4],"*(",setComps[,1],")",sep=""),paste(setComps[,4],"*(",setComps[,2],")",sep=""),sep="-"),sep="=")
#	cont.matrix<-makeContrasts(contrasts=contrast.names,levels=levels(conds))
	contrast.names <- paste(setComps[,3],paste(paste(setComps[,4],"*(",setComps[,1],")",sep=""),paste(setComps[,5],"*(",setComps[,2],")",sep=""),sep="-"),sep="=")
	cont.matrix<-makeContrasts(contrasts=contrast.names,levels=colnames(design))
	colnames(cont.matrix) <- setComps[,3]
	fit2<-contrasts.fit(fit,cont.matrix)
	fit2<-eBayes(fit2)
	###############################################
	
	#### ANOVA Test ###############################
	message("Run ANOVA Test")
	deg.anova=fit2$genes[which(fit2$F.p.value<=fdr)]
	message(paste("ANOVA Test:",length(deg.anova)),"\n\n")
	write.table(data.frame(as.matrix(deg.anova),counts[as.matrix(deg.anova),],
					ann.all[match(as.matrix(deg.anova),ann.all[,paste("ensembl",element.type,"id",sep="_")]),],celement.id[as.matrix(deg.anova),celement.type]),
			paste(dir.results,"/Anova.txt",sep=""),quote=FALSE,sep="\t" ,
			row.names=FALSE)
	cat(paste("ANOVA Test:",length(deg.anova)),"\n\n")
	###############################################
	
	#### DE Tests #################################
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
		
		x<-as.data.frame(topTable(fit2,coef=i,genelist=element.id,number=length(element.id),adjust="BH",sort.by="P"))
		
		fc <- ifelse(x[,"logFC"]>0,2^x[,"logFC"],-(2^abs(x[,"logFC"])))
		
		xout <- data.frame(ID=x$ID, AvrExp=ave.expr[x$ID,"AvrExp"])
		
		x.header <- c("ID","AvrExp")
		
		for (cn in c(1,2)) {
			for (co in unlist(strsplit(as.character(setComps[i,cn]),"[+-]",perl=T))) {
				co <- gsub("\\w*\\(","",co,perl=T)
				co <- gsub("\\).*","",co,perl=T)
				for (sa in unlist(samples.by.cond[which(samples.by.cond[,"cond"] == co),"sample"])) {
					xout <- data.frame(xout, norm.counts[x$ID,sa] )
					xout <- data.frame(xout, counts[x$ID,sa] )
					x.header <- c(x.header,paste("Norm",sa,sep="_"))
					x.header <- c(x.header,paste("Raw",sa,sep="_"))
				}
			}
		}
		
		for (cn in c(1,2)) {
			for (co in unlist(strsplit(as.character(setComps[i,cn]),"[+-]",perl=T))) {
				co <- gsub("\\w*\\(","",co,perl=T)
				co <- gsub("\\).*","",co,perl=T)
				xout <- data.frame(xout, ave.conds[x$ID,paste("ave_",co,sep="")] )
				x.header <- c(x.header, paste("ave_",co,sep=""))
			}
		}
		
		if (analysis.type == "miRNAs") {
			xout <- data.frame(xout,
					foldChange=fc,
					x[,c("logFC","P.Value","adj.P.Val")])
			x.header <- c(x.header,"foldChange","logFC","P.Value","adj.P.Val")
			
		} else {
			
			xout <- data.frame(xout,
					foldChange=fc,
					x[,c("logFC","P.Value","adj.P.Val")],
					ann.all[match(x$ID,ann.all[,paste("ensembl",element.type,"id",sep="_")]),1:(ncol(ann.all)-1)],celement.id[x$ID,celement.type])
			
			x.header <- c(x.header, "foldChange","logFC","P.Value","adj.P.Val",
					names(ann.all)[1:(ncol(ann.all)-1)],paste("ensembl",celement.type,"id",sep="_"))
		}	
		colnames(xout) <- x.header
		
		write.table(xout,file.names[i],row.names=FALSE,quote=FALSE,sep="\t")
		WriteXLS("xout",ExcelFileName=paste(file.names[i],".xlsx",sep=""),SheetName=c(as.character(setComps[i,3])))
		
		#x.deg=topTable(fit2,coef=i,genelist=element.id,number=length(element.id),adjust="BH",p.value=fdr,sort.by="logFC")
		x.deg <- xout[which(xout$adj.P.Val<=fdr),]
		x.deg <- x.deg[order(-abs(x.deg$logFC)),]
		contrastName <- as.character(setComps[i,3])
		setComps[i,"DEGenes"] <- dim(x.deg)[1]
		cat(paste(setComps[i,3],length(x.deg$ID),sep="\t"),"\n")
		
		gid.deg <- rbind(as.data.frame(gid.deg),as.data.frame(x.deg$ID))
		
		assign(contrastName,x.deg)
		
		plot(log2(as.numeric(xout$AvrExp)+0.01),xout$logFC,type="n",main=paste(contrastName,"SmearPlot"),
				xlab="Log2(AvrExp)",ylab="log2FC")
		grid(col="cadetblue1")
		points(log2(as.numeric(xout$AvrExp)+0.01),xout$logFC,pch = ".")
		points(log2(as.numeric(x.deg$AvrExp)+0.01),x.deg$logFC,col="red",pch = ".")
		abline(h = c(-1, 1), col = "blue")
		
	}  
	dev.off()
	
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
			Disp=rep("NA",dim(setComps)[1]))
	
	write.table(de.stats,paste(dir.results,"/DEG_AllContrasts.stats",sep=""),col.names=T, row.names=F, quote=F, sep="\t")
	
	if (length(gid.deg[,]) >= 2) {
		message(paste("Plot HeatMap with Differentialy Expressed",analysis.type))
		a <- as.matrix(log2(norm.counts[unique(gid.deg[,1]),]+0.01))
		a <- t(scale(t(a), center = TRUE, scale = TRUE))
		colv<- as.dendrogram(hclust(as.dist(sqrt(2*dim(a)[1]*(1-cor(a)))), method="average"))
		rowv<- as.dendrogram(hclust(as.dist(sqrt(2*dim(a)[2]*(1-cor(t(a))))), method="average"))
		jpeg(filename = paste(dir.diagnostics,"/heatmapSamples_DE_LogNormCounts.jpeg",sep=""), width =1000, height = 1000,quality=100,units = "px", pointsize = 12, bg="white",res = NA)
		heatmap.2(a, Rowv=rowv, Colv=colv, col=redgreen(75), scale="none", ColSideColors=as.character(cond.color[targets$cond,"color"]),key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, margins=(c(0.7*max(nchar(rownames(targets))),6)))
		dev.off()
	}
	
	if (length(deg.anova) >= 2) {
		message(paste("Plot HeatMap with ANOVA detected",analysis.type))
		a <- as.matrix(log2(norm.counts[as.matrix(deg.anova),]+0.01))
		a <- t(scale(t(a), center = TRUE, scale = TRUE))
		colv<- as.dendrogram(hclust(as.dist(sqrt(2*dim(a)[1]*(1-cor(a)))), method="average"))
		rowv<- as.dendrogram(hclust(as.dist(sqrt(2*dim(a)[2]*(1-cor(t(a))))), method="average"))
		jpeg(filename = paste(dir.diagnostics,"/heatmapSamples_Anova_LogNormCounts.jpeg",sep=""), width =1000, height = 1000,quality=100,units = "px", pointsize = 12, bg="white",res = NA)
		heatmap.2(a, Rowv=rowv, Colv=colv, col=redgreen(75), scale="none",ColSideColors = as.character(cond.color[targets$cond,"color"]),key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,margins=(c(0.7*max(nchar(rownames(targets))),6)))
		dev.off()
	}
}

