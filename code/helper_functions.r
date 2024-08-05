prepBias <- function(bindev.results, n.top=10000, sd.interval=5, ...) {
        #description
        ### bindev.results = dataframe with gene_id, gene_name, rank_default, rank_[x] where x in the batch effect 
        ### n.top = maximum rank (default) of features to include when calculating SD
        ### sd.interval = bin width for nSD(rank difference). value n determines upper open limit: [0,n), [n,2n), [2n,3n) ... 
        
        #error conditions 
        if(!is.numeric(sd.interval)) {stop("`sd.interval` determinning bin width must be a single numeric value.")}
        tmp = grep("rank_", colnames(bindev.results), value=TRUE)
        if(length(tmp)!=2) {stop("Dataframe of binomial deviance results must have two rank columns.\nColumns must be named `rank_default` and `rank_[batch]`.\nRe-run binomial deviance or reformat results table.\n")}
        if(length(intersect(tmp,"rank_default"))==0) {stop("One of the rank columns must be named `rank_default` in order to compare influence of rank with batch effect.\nReformat results table.\n")}
        
        #clarification message
        batch = setdiff(tmp, "rank_default")
        cat("Detected",batch,"as batch variable for bias detection\n")

        #limit search for biased features to top 10k ranked genes
	bindev.results = filter(bindev.results, rank_default<=n.top)
        ## more complicated version to check for change in nSD dist with increasing rank
#	seq1 = seq(0,n.top, by=1000)
#        bindev.results$rank_default_bin = NA
#        for (i in 2:length(seq1)) {
#                bindev.results = mutate(bindev.results,
#                        rank_default_bin=if_else(between(rank_default, seq1[i-1]+1, seq1[i]),
#                                paste0("[",seq1[i-1]+1,",",seq1[i],"]"),
#                                rank_default_bin))
#        }
#        bindev.results = filter(bindev.results, !is.na(rank_default_bin))

        #determine SD based on all entries
        bindev.results$r.diff = bindev.results[,batch]-bindev.results[,"rank_default"]
        mean1 = mean(bindev.results$r.diff); sd1 = sd(bindev.results$r.diff)
        cat("\npooled mean rank diff.=",mean1);cat("\npooled SD rank diff.=",sd1,"\n")
        bindev.results$nSD = (bindev.results$r.diff-mean1)/sd1
        bindev.results$nSD.bin = cut(abs(bindev.results$nSD), right=FALSE,
                breaks=seq(0,max(bindev.results$nSD)+sd.interval, by=sd.interval), include.lowest=TRUE)
        #if working with bindev results that were looped over slides (each slide analyzed separately), also calculate per-slide SD
        if(length(grep("^slide$", colnames(bindev.results)))==1) {
                #determine per slide SD
                cat("\nper-slide mean and SD rank diff.\n")
                group_by(bindev.results, slide) %>% summarise(avg.r.diff=mean(r.diff), sd.r.diff=sd(r.diff))
                bindev.results = group_by(bindev.results, slide) %>% mutate(mean_slide=mean(r.diff), sd_slide=sd(r.diff), nSD_slide=(r.diff-mean_slide)/sd_slide) %>% ungroup()
                bindev.results$nSD.bin_slide = cut(abs(bindev.results$nSD_slide), right=FALSE,
                        breaks=seq(0,max(bindev.results$nSD_slide)+sd.interval, by=sd.interval), include.lowest=TRUE)
        }
        return(bindev.results)
}

findBiasedFeatures <- function(prepped.df, sd.safe=list(c("[0,5)"))) {
        #description
        ### prepped.df = dataframe produced by prepBias function
        ### sd.safe = character string of binned SD values that cover non-biased genes; all other values of binned SD will be considered biased
        ###### sd.safe is a list where the first element indicates the values to use for pooled SD approach and the second list element (if present) sets a different cutoff for the per-slide SD approach
        
        #error conditions 
        if(!is.list(sd.safe) & !is.character(sd.safe)) {stop("`sd.safe` must be character string specifying values of `nSD.bin` to mark safe from biased feature list")}
	if(length(sd.safe)>2) {stop("More than 2 list elements detected in `sd.safe`.\nReformat `sd.safe` so that the first list element is a character vector for pooled SD threshold and the second (optional) list element is a character vector for per-slide SD threshold.\n")} 
	check.missing = setdiff(unlist(sd.safe), unique(prepped.df$nSD.bin))
	if(length(check.missing)>0) {stop("One or more elements `sd.safe` are not valid `nSD.bin` values: ",check.missing,"\n")}
        
        tmp = grep("nSD\\.bin", colnames(prepped.df))
        if(length(tmp)==0) {stop("Supplied dataframe does not contain `nSD.bin` values. Make sure prepBias has been run.")}
        if(length(tmp)==1 & is.list(sd.safe) & length(sd.safe)==2) {
                cat("Only",tmp,"detected in supplied dataframe, but `sd.safe` expects `nSD.bin` and `nSD.bin_slide` to both be present.\n")
                stop("Either re-run prepBias to produce `nSD.bin_slide` or reformat `sd.safe` to be a single character vector.")
        }

        #formatting sd.safe and sd.safe messages
        if(length(sd.safe)==1) {
		if(length(tmp)==2) {
                        cat("Only 1 SD threshold supplied but detected pooled and per-slide SD approaches in results.\nApplying same threshold to both.\n")
                        cat("*** To identify biased features for pooled SD approach only, set second list element to NULL\n")
			cat("Pooled SD calculations: SD bins greater than",sd.safe[[1]],"will be considered biased.\n")
                        cat("Per-slide SD calculations: SD bins greater than",sd.safe[[1]],"will be considered biased.\n")
			if(is.character(sd.safe)) {cat("Coercing `sd.safe` to list...\n"); sd.safe = list(sd.safe, sd.safe)}
			else {sd.safe = list(sd.safe[[1]], sd.safe[[1]])}
		} 
		else {
			cat("Applying SD bin threshold to pooled SD calculation only.\n")
			cat("SD bins greater than",sd.safe[[1]],"will be considered biased.\n")
                        if(is.character(sd.safe)) {
				cat("Coercing `sd.safe` to list...\n")
				sd.safe = list(sd.safe)
			}
                }
        }
        else {
                if(is.list(sd.safe)) {
			if(sum(sapply(sd.safe, is.null))==0) {
				cat("Pooled SD calculations: SD bins greater than",sd.safe[[1]],"will be considered biased.\n")
				cat("Per-slide SD calculations: SD bins greater than",sd.safe[[2]],"will be considered biased.\n")
			}
			else {
				cat("NULL list element detected.\n")
				sd.safe = sd.safe[[1]]
				if(is.null(sd.safe)) {stop("NULL pooled SD threshold detected. Calculating per-slide bias only is not currently supported.")}
			}
		}
                if(is.character(sd.safe)) {
                        cat("Applying SD bin threshold to pooled SD calculation only.\n")
                        cat("SD bins greater than",sd.safe,"will be considered biased.\n")
                        cat("Coercing `sd.safe` to list...\n")
                        sd.safe = list(sd.safe)
                }
        }

        #label outliers
        prepped.df = ungroup(prepped.df)
        out.df = mutate(prepped.df, nSD.outlier=nSD.bin %in% setdiff(prepped.df$nSD.bin, sd.safe[[1]]))
        if(length(sd.safe)==2) {
                out.df = mutate(out.df, nSD.outlier_slide=nSD.bin_slide %in% setdiff(prepped.df$nSD.bin_slide, sd.safe[[2]]),
                        outlier.group= factor(paste(nSD.outlier, nSD.outlier_slide),
                                levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
                                labels=c("none","per-slide only","pooled only","both")
				)
                        )
        }
        return(out.df)

}

dotplotDF <- function(plot.genes, spe, norm.to.mbp=TRUE, order.by.rank=FALSE, rank.df=NULL, ...) {
	#description
	### plot.genes: list of ensembl IDs of genes for which expression will be extracted from spe object
	### spe: spatial experiment object with logcounts assay for plotting gene expression
	### norm.to.mbp: logical indicating whether or not to divide image avg gene expression by mbp
	### order.by.rank: logical indicating if the y axis of the resulting dotplot should be organized by genes with decreasing max rank (if FALSE, genes will be in default anti-alphabetical order)
	### rank.df: dataframe required if order.by.rank==TRUE, must provide bindev csv results here

	#error conditions
	if(order.by.rank==TRUE & is.null(rank.df)) {stop("order.by.rank = TRUE but no dataframe provided to pull rank information from")}
	if(length(grep("^ENSG",plot.genes))!=length(plot.genes)) {stop("ENSG names not detected in gene list provided. Concert from gene names to ensembl IDs")}
	
	#reminder messages of what mbp norm means for output
	if(norm.to.mbp==TRUE) {
		cat("norm.to.mbp=TRUE : avg expression values in raw counts divided by raw counts MBP\n")
		mbp.gene = rownames(spe)[rowData(spe)$gene_name=="MBP"]
		if(length(grep(mbp.gene, plot.genes))==0) plot.genes= c(plot.genes, mbp.gene)
	}

	if(norm.to.mbp==FALSE) cat("norm.to.mbp=FALSE : avg expression values in logNormCounts\n")

	#subset spe and build spe list
	spe2 = spe[plot.genes,]
	l4 = unique(spe2$sample_id)
	names(l4) = lapply(l4, function(x) unique(colData(spe2)[spe2$sample_id==x,"brain"]))
	l4 = lapply(l4, function(x) spe2[,colData(spe2)$sample_id==x])

	#build dotplot dataframe
	plot.genes.df = do.call(rbind, lapply(seq_along(l4), function(x) {
		if(norm.to.mbp==TRUE) {m1 = cbind("avg"=rowMeans(counts(l4[[x]])),
			"n"=rowSums(counts(l4[[x]])>0)/ncol(l4[[x]])) %>% as.data.frame() %>% tibble::rownames_to_column(var="gene")}
		if(norm.to.mbp==FALSE) {m1 = cbind("avg"=rowMeans(logcounts(l4[[x]])),
			"n"=rowSums(logcounts(l4[[x]])>0)/ncol(l4[[x]])) %>% as.data.frame() %>% tibble::rownames_to_column(var="gene")}
		
		m1$gene_name = rowData(l4[[x]])[m1$gene,"gene_name"]
		m1$brain = names(l4)[x]
		m1$position = unique(l4[[x]]$position)
		m1$slide = unique(l4[[x]]$slide)
		return(m1)
		})
	)
	
	#perform mbp norm if needed
	if(norm.to.mbp==TRUE) {
		mbp.df = filter(plot.genes.df, gene_name=="MBP")
		plot.genes.df = left_join(filter(plot.genes.df, gene_name!="MBP"), 
			mbp.df[,c("avg","brain","position","slide")], 
			by=c("brain", "position", "slide"), suffix=c("","_mbp")) %>% 
		mutate(expr_norm_mbp = avg/avg_mbp) %>% group_by(gene) %>% mutate(scaled_expr_norm_mbp= as.numeric(scale(expr_norm_mbp)))
	}
	else {
                plot.genes.df = group_by(plot.genes.df, gene) %>% mutate(scaled_avg= as.numeric(scale(avg)))
        }

	#add in rank information if requestioned
	if(order.by.rank==TRUE | !is.null(rank.df)) {
		#recode rank_batch so that this works with any batch variable
		## first isolate name of batch var
		tmp = setdiff(grep("rank_",colnames(rank.df), value=T), "rank_default")
		## then find column
		tmp.loc = grep(tmp, colnames(rank.df))
		## then rename
		colnames(rank.df)[tmp.loc] = "rank_batch"
		#now proceed
		best.rank.all.df = filter(rank.df, gene %in% plot.genes) %>% group_by(gene, gene_name) %>%
                        summarize(best.rank.all = min(rank_batch), .groups="drop") %>%
                        arrange(best.rank.all) %>% mutate(index=row_number(), i2=length(plot.genes)-index, ytext=paste(best.rank.all,gene_name, sep=" - ")) %>%
                        arrange(i2) %>% mutate(ylabel=factor(i2, levels=i2, labels=ytext))
		# best.rank.slide.df was only generated in case I wanted to see which slide was driving the result
		# i think i can comment this out with no downstream consequences
		#extra steps if working with per-slide loop
#		if(length(grep("^slide$", colnames(rank.df)))>0) {
#			best.rank.slide.df = filter(rank.df, gene %in% plot.genes) %>% group_by(slide, gene, gene_name) %>%
#				summarize(best.rank.slide=min(rank_batch), .groups="drop")
#			plot.genes.df <- left_join(plot.genes.df, best.rank.slide.df, by=c("slide","gene","gene_name")) %>%
#				left_join(best.rank.all.df, by=c("gene","gene_name")) %>%
#				mutate(xlabel=paste(position, brain))
#		}
		# if i do end up including above, i will need to add an else{ statement here
		#add rank information
		plot.genes.df <- left_join(plot.genes.df, best.rank.all.df, by=c("gene","gene_name")) %>%
			mutate(xlabel=paste(position, brain))
	}
	else {plot.genes.df <- mutate(plot.genes.df, xlabel=paste(position, brain), ylabel=gene_name)}
	return(plot.genes.df)
}

cols_7_legend <- function() {
	ggplot(cbind.data.frame("clus"=factor(c("L1","L2","L3","L5","L6","WM","bad"),
				levels=c("L1","L2","L3","L5","L6","WM","bad")), 
			"colpal"=cols_7), 
		aes(x="",y=11:5,fill=clus,label=clus))+
	geom_label(size=3, hjust=0)+
	ylim(0,15)+
	scale_fill_manual(values=cols_7)+
	theme_void()+theme(legend.position="none", plot.margin=unit(c(1,1,0,0), units="cm"))
}

cols_7o_legend <- function() {
	ggplot(cbind.data.frame("clus"=factor(c("L1","L1 (2)","L2","L3","L5/6","WM","bad"),
				levels=c("L1","L1 (2)","L2","L3","L5/6","WM","bad")), 
			"colpal"=cols_7_other), 
		aes(x="",y=11:5,fill=clus,label=clus))+
	geom_label(size=3, hjust=0)+
	ylim(0,15)+
	scale_fill_manual(values=cols_7_other)+
	theme_void()+theme(legend.position="none", plot.margin=unit(c(1,1,0,0), units="cm"))
}

cols_7o2_legend <- function() {
	ggplot(cbind.data.frame("clus"=factor(c("L1","L1/2/3","L2","L3",
	                                        "L5/6","WM","bad"),
				levels=c("L1","L1/2/3","L2","L3","L5/6","WM","bad")), 
			"colpal"=cols_7_other2), 
		aes(x="",y=11:5,fill=clus,label=clus))+
	geom_label(size=3, hjust=0)+
	ylim(0,15)+
	scale_fill_manual(values=cols_7_other2)+
	theme_void()+theme(legend.position="none", plot.margin=unit(c(1,1,0,0), units="cm"))
}

cols_8_legend <- function() {
	ggplot(cbind.data.frame("clus"=factor(c("L1","L1 (2)","L2","L2/3","L3","L5/6","WM","bad"),
				levels=c("L1","L1 (2)","L2","L2/3","L3","L5/6","WM","bad")), 
			"colpal"=cols_8), 
		aes(x="",y=seq(11.5,5.5, length.out=8),fill=clus,label=clus))+
	geom_label(size=3, hjust=0)+
	ylim(0,15)+
	scale_fill_manual(values=cols_8)+
	theme_void()+theme(legend.position="none", plot.margin=unit(c(1,1,0,0), units="cm"))
}

addSmallLegend <- function(myPlot, pointSize = 1, textSize = 6, spaceLegend = 0.1) {
    myPlot +
        guides(color = guide_legend(override.aes = list(size = pointSize, alpha=1))) +
        theme(legend.title = element_text(size = textSize), 
              legend.text  = element_text(size = textSize),
              legend.key.size = unit(spaceLegend, "lines"))
}
