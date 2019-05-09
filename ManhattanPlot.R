manhattan <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", xaxis.cex=1, limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), annotate=NULL, ...) {
 
    d=dataframe
    
    #throws error if you don't have columns named CHR, BP, and P in your data frame.
	if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
    
	# limits chromosomes to plot. (23=x, 24=y, 25=par?, 26=mito?)
    if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
	
	 # remove na's, sort by CHR and BP, and keep snps where 0<P<=1
    d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1))
	
	# -log10(p-value)
	d$logp = -log10(d$P)
	
	# sets colors based on colors argument.
    colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
    
	# sets the maximum value on the y axis (on the -log10p scale).
    if (ymax=="max") ymax<-ceiling(max(d$logp))
    if (ymax<8) ymax<-8
    
	# creates continuous position markers for x axis for entire chromosome. also creates tick points.
    d$pos=NA
    ticks=NULL
	lastbase=0
    numchroms=length(unique(d$CHR))
    if (numchroms==1) {
        d$pos=d$BP
        ticks=floor(length(d$pos))/2+1
    } else {
        for (i in unique(d$CHR)) {
        	if (i==1) {
    			d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
    		} else {
    			lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
    			d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
    		}
    		ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    	}
    }
    
	# create the plot
    if (numchroms==1) {
		# if you only have a single chromosome, the x axis is the chromosomal position
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
    }	else {
		# if you have multiple chromosomes, first make the plot with no x-axis (xaxt="n")
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", cex=0.3, ...))
		# then make an axis that has chromosome number instead of position
        axis(1, at=ticks, lab=unique(d$CHR), cex.axis=xaxis.cex)
        icol=1
        for (i in unique(d$CHR)) {
            with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], cex=0.3, ...))
            icol=icol+1
    	}
    }
    
	# create a new data frame with rows from the original data frame where SNP is in annotate character vector.
	# then plot those points over the original graph, but with a larger point size and a different color.
    if (!is.null(annotate)) {
        d.annotate=d[which(d$SNP %in% annotate), ]
        with(d.annotate, points(pos, logp, col="red", cex=0.5, ...))
    }
    
	# add threshold lines
    if (suggestiveline) abline(h=suggestiveline, col="blue")
    if (genomewideline) abline(h=genomewideline, col="red")
}
