#-------------------------------------------------------------------------------#
# 
# Input parameters to set (n=6)
# 1. hemisphereParameter: 'lh' (left hemisphere) or 'rh' (right hemisphere), 
# 2. ReferenceGeneListFile: Path to the file contatining Gene Symbols, correlation 
#                           between BrainSpan and Allen expression levels and Zeisel 
#                           Cell types for the genes 
#    remaining after 2-step QC. 
# 3. wd: working directory where the downloaded files are loacated.
# 4. nreps: number of re-sampling replications for association test
# 5. fig_cex.main: for graphical display: scale of the title for each panel showing 
#                  the empirical density of cell-type-specific expression-phenotype  
#                  correlation coefficients.
# 6. fig_ylim: NULL or c(min,max), where min and max correspnd to minimum and maximum
#              values of y-axis. If NULL, each panel uses its local minimum and maximum
#			   estimated density values for each cell type; consequently, the y-axis scale 
#              can vary by cell type. Otherwise, all the panel will use the same y-axis 
#              scale specified by c(min,max).
# **** Example: ****
# hemisphereParameter = "lh"
# UserProfileFile = 'Profile_males_corcoef_thickness_age.csv'
# wd = "/Users/jshinb/Downloads/4752955/" 
# nreps = 10000
# siglevel = 0.05
# fig_cex.main = 1  
# fig_ylim = c(0,1.6) # set a global y-axis scale
#
# Notes
# UserPhenoFile must be clneaned and formatted before loading:
# * The file must be saved as a csv.
# * The file must contain a column vector contatining 34 cortical-region-specific values 
#   with rownames containing the 34 regions.
# * It is assumed that these values are based on QC-ed data set 
#   (without the MRI-QC-failed, outlier measurements).  
# 
# Outputs (n=2): two files containing the folloiwng plot and table are created in 
#                'output' directory 
# 1. 3-row-by-3-column Figure showing empirical distributions for expression-phenotype 
#    correlation coefficients and the correponding (1-'siglevel')*100% critical values 
#    under null hypothesis of no association cell-type phenotype profile obtained from 
#    re-sampling distributions with 'nreps' replications
# 2. Re-sampling based test results, where the test statistic is mean expression-phenotype
#    coefficients (cell-type-specific)
#-------------------------------------------------------------------------------#

#------------------------- USER INPUTS #1 - #6 --------------------------#
# Users can complete this part and source the file (see the **example** above)
# To source file, Open R (or R studio) and type the following command line:
# source('CreateCorrelationCoefficients.r', echo=TRUE)
#
hemisphereParameter = "lh"
UserProfileFile = 'Profile_males_corcoef_thickness_age.csv'
wd = "/Users/jshinb/Downloads/4752955/"
nreps = 10000
siglevel = 0.05
fig_cex.main = 1
fig_ylim = c(0,1.6) # set a global y-axis scale
#------------------------------------------------------------------------#

## Load necessary libraries 
## 
if("stringr" %in% rownames(installed.packages()) == FALSE) {install.packages("stringr")}
if("scales" %in% rownames(installed.packages()) == FALSE) {install.packages("scales")}

require(stringr)
require(scales)

## Set working directory to where the files are downloaded
setwd(wd) 
dir.create('output')

## Load and format profile data (from the user)
ctx.pheno.profile = read.csv(UserProfileFile,stringsAsFactors = F,row.names = 1)
if(all(str_detect(rownames(ctx.pheno.profile),"ctx.lh") | str_detect(rownames(ctx.pheno.profile),"ctx.rh")) ){
  ctx.pheno.profile = ctx.pheno.profile[str_detect(rownames(ctx.pheno.profile),
                                                   paste("ctx.",hemisphereParameter,sep="")),,drop=F]
}else{# assume it has the roi names of the correct hemisphere
  rownames(ctx.pheno.profile) <- paste("ctx",hemisphereParameter,rownames(ctx.pheno.profile),sep=".")
}

## Load and format rownames (i.e., cortical region names) for the Allen Gene expression data
geneExpressionMatrix <- read.table("AllenHBA_DK_ExpressionMatrix.tsv")
regionMapping <- read.table("DKRegionStatistics.tsv")
regionMapping <- subset(regionMapping, Hemisphere == hemisphereParameter)
rownames(regionMapping) <- gsub("-", ".", rownames(regionMapping))

## Load gene symbols and other information for the reference panel 
refpanel <- read.csv('Reference_Consistent_Genes_ObtainedBy2StageFiltering.tsv',
                     sep="\t",stringsAsFactors = F)

## Extract gene expression levels for the 2511 genes in refpanel
expressionForGene <- geneExpressionMatrix[refpanel$GeneSymbol,rownames(regionMapping)]
expressionForGene <- t(expressionForGene)
expressionForGene.df <- data.frame(expressionForGene,stringsAsFactors = F)
names(expressionForGene.df) <- colnames(expressionForGene)
expressionForGene <- expressionForGene.df
rm(expressionForGene.df)

## Match the column names between expression and cortical phenotype data
if(sum(rownames(expressionForGene) == rownames(ctx.pheno.profile)) < 34){
  o <- match(rownames(expressionForGene),rownames(ctx.pheno.profile))
  ctx.pheno.profile <- ctx.pheno.profile[o,,drop=F] 
  if(sum(rownames(expressionForGene) == rownames(ctx.pheno.profile))){
    cat("The rows are matched with respect to the cortical regions.\n")
  }
}else if( sum(rownames(expressionForGene) == rownames(ctx.pheno.profile)) == 34 ){
  cat("The rows between Allen expression and cortical phenotype profile data are already matched with respect to the cortical regions - no need to do anything.\n")
}

## Calculate correlation coefficients between gene-expression and cortical-phenotype profiles
cor_ReferencePanel <- NULL
for(i in 1:ncol(expressionForGene)){
  cor_ReferencePanel <- c(cor_ReferencePanel,cor(expressionForGene[,i],ctx.pheno.profile$x,method="pearson"))
}
cor_ReferencePanel <- data.frame(Gene=names(expressionForGene),
                                 expression_phenotype_r = cor_ReferencePanel,
                                 stringsAsFactors = F)
cor_ReferencePanel.temp <- merge(cor_ReferencePanel,refpanel,by.x="Gene",by.y="GeneSymbol")
cor_ReferencePanel <- cor_ReferencePanel.temp
rm(cor_ReferencePanel.temp)

## Test association between cell-type and profile-correlation and create Figure and Table
## generate the names of the empirical density plot and the re-sampling test result table files
if(str_sub(UserProfileFile,(str_length(UserProfileFile)-3),str_length(UserProfileFile))==".csv"){
  CorPlotFile = str_replace(UserProfileFile,'.csv','.png')
  CorPlotFile = paste('output/',CorPlotFile,sep="")
}else{
  stop("\'UserProfileFile\' must be a csv file. Please format it as instructed.\n")
}
tab_filename = paste('output/Table_resampling_res_with_',nreps,'_replications_',
                     UserProfileFile,sep="")

## plotting begins
png(CorPlotFile,width=10,height=8.5,pointsize=17,bg="white",units="in",res=250)
par(mfrow=c(3,3),mar=c( 2.6, 3.6, 2.6, 1.1),oma=c(2,0,0,0),lwd=3)

celltypes = sort(unique(cor_ReferencePanel$CellType)[!is.na(unique(cor_ReferencePanel$CellType))])
r_resampling <- vector(mode = "list",length=length(celltypes))
r.T.celli <- r.T.celli.p <- NULL
for(i in 1:9){
  celltype = celltypes[i]
  gene_ind = !is.na(cor_ReferencePanel$CellType) & cor_ReferencePanel$CellType==celltype
  celltype_genes = cor_ReferencePanel$Gene[gene_ind] 
  
  # print the range of donor to median correlation
  range.r = round(range(geneExpressionMatrix[celltype_genes,"Average.donor.correlation.to.median"]),2)
  cat(paste(celltype," - Average correlation between donors and median expression: (",range.r[1],",",range.r[2],")\n",
            sep=""))
  
  ## perform re-sampling test (2-sided)
  size.i=sum(gene_ind)
  r_resampling[[i]] <- matrix(NA,nrow=size.i,ncol=nreps)
  x=1:length(cor_ReferencePanel$Gene)
  for(j in 1:nreps){
    row.ind = sample(x,size = size.i,replace = F)
    r_resampling[[i]][,j] <- cor_ReferencePanel$expression_phenotype_r[row.ind]
  }
  names(r_resampling)[i] <- celltype
  
  ## test statistics (T): mean correlation coefficient
  avgr = mean(cor_ReferencePanel$expression_phenotype_r[cor_ReferencePanel$Gene %in% celltype_genes])
  avgr_reps = apply(r_resampling[[i]],2,mean)
  avgr_reps = c(avgr,avgr_reps)
  p = sum(abs(avgr_reps) >= abs(avgr))/length(avgr_reps)
  r.T.celli <- c(r.T.celli,avgr)
  r.T.celli.p <- c(r.T.celli.p, p)
  
  ## plotting empirical density of expression-phenotype correlation coefficients for each cell type
  ntest=1
  unAdjCI = quantile(avgr_reps,probs = c((siglevel/2)/ntest,(1-(siglevel/2)/ntest)))
  res_r=cor_ReferencePanel$expression_phenotype_r[gene_ind]
  d_res=density(res_r)
  xlim=range(density(cor_ReferencePanel$expression_phenotype_r)$x)
  plot(d_res,main=celltype,las=1,cex.main=fig_cex.main,
       xlab="",yaxs="i",ylab="",xlim=xlim,ylim=fig_ylim)
  abline(v=mean(res_r),col="black",lwd=3,lty=2) ## fixed
  myCI = unAdjCI
  if(is.null(fig_ylim)){
    polygon(c(myCI,rev(myCI)),c(c(0,0),c(max(d_res$y),max(d_res$y))),
            col=scales::alpha("grey",0.5),border=NA)
  }else{#!is.null(ylim)
    polygon(c(myCI,rev(myCI)),c(c(0,0),c(max(fig_ylim),max(fig_ylim))),
            col=alpha("grey",0.5),border=NA)
  }
  
}
title(xlab="Correlation coefficient (r)",outer=T,line=0,cex.lab=1.25)
title(ylab="Density",outer=T,line=-1,cex=0.9,cex.lab=1.25)
dev.off()

## Construct the result table
res_tab <- data.frame(Cell_Type=celltypes,
                      nGenes = sapply(r_resampling,nrow),
                      avgr = round(r.T.celli,3),
                      p = r.T.celli.p,
                      fdr.p = p.adjust(r.T.celli.p,method="fdr"), 
                      stringsAsFactors=F)
print(res_tab)

## write out the result table to 'res_tab,tab_filename'
write.csv(res_tab,tab_filename,quote=F,row.names=F)
