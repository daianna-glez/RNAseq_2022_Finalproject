## R CODE

## Load recount3 R/Bioconductor package
library("recount3")
## Get all human projects available in recount3
human_projects <- available_projects()
## Find interest project
proj_info <- subset(
  human_projects,
  project == "SRP095512" & project_type == "data_sources"
)
## Create a RangedSummarizedExperiment (RSE) object
rse_gene_SRP095512 <- create_rse(proj_info)

## Convert raw counts to read counts 
assay(rse_gene_SRP095512, "counts") <- compute_read_counts(rse_gene_SRP095512)

## Check samples’ attributes
rse_gene_SRP095512$sra.sample_attributes
## [1] "cell type;;endothelial cell|disease state;;Healthy control|gender;;female|source_name;;dermal blood endothelial cell" 
## [2] "cell type;;endothelial cell|disease state;;Diabetic Patient|gender;;male|source_name;;dermal blood endothelial cell"  
## ...
##[10] "cell type;;endothelial cell|disease state;;Healthy control|gender;;female|source_name;;dermal blood endothelial cell" 

## Expand samples' attributes to access them
rse_gene_SRP095512 <- expand_sra_attributes(rse_gene_SRP095512)

## Check samples' information
colnames(colData(rse_gene_SRP095512))

## Check samples’ attributes
rse_gene_SRP095512$sra.sample_attributes
## [1] "cell type;;endothelial cell|disease state;;Healthy control|gender;;female|source_name;;dermal blood endothelial cell" 
## [2] "cell type;;endothelial cell|disease state;;Diabetic Patient|gender;;male|source_name;;dermal blood endothelial cell"  
## ...
##[10] "cell type;;endothelial cell|disease state;;Healthy control|gender;;female|source_name;;dermal blood endothelial cell" 

## Calculate gene assigned reads proportion for each sample
rse_gene_SRP095512$assigned_gene_prop <- rse_gene_SRP095512$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP095512$recount_qc.gene_fc_count_all.total
## Calculate the minimum proportion accepted: Median -3(Standard Deviation)
median(rse_gene_SRP095512$assigned_gene_prop)-3*sd(rse_gene_SRP095512$assigned_gene_prop)
## [1] 0.3239964

## General information of samples proportions
summary(rse_gene_SRP095512$assigned_gene_prop)
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
## 0.3665  0.4902  0.5349  0.5189  0.5608  0.6074

## Visualize graphically the frequency of samples proportions
hist(rse_gene_SRP095512$assigned_gene_prop)

## Samples with gene assigned proportion smaller than 0.32399
table(rse_gene_SRP095512$assigned_gene_prop < 0.32399)
## FALSE
## 10

## Gene assigned reads proportion for control and cases samples
tapply(rse_gene_SRP095512$assigned_gene_prop, rse_gene_SRP095512$sra_attribute.disease_state, summary)
## $`Diabetic Patient`
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.4544  0.5260  0.5503  0.5298  0.5541  0.5641 

## $`Healthy control`
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.3665  0.4902  0.5156  0.5116  0.5659  0.6074


## Boxplot of gene assigned reads proportion of controls and cases samples
library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP095512)), aes(y = assigned_gene_prop, x = sra_attribute.disease_state)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Group")

## Mean of expression levels of each gene
gene_means <- rowMeans(assay(rse_gene_SRP095512, "counts"))
summary(gene_means)
##  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##  0.0      0.0      0.6    277.5     35.2 884351.6 

## Number of genes with the first expression values
head(table(gene_means))
## gene_means
##     0   0.1   0.2   0.3   0.4   0.5 
## 22606  3700  2189  1403  1141   879

## Hold genes with mean expression greater than 0
rse_gene_SRP095512 <- rse_gene_SRP095512[gene_means >= 0.1, ]
## New dimensions
dim(rse_gene_SRP095512)
## [1] 41250    10

## Load the necessary library
library("edgeR")
## Convert RSE object into a DGElist object to analyze it through edgeR
dge <- DGEList(
  counts = assay(rse_gene_SRP095512, "counts"),
  genes = rowData(rse_gene_SRP095512)
)
## Calculate scaling factors to convert raw library sizes into effective library sizes.
dge <- calcNormFactors(dge)
## Check samples variables integrity

table(rse_gene_SRP095512$sra_attribute.gender)
## female   male 
##     8      2 
table(rse_gene_SRP095512$sra_attribute.cell_type)
## endothelial cell 
##              10 
table(rse_gene_SRP095512$sra_attribute.disease_state)
## Diabetic Patient  Healthy control 
##               4                6 
table(rse_gene_SRP095512$sra_attribute.source_name)
## dermal blood endothelial cell 
##                           10
## Create the model Y ~ disease state + gender + intercept
mod <- model.matrix(~ sra_attribute.disease_state + sra_attribute.gender,
                    data = colData(rse_gene_SRP095512)
)

## Matrix rank must equal the number of columns
qr(mod)$rank==ncol(mod)
## [1] TRUE

## Visualize the design of the matrix model 
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = colData(rse_gene_SRP095512),
  designFormula = ~ sra_attribute.disease_state + sra_attribute.gender,
  textSizeFitted = 4
)
cowplot::plot_grid(plotlist = vd$plotlist)

library("limma")
## Estimate and plot mean-variance relation
vGene <- voom(dge, mod, plot = TRUE)

## Linear model fit and t values for genes
eb_results <- eBayes(lmFit(vGene))
## All DE genes with higher t and p values for interest condition
de_results <- topTable(
  eb_results,
  ## Index of the interest coefficient (disease state)
  coef = 2,
  number = nrow(rse_gene_SRP095512),
  ## Conserve the original order of genes
  sort.by = "none")
dim(de_results)
## [1] 41250    16

## MA plot takes gene expression and disease state 
limma::plotMA(eb_results, coef = 2)
## Hold only DE genes
limma::plotMA(eb_results[which(de_results$adj.P.Val < 0.05),], coef=2)

## Volcano plot showing the top 3 genes with lowest p values 
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name, cex=0.5, hl.col = "orchid4")
## Hold only DE genes
volcanoplot(eb_results[which(de_results$adj.P.Val < 0.05),], coef = 2, highlight = 3, names = de_results[which(de_results$adj.P.Val < 0.05), "gene_name"], cex=0.5, hl.col = "orchid4", ylim=c(0,7))


## Expression of DE genes along the 10 samples
exprs_heatmap <- vGene$E[de_results$adj.P.Val <= 0.05, ]
## DE genes in de_results and vGene are in the same order
identical(rownames(exprs_heatmap), de_results[de_results$adj.P.Val <=0.05, "gene_id"])

## DE genes names as row names of heatmap
rownames(exprs_heatmap) <- de_results[de_results$adj.P.Val <=0.05, "gene_name"]
## Data frame with interest variables of samples
col_df <- as.data.frame(colData(rse_gene_SRP095512)[,c("sra_attribute.disease_state","sra_attribute.gender")])
colnames(col_df) <- c("Disease state", "Gender")

## Heatmap with clustered rows and columns
library("pheatmap")
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  fontsize_row = 5,
  show_colnames = FALSE,
  annotation_col = col_df
)

## Reproduce code in 120 chars
options(width = 120)
sessioninfo::session_info()
