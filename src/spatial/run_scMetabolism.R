library(scMetabolism)
sc.metabolism <- function(countexp, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG") {
  signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolism")
  signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", package = "scMetabolism")

  if (metabolism.type == "KEGG")  {gmtFile<-signatures_KEGG_metab; cat("Your choice is: KEGG\n")}
  if (metabolism.type == "REACTOME")  {gmtFile<-signatures_REACTOME_metab; cat("Your choice is: REACTOME\n")}

  #imputation
  if (imputation == F) {
    countexp2<-countexp
  }
  if (imputation == T) {
    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]; row.names(countexp2) <- row.names(countexp)
  }

  #signature method
  cat("Start quantify the metabolism activity...\n")

  #VISION
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2) / n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile,projection_methods = "UMAP")
    options(mc.cores = ncores)
    vis <- analyze(vis)
    signature_exp<-data.frame(t(vis@SigScores))
  }
  
  #AUCell
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), nCores=ncores, plotStats=F) #rank
    geneSets <- getGmt(gmtFile) #signature read
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #calc
    signature_exp <- data.frame(getAUC(cells_AUC))
  }
  #ssGSEA
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }
  #GSVA
  if (method == "GSVA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("gsva"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }
  signature_exp
}
# Load single cell transcriptomic data (e.g., Mouse Liver, Melanoma or HNSCC)
# The loaded object "single cell data" should be a matrix with genes in rows and samples in columns.
load("/single_cell.rda")

# Infer metabolic pathway activity using scMetabolism
# method can be one of: "AUCell", "GSVA", "ssGSEA", or "VISION"
metabolism.matrix <- sc.metabolism(
    countexp = single_cell,
    method = "AUCell",          # change as needed
    imputation = F,
    ncores = 2,
    metabolism.type = "KEGG"
)
