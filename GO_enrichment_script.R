#Libraries####
BiocManager::install("org.Hs.eg.db")
BiocManager::install("UniProt.ws")
BiocManager::install("enrichplot")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(GOSemSim)
library(enrichplot)

#Query protein list file opening####
human<-read.csv2("ChIRP.csv", header= TRUE)

#Uniprot accession numbers conversion to ENTREZID####

#H. sapiens proteome loading
library(UniProt.ws)
up <- UniProt.ws(taxId=9606)

#Conversion
ids_human <- select(up,
                    keys=human,
                    columns=c("xref_geneid","UniProtKB"),
                    keytype="UniProtKB")

#GeneID trimming
entrez_genes_human <- ids_human$GeneID
entrez_genes_human <- entrez_genes_human[!is.na(entrez_genes_human)]
entrez_genes_human <- gsub(";", "", entrez_genes_human)
entrez_genes_human <- trimws(entrez_genes_human)  

#GO term enrichment with enrichGO()####
ego_BP <- enrichGO(gene          = entrez_genes_human,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ego_MF <- enrichGO(gene          = entrez_genes_human,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ego_CC <- enrichGO(gene          = entrez_genes_human,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

#Enrichment tables export####
write.table(ego_BP, file = "BP_human_clusterProfiler.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ego_MF, file = "MF_human_clusterProfiler.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ego_CC, file = "CC_human_clusterProfiler.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

#Retrieval of all level 4 (and higher) GO terms containing the query genes with groupGO()####
ggo_BP <- groupGO(gene           = entrez_genes_human,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENTREZID',
                  ont           = "BP",
                  level         = 4,
                  readable      = TRUE)

ggo_MF <- groupGO(gene           = entrez_genes_human,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENTREZID',
                  ont           = "MF",
                  level         = 4,
                  readable      = TRUE)

ggo_CC <- groupGO(gene           = entrez_genes_human,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENTREZID',
                  ont           = "CC",
                  level         = 4,
                  readable      = TRUE)

#Count == 0 filtering to delete empty GO terms (containing none of the query genes)
ggo_BP_filtered <- ggo_BP@result[ggo_BP@result$Count > 0, ]
ggo_MF_filtered <- ggo_MF@result[ggo_MF@result$Count > 0, ]
ggo_CC_filtered <- ggo_CC@result[ggo_CC@result$Count > 0, ]

#"All GO terms" tables export####
write.table(ggo_BP_filtered, file = "BP_human_clusterProfiler_all.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ggo_MF_filtered, file = "MF_human_clusterProfiler_all.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ggo_CC_filtered, file = "CC_human_clusterProfiler_all.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

#Semantic similarity reduction####

#Data frame construction
d <- GOSemSim::godata('org.Hs.eg.db', ont="BP")
e <- GOSemSim::godata('org.Hs.eg.db', ont="MF")
f <- GOSemSim::godata('org.Hs.eg.db', ont="CC")

#Reduction
ego_bp2 <- clusterProfiler::simplify(ego_BP, cutoff=0.4, by="p.adjust", select_fun=min)
ego_bp_reduced <- pairwise_termsim(ego_bp2, method = "Wang", semData = d)

ego_mf2 <- clusterProfiler::simplify(ego_MF, cutoff=0.5, by="p.adjust", select_fun=min)
ego_mf_reduced <- pairwise_termsim(ego_mf2, method = "Wang", semData = e)

ego_cc2 <- clusterProfiler::simplify(ego_CC, cutoff=0.6, by="p.adjust", select_fun=min)
ego_cc_reduced <- pairwise_termsim(ego_cc2, method = "Wang", semData = f)

#"Reduced GO terms" tables export####
write.table(ego_bp_reduced, file = "BP_reduced_human_clusterProfiler.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ego_mf_reduced, file = "MF_reduced_human_clusterProfiler.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ego_cc_reduced, file = "CC_reduced_human_clusterProfiler.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

#BP dendrogram construction####
p <- treeplot(ego_bp_reduced,
              fontsize_tiplab = 4.5,
              fontsize_cladelab = 4.5,
              showCategory = nrow(ego_bp_reduced@result), 
              tiplab_offset = rel(0.226),   #tunes the offset of labels from the points
              cladelab_offset = rel(12.3))  #tunes the offset of clade bars and labels from the points

p$layers[[4]]$aes_params$linewidth <- 1.05  #tunes the thickness of clade bars
p$layers[[3]]$data$x <- 18.6                #tunes the offset of clade labels from the clade bars

p + theme(
  legend.position   = "bottom", 
  legend.box = "horizontal",    
  legend.spacing.x   = unit(2, "cm"),   
  legend.box.spacing = unit(0, "cm")    
)+ guides(
  size = guide_legend( 
    order = 1,         
    title = "Count",
    theme = theme(
      legend.text = element_text(size = 10),
      legend.title = element_text(
        size=13,
        margin = margin(t=-3, r=3)
      ))),
  colour = guide_colorbar( 
    order = 2,
    title = "Adjusted pvalue",
    barwidth = unit(6,"cm"), 
    theme = theme(
      legend.text = element_text(size = 10), 
      legend.title = element_text(
        size=13,                               
        margin = margin(b=5, t=-10, r=5))      
    )),
  colour_ggnewscale_1 = "none" 
)+
  xlim(0,23) #tunes the branches length of the dendrogram 

ggsave("BP_human_clusterprofiler.png", dpi = 800)
dev.off()


#MF dendrogram construction####
p <- treeplot(ego_mf_reduced,
              fontsize_tiplab = 4.5,
              fontsize_cladelab = 4.5,
              showCategory = nrow(ego_mf_reduced@result), 
              tiplab_offset = rel(0.226),   
              cladelab_offset = rel(9.6))   

p$layers[[4]]$aes_params$linewidth <- 1.05  
p$layers[[3]]$data$x <- 15.86              

p + theme(
  legend.position   = "bottom", 
  legend.box = "horizontal",    
  legend.spacing.x   = unit(2, "cm"),   
  legend.box.spacing = unit(0, "cm")    
)+ guides(
  size = guide_legend( 
    order = 1,         
    title = "Count",
    theme = theme(
      legend.text = element_text(size = 10),
      legend.title = element_text(
        size=13,
        margin = margin(t=-3, r=3)
      ))),
  colour = guide_colorbar( 
    order = 2,
    title = "Adjusted pvalue",
    barwidth = unit(6,"cm"), 
    theme = theme(
      legend.text = element_text(size = 10), 
      legend.title = element_text(
        size=13,                               
        margin = margin(b=5, t=-10, r=5))      
    )),
  colour_ggnewscale_1 = "none" 
)+
  xlim(0,23)  

ggsave("MF_human_clusterprofiler.png", dpi = 800)
dev.off()

#CC dendrogram construction####
p <- treeplot(ego_cc_reduced,
              fontsize_tiplab = 4.5,
              fontsize_cladelab = 4.5,
              showCategory = nrow(ego_cc_reduced@result), 
              tiplab_offset = rel(0.29),   
              cladelab_offset = rel(8.4))   

p$layers[[4]]$aes_params$linewidth <- 1.05  
p$layers[[3]]$data$x <- 16.75                 

p + theme(
  legend.position   = "bottom", 
  legend.box = "horizontal",    
  legend.spacing.x   = unit(2, "cm"),   
  legend.box.spacing = unit(0, "cm")    
)+ guides(
  size = guide_legend( 
    order = 1,         
    title = "Count",
    theme = theme(
      legend.text = element_text(size = 10),
      legend.title = element_text(
        size=13,
        margin = margin(t=-3, r=3)
      ))),
  colour = guide_colorbar( 
    order = 2,
    title = "Adjusted pvalue",
    barwidth = unit(6,"cm"), 
    theme = theme(
      legend.text = element_text(size = 10), 
      legend.title = element_text(
        size=13,                               
        margin = margin(b=5, t=-10, r=5))      
    )),
  colour_ggnewscale_1 = "none"
)+
  xlim(0,28) 

ggsave("CC_human_clusterprofiler.png", dpi = 800)
dev.off()
