library(tidyverse)
library(limma)
library(edgeR)
library(WGCNA)
library(ggplot2)
library(WGCNA)
library(ggplot2)
library(igraph)
# BiocManager::install("pSI")
library(pSI)
library(nlme)
library(gridExtra)
library(plyr)
enableWGCNAThreads()
allowWGCNAThreads()

dir.create("./Module_Plots")

load("./P35_CX_norm_counts_bicor20_DeepSplit_0sft_power18.RData")
geneTree = networks$tree
datExpr=networks$datExpr
merged = networks$merged
modules = merged$colors
genes=rownames(datExpr)
MEs = networks$MEs
kMEtable = networks$kMEtable
CX_meta=readxl::read_xlsx("./meta_group_P35.xlsx", sheet = 1)
datMeta=CX_meta[-c(9:10),]
datMeta <- datMeta %>% arrange(SampleGroup)
datMeta$SampleGroup=factor(datMeta$SampleGroup,levels = c("WT","Cul3 HET"))


modTrait=data.frame()
for(i in 2:length(unique(modules))) {
  me = MEs$eigengenes[,i]
  moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
  moduleColor = labels2colors(moduleNumber)
  s = summary(lm(me ~ SampleGroup + Gender,data=datMeta))$coefficients
  for(grp in c("Cul3 HET")){
    
    rowID = paste0("SampleGroup", grp)
    
    modTrait = rbind.fill(modTrait,
                          data.frame(Module=moduleColor, moduleNumber= moduleNumber, Group=grp,
                                     beta = s[rowID, "Estimate"], SE = s[rowID, "Std. Error"], t=s[rowID, "t value"], p=s[rowID, "Pr(>|t|)"]))}
  
  modTrait$fdr=p.adjust(modTrait$p,method = "fdr")
  modTrait$signedLog10fdr = -log10(modTrait$fdr) * sign(modTrait$beta)
  modTrait$signedLog10fdr[modTrait$fdr > .05] = 0
  modTrait$text = signif(modTrait$beta, 1)
  print(modTrait)
  
}
#writexl::write_xlsx(modTrait,"Module_Trait_E17_5_CX.xlsx")

sft <- 18
significant_modules <- c()
for (i in seq(length(readxl::excel_sheets("./P35_GOEnrich_KME_sft18_signed_norm_counts_bicor_all_moderate_0_MM_20.xlsx")))){
  if (nrow(readxl::read_xlsx("./P35_GOEnrich_KME_sft18_signed_norm_counts_bicor_all_moderate_0_MM_20.xlsx",sheet = i)) > 0){
    significant_modules <- c(significant_modules, i+1)
  }
  
}


module_colors <- list()
module_numbers <- list()
module_genes <- list()
for(i in significant_modules) {
  
  print(i)
  me = MEs$eigengenes[,i]
  moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
  moduleColor = labels2colors(moduleNumber)
  module_colors <- append(module_colors,moduleColor)
  module_numbers <- append(module_numbers,moduleNumber)
  moduleGenes = genes[modules==moduleNumber]
  
  pdf(file=paste0("./Final_sec_round/", moduleNumber, moduleColor, "_P35_CX_WPCNA.pdf"), width=20, height=11, useDingbats = FALSE)
  s = summary(lm(me ~ SampleGroup,data=datMeta))$coefficients
  dat2=data.frame()
  for(grp in c("Cul3 HET")) {
    rowID = paste0("SampleGroup", grp)
    dat2 = rbind.fill(dat2,
                 data.frame(Module=moduleColor, moduleNumber= moduleNumber, Group=grp,
                            beta = s[rowID, "Estimate"], SE = s[rowID, "Std. Error"], t=s[rowID, "t value"], p=s[rowID, "Pr(>|t|)"]))
  }
  dat2$fdr = modTrait$fdr[modTrait$moduleNumber == moduleNumber]
  dat2$p.symbol = ""; dat2$p.symbol[dat2$fdr<.05] ="*"
  g1=ggplot(dat2, aes(x=Group, y=beta, label=p.symbol)) + geom_bar(stat="identity",fill=moduleColor) +  
    geom_errorbar(aes(ymin=(beta - SE), ymax=(beta + SE)), width=0.25,size=0.25) +
    geom_text(color="red",size=8,aes(y=beta+ 1.3*sign(beta)*SE ))+
    ggtitle("Module-Trait Association")+
    xlab("")+
    ylab("Linear Regression Beta")
  
  dat=data.frame(Eigengene=me, Gender=datMeta$Gender,Group=datMeta$SampleGroup)
  g2=ggplot(dat,aes(x=Group,y=Eigengene,color=Group)) + geom_point(size=5)+    
    scale_color_manual(values = c("black","red","blue"))
  
  go=readxl::read_xlsx("./P35_GOEnrich_KME_sft18_signed_norm_counts_bicor_all_moderate_0_MM_20.xlsx", sheet = i-1)
  go=go%>%
    top_n(15)
  g3=ggplot(go, aes(x=reorder(term.name, -log10(p.value)), y=-log10(p.value))) + geom_bar(stat="identity", fill="royalblue") + 
    coord_flip() + geom_hline(yintercept=-log10(0.05), lty=2, color="red")
  
  hubGenes = moduleGenes[order(kMEtable[moduleGenes,i], decreasing = T)[1:min(20,length(moduleGenes))]]
  annotation<-readxl::read_xlsx("./P35_CX_combat_norm_data_log.xlsx", sheet = 1)
  annotation<-as.data.frame(annotation)
  hubGene.symbols=annotation[match(hubGenes, annotation$Protein),]$Gene
  
  gene_symbols <- annotation[match(moduleGenes,annotation$Protein),]$Gene
  
  module_genes <- append(module_genes,list(data.frame("Protein_ID" = moduleGenes,
                                                      "Gene_name" = gene_symbols)))
  
  adjMat = adjacency(t(datExpr[hubGenes,]),type = "signed",corFnc = "cor", power=sft)
  adjMat[adjMat < quantile(adjMat,0.1)]=0
  graph <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=T,diag=F)
  plotcord= data.frame(layout_with_fr(graph))
  colnames(plotcord) = c("X1","X2")
  edgelist <- get.edgelist(graph,names = F)
  edges <- data.frame(plotcord[edgelist[,1],], plotcord[edgelist[,2],])
  colnames(edges) <- c("X1","Y1","X2","Y2")
  plotcord = cbind(plotcord, data.frame(gene=hubGene.symbols))
  
  g4=ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), 
                             data=edges, 
                             size = 0.5, colour="grey") + 
    geom_point(aes(X1, X2), data=plotcord,color=moduleColor,size=4) + geom_text(aes(x=X1, y=X2+.2, label=gene),fontface="bold",size=6,data=plotcord) +
    theme_classic() +
    theme_void()
  
  
  grid.arrange(grobs=list(g1,g2,g3,g4), layout_matrix=rbind(c(1,1,2,2,2,5,5),c(3,3,3,4,4,4,4)),
               top=paste0("Module ", moduleNumber, " (", moduleColor, ")"),padding=unit(2, "line"))
               dev.off()
}

# module_genes <- set_names(module_genes, paste0("M",module_numbers))
# writexl::write_xlsx(module_genes, "./module_genes_P7_CB_num.xlsx")


