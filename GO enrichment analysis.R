
library(topGO)
library('org.Hs.eg.db')

###Matrices of differentiallly expressed genes were created using the script 'Participant Clustering and DEG'

#Create lists of universal gene IDs for Gene Ontology enrichment analysis
up_genes_SARP<-list(rownames(pc1.up.genes), rownames(pc2.up.genes), rownames(pc3.up.genes), rownames(pc4.up.genes))
down_genes_SARP<-list(rownames(pc1.down.genes), rownames(pc2.down.genes), rownames(pc3.down.genes), rownames(pc4.down.genes))

#Set node_num for the number of included terms in bar plotting
#For further information on TopGO, please visit: https://bioconductor.org/packages/release/bioc/html/topGO.html

node_num<-10
mypal <- colorRampPalette( c( "darkblue", "dodgerblue", "white", "red" ) )( 200 )

for(n in 1:length(up_genes_SARP)){
  
  if(length(up_genes_SARP[[n]]>1)){
    GOI<-unique(mapIds(org.Hs.eg.db, up_genes_SARP[[n]], 'ENTREZID', 'SYMBOL'))
    background_names<-unique(mapIds(org.Hs.eg.db, rownames(SARP_BEC_data), 'ENTREZID', 'SYMBOL'))
    geneList <- factor(as.integer(background_names %in% GOI))
    names(geneList) <- background_names
    GOdata <- new("topGOdata", allGenes = geneList,
                  nodeSize = 25, ontology = "BP", annot = annFUN.org, mapping = "org.Hs.eg.db")
    resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
    upRes <- GenTable(GOdata, classicFisher = resultFisher, ranksOf = "classicFisher", topNodes = node_num, numChar = 500)
    upRes$classicFisher[which(upRes$classicFisher=="<1e-30")]<-1*10^-30
    GOI<-unique(mapIds(org.Hs.eg.db, down_genes_SARP[[n]], 'ENTREZID', 'SYMBOL'))
    background_names<-unique(mapIds(org.Hs.eg.db, rownames(SARP_BEC_data), 'ENTREZID', 'SYMBOL'))
    geneList <- factor(as.integer(background_names %in% GOI))
    names(geneList) <- background_names
    GOdata <- new("topGOdata", allGenes = geneList,
                  nodeSize = 25, ontology = "BP", annot = annFUN.org, mapping = "org.Hs.eg.db")
    resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
    downRes <- GenTable(GOdata, classicFisher = resultFisher, ranksOf = "classicFisher", topNodes = node_num, numChar = 500)
    downRes$classicFisher[which(downRes$classicFisher=="<1e-30")]<-1*10^-30
    
    upRes_order<-upRes[order(-log10(as.numeric(upRes$classicFisher)), decreasing = F),]
    allRes<-rbind(upRes_order, downRes)
    plot_data<- -log10(as.numeric(allRes$classicFisher))
    plot_data[(length(plot_data)-node_num+1):length(plot_data)]<-plot_data[(length(plot_data)-node_num+1):length(plot_data)]*-1
    
    png(paste0("GO_barplot_SARP_cluster_", n, ".png"), height = 485*5, width = 1100*5, res = 300, pointsize = 16)
    par(mar=c(5.1,40.1,4.1,2.1))
    barplot(plot_data, horiz = T, axes = F, col = map2color(plot_data, mypal),
            names.arg = allRes$Term, axisnames = T, las = 1, cex.axis = 1.25)
    abline(v = 0, lty = 1, col = "black")
    dev.off()
    
  }
}

