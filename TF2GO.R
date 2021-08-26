###Transcription Factor to Gene Ontology Term Mapping
###Accepts an input gene list as Universal Gene Symbols and creates a network 
###connectivity map between the two based on mutual enrichment

TF2GO<-function(genes_of_interest, all_genes, node_size, p_cutoff,
                TF_TERM2GENE){
  
  ###Perform GO Enrichment with topGO
  
  GOI<-unique(mapIds(org.Hs.eg.db, genes_of_interest, 'ENTREZID', 'SYMBOL'))
  background_names<-unique(mapIds(org.Hs.eg.db, all_genes, 'ENTREZID', 'SYMBOL'))
  geneList <- factor(as.integer(background_names %in% GOI))
  names(geneList) <- background_names
  GOdata <- new("topGOdata", allGenes = geneList,
                nodeSize = node_size, ontology = "BP", annot = annFUN.org, mapping = "org.Hs.eg.db")
  resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  upRes <- GenTable(GOdata, classicFisher = resultFisher, ranksOf = "classicFisher", numChar = 500)
  upRes$classicFisher[which(upRes$classicFisher=="<1e-30")]<-1*10^-30
  
  ###Create GO Term library 
  
  ensembl_symbol_GO = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
  connect_list<-list()
  for(t in 1:nrow(upRes)){
    genes<-biomaRt::getBM(attributes="hgnc_symbol", filters = "go", values = upRes$GO.ID[t],  mart = ensembl_symbol_GO)
    connect_list[[t]]<-genes[t(genes)%in%genes_of_interest,]
  }

  ###Connect Enriched Terms to Genes 
  
  GO_to_genes<-data.frame()
  
  for(g in 1:length(connect_list)){
    
    term_spec_genes<-as.vector(connect_list[[g]])
    connect_frame<-data.frame()
    
    if(length(term_spec_genes)>1){
      connect_frame<-cbind(rep(upRes$GO.ID[g], length(term_spec_genes)),term_spec_genes)
      GO_to_genes<-rbind(GO_to_genes, connect_frame)
    }
    
  }
     
  colnames(GO_to_genes)<-c("Term", "gene")
   
  ###Identify enriched TFs and connect to genes
  
  egmt <- enricher(genes_of_interest, TERM2GENE=rbind(TF_TERM2GENE), pAdjustMethod = "fdr", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
  egmt_sig<-egmt[which(egmt@result$p.adjust<p_cutoff),]  
  
  TF_to_genes<-data.frame()
  
  for(g in 1:nrow(egmt_sig)){
    
    TF<-egmt_sig$ID[g]
    genes<-unlist(stringr::str_split(egmt_sig$geneID[g], "/"))
    res<-cbind(TF, genes)
    TF_to_genes<-rbind(TF_to_genes, res)
  }
    
  ###Link GO terms to TF via common Genes

  unique_GO<-as.character(unique(GO_to_genes$Term))
  GO_to_TF_frame<-data.frame()
  
  for(x in 1:length(unique_GO)){
  
  genes_to_match_GO<-as.character(GO_to_genes$gene[GO_to_genes$Term%in%unique_GO[x]])
  TF_grow<-vector()
  
    for(y in 1:length(genes_to_match_GO)){
      
      query<-genes_to_match_GO[y]
      TF_grow<- c(TF_grow, unique(TF_to_genes$TF[which(TF_to_genes$genes%in%query)]))
      
    }
  
  GO_to_TF_frame<-rbind(GO_to_TF_frame, cbind(TF_grow, rep(as.character(unique_GO[x]), length(TF_grow))))
  
  }

  ###Convert GO IDs into names
  
  colnames(GO_to_TF_frame)<-c("TF", "GO_ID")
  
  GO_to_TF_frame$GO_ID<-upRes$Term[match(GO_to_TF_frame$GO_ID, upRes$GO.ID)]

  lc <- linkcomm::getLinkCommunities(GO_to_TF_frame, hcmethod = "average", plot = F, removetrivial = F, check.duplicates = T)
  lc2 <- linkcomm::newLinkCommsAt(lc, cutat =0.85)

  png("TF2GO Results.png", height = 2400, width = 2400, res = 300, pointsize = 12)
  plot(lc2, type = "graph", layout = igraph::layout.fruchterman.reingold, shownodesin = 1, scale.vertices = 0.1)#, node.pies = F)
  dev.off()

  plot(lc2, type = "graph", layout = "spencer.circle")
  
}

###Set your gene of interest list.  Must be a vector of universal gene symbols
genes_of_interest<-rownames(pc4.up.genes)

###Load a GMT file downloaded separately from mSigDB
c3.tft <- read.gmt("c3.tft.v7.1.symbols.gmt")
TF_TERM2GENE<-c3.tft

###Run the function
TF2GO(genes_of_interest = genes_of_interest,
      all_genes = all_genes<-rownames(SARP_BEC_data),
      node_size = 15,
      p_cutoff = 0.0001,
      TF_TERM2GENE = c3.tft)
      
      
