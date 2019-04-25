# this document contains functions created for the microbiomeSeq package https://github.com/umerijaz/microbiomeSeq
# make sure to cite the microbiomeSeq package if you use any of these functions

KDA <- function(physeq, grouping_column, analyse="abundance", method="dscore",
                
                p.adjust.method="BH",corr=0.99,adjusted.p.value.threshold= 0.05, select.variables=NULL, ...){
  
  
  
  
  meta_table <- data.frame(sample_data(physeq))
  
  
  
  
  # choose whether to analyse taxa differential abundance or differential variation
  
  # among measured environmental variables
  
  analyse <- match.arg(analyse,c("abundance","meta"))
  
  
  
  
  if(analyse=="abundance"){
    
    #enforce orientation of Otu_table , this is beacause dscore/sscore requiresamples to be rows.
    
    if(!taxa_are_rows(physeq)){
      
      
      
      
      abund_table <- t(abund_table)
      
    }
    
    analyse.data <- data.frame(abund_table)
    
  }
  
  else if(analyse=="meta"){
    
    
    
    
    analyse.data <- (get.num.variables(physeq))$num.variables
    
    
    
    
    if(!is.null(select.variables)){
      
      analyse.data <- analyse.data[,colnames(analyse.data)%in%select.variables]
      
    }
    
    
    
    
    analyse.data <- t(analyse.data) #transpose such that samples are rows
    
  }
  
  # pick out grouping variable
  
  meta_table$Groups <- meta_table[,grouping_column]
  
  # get levels of grouping variable
  
  group_levels<-levels(meta_table$Groups)
  
  # change the variables to 1 or 0 as required by the kmda procedure
  
  meta_table$Groups <- sapply(meta_table$Groups, function(x) ifelse(x == group_levels[1], 1, 0))
  
  # group taxa/environmental variables into sets using pearson correlation coefficients
  
  Sets<- KMDA::pearson.group(as.matrix(analyse.data),corr)
  
  Sets_mapping<-data.frame(row.names=rownames(analyse.data),Set=Sets)
  
  #obtain distance based kernel scores of each otu set
  
  #specify the lower and upper bounds of kernel parameter and number of gridpoints
  
  #selected in the interval
  
  KMtable<-NULL
  
  for(i in unique(Sets)){
    
    x=analyse.data[Sets==i,,drop=F]
    
    if(method=="dscore"){
      
      ks=KMDA:: dscore(x,meta_table$Groups,lower=1,upper=10,m=3)
      
    }
    
    if(method=="sscore"){
      
      ks=KMDA::sscore(x,meta_table$Groups,lower=10^-3,upper=10^3,m=10)
      
    }
    
    tmp<-data.frame(row.names=i,kscore=ks)
    
    ifelse(is.null(KMtable),KMtable<-tmp,KMtable<-rbind(KMtable,tmp))
    
  }
  
  #Adjust p-values
  
  if(is.null(p.adjust.method)){
    
    p.adjust.method <- "BH"
    
  }
  
  KMtable$kscore.adjusted<-p.adjust(KMtable$kscore,method=p.adjust.method,n=dim(KMtable)[1])
  
  
  
  
  #Select the sets of otus/env_variables that are differentially expressed
  
  kscore_selected<-rownames(KMtable)[KMtable$kscore.adjusted<=adjusted.p.value.threshold]
  
  #organise the data in a format accepted by ggplot2
  
  df<-NULL
  
  for(i in kscore_selected){
    
    tmp<- reshape2::melt(as.matrix(analyse.data[Sets==i,,drop=F]))
    
    colnames(tmp)<-c("Set_variable","Sample","Value")
    
    tmp$Groups<-meta_table[as.character(tmp$Sample),"Groups"]
    
    tmp$Set<-Sets_mapping[as.character(tmp$Set_variable),]
    
    tmp$Set<-paste("Set",tmp$Set,sprintf("Adj.p = %0.5g",KMtable[as.character(tmp$Set),"kscore.adjusted"]),
                   
                   cut(KMtable[as.character(tmp$Set),"kscore.adjusted"],breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))
    
    ifelse(is.null(df),df<-tmp,df<-rbind(df,tmp))
    
  }
  
  #change the grouping values back to original representation
  
  df$Groups<-sapply(df$Groups, function(x) ifelse(x == 1, group_levels[1], group_levels[2]))
  
  
  
  
  out <- list("plotdata"=df, "kscore_table"=KMtable)
  
  return(out)
  
}




plot_kda <- function(df){
  
  p<-ggplot2::ggplot(aes(x=Set_variable,y=Value,color=Groups),data=df)+theme_bw()+geom_boxplot(outlier.size = NA)+facet_grid(. ~ Set, drop=TRUE,scales="free",space="free_x")
  
  p<-p+ ggplot2::theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+theme(strip.text.x = element_text(size = 12, colour = "black", angle = 90,vjust = 0))
  
  p<-p+ ggplot2::theme(strip.background = element_rect(fill = "white"))+ylab("")+xlab("")
  
  return(p)
  
}


alpha_div <- function(physeq,method){
  
  #==check for validity of selected methods
  
  method<- match.arg(method,c("richness", "fisher", "simpson", "shannon", "evenness","pd"), several.ok = TRUE)
  
  
  
  
  abund_table <- otu_table(physeq)
  
  df <- NULL
  
  if("richness"%in%method){
    
    R<- vegan::rarefy(abund_table,min(rowSums(abund_table)))
    
    df_R<-data.frame(sample=names(R),value=R,measure=rep("Richness",length(R)))
    
    if(is.null(df)){
      
      df<-df_R}
    
    else {
      
      df<-rbind(df,df_R)}
    
  }
  
  if("fisher"%in%method){
    
    alpha <- vegan::fisher.alpha(abund_table)
    
    df_alpha<-data.frame(sample=names(alpha),value=alpha,measure=rep("Fisher alpha",length(alpha)))
    
    if(is.null(df)){
      
      df<-df_alpha}
    
    else {
      
      df<-rbind(df,df_alpha)}
    
  }
  
  if("simpson"%in%method){
    
    simp <- vegan::diversity(abund_table, "simpson")
    
    df_simp<-data.frame(sample=names(simp),value=simp,measure=rep("Simpson",length(simp)))
    
    if(is.null(df)){
      
      df<-df_simp}
    
    else {
      
      df<-rbind(df,df_simp)}
    
  }
  
  if("shannon"%in%method){
    
    H<- vegan::diversity(abund_table)
    
    df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon",length(H)))
    
    if(is.null(df)){
      
      df<-df_H}
    
    else {
      
      df<-rbind(df,df_H)}
    
  }
  
  if("evenness"%in%method){
    
    H<-vegan::diversity(abund_table)
    
    S <- specnumber(abund_table)
    
    J <- H/log(S)
    
    df_J<-data.frame(sample=names(J),value=J,measure=rep("Pielou's evenness",length(J)))
    
    if(is.null(df)){
      
      df<-df_J}
    
    else {
      
      df<-rbind(df,df_J)}
    
  }
  
  if("pd"%in%method){
    
    otu_tree <- phyloseq::phy_tree(physeq)
    
    PD <- pd(abund_table, otu_tree ,include.root = TRUE)
    
    df_PD<-data.frame(sample=names(PD),value=PD,measure=rep("PD",length(PD)))
    
    if(is.null(df)){
      
      df<-df_PD}
    
    else {
      
      df<-rbind(df,df_PD)}
    
  }
  
  return(df)
  
}

co_occurence_network <- function(physeq ,qval_threshold=0.05, grouping_column,rhos=c(-0.75,-0.5,0.5,0.75),select.condition=NULL,method="cor", filename=NULL, ...){
  
  
  
  
  corr_results <- network_correlation(physeq, grouping_column, select.condition, method, filename)
  
  
  
  
  comm_results<-comm_stats(corr_results$results, rhos)
  
  
  
  
  edge_lists<-co_occur_pairs(corr_results$results, rhos)
  
  
  
  
  edge_lists<-edge_lists[edge_lists$qval<=qval_threshold,]
  
  
  
  
  net <- construct_network(edge_lists, comm_results, ...)
  
  
  
  
  out <- list(net=net, comm_data=corr_results$comm.data, meta_table=corr_results$meta_table)
  
  
  
  
  return(out)
  
}




#computing correlation of taxa/features

network_correlation <- function(physeq, grouping_column, select.condition, method, filename){
  
  abund_table <- as.data.frame(otu_table(physeq))
  
  meta_table <- data.frame(sample_data(physeq))
  
  meta_table$Groups <- meta_table[,grouping_column]
  
  
  
  
  if(!is.null(select.condition)){meta_table <- subset(meta_table, Groups==select.condition)}
  
  
  
  
  # === rarefy
  
  comm.data<-vegan::rrarefy(abund_table,min(rowSums(abund_table)))
  
  trts<-unique(meta_table$Groups)
  
  
  
  
  results<-matrix(nrow=0,ncol=7)
  
  for(a in 1:length(trts)){
    
    trt.temp<-trts[a]
    
    #subset the dataset for those treatments
    
    temp<-subset(comm.data, meta_table$Groups==trt.temp)
    
    #remove empty species
    
    temp<-temp[,colSums(temp)>0]
    
    
    
    
    new.row<-NULL
    
    if(method=="bicor"){
      
      #Biweight midcorrelation
      
      res<- WGCNA::bicorAndPvalue(temp,use="p")
      
      
      
      
      res$bicor[upper.tri(res$bicor)] <- NA
      
      diag(res$bicor)<-NA
      
      res$p[upper.tri(res$p)]<-NA
      
      diag(res$p)<-NA
      
      
      
      
      new.row<-data.frame(rep(trts[a],dim(temp)[2]),reshape2::melt(as.matrix(res$bicor)),reshape2::melt(as.matrix(res$p))[,3])
      
      
      
      
    } else if (method=="cor"){
      
      #Pearson correlation
      
      res<- WGCNA::corAndPvalue(temp,use="p")
      
      
      
      
      res$bicor[upper.tri(res$cor)] <- NA
      
      diag(res$cor)<-NA
      
      res$p[upper.tri(res$p)]<-NA
      
      diag(res$p)<-NA
      
      
      
      
      new.row<-data.frame(rep(trts[a],dim(temp)[2]), reshape2::melt(as.matrix(res$cor)),reshape2::melt(as.matrix(res$p))[,3])
      
    }
    
    colnames(new.row)<-c("trt","taxa1","taxa2","rho","p.value")
    
    new.row<-data.frame(new.row,ab1=as.vector(sapply(as.character(new.row$taxa1),function(x){colSums(temp)[x]})),ab2=as.vector(sapply(as.character(new.row$taxa2),function(x){colSums(temp)[x]})))
    
    results<-rbind(results,new.row)
    
  }
  
  
  
  
  results$taxa1<-as.character(results$taxa1)
  
  results$taxa2<-as.character(results$taxa2)
  
  #We remove the upper triangle and diagonal from the correlation just calculated using NA assignments done before
  
  results<-results[complete.cases(results),]
  
  results[(results$ab1 <= 1) | (results$ab2 <= 1),"rho"]<-0
  
  results[(results$ab1 <= 1) | (results$ab2 <= 1),"p.value"]<-1
  
  
  
  
  if(!is.null(filename)){write.csv(results,paste(filename,"_co-occurence","_",paste(unique(meta_table$Groups),collapse="_vs_"),"_correlations.csv",sep=""))}
  
  
  
  
  out <- list(comm.data=comm.data, results=results, meta_table=meta_table)
  
  return(out)
  
}

co_occur_pairs<-function(dataset, rhos){
  
  
  
  final.results<-data.frame()
  
  
  
  
  trts<-as.vector(unique(dataset$trt))
  
  
  
  for(t in 1:length(trts)){
    
    
    
    dataset_trt<-subset(dataset, dataset$trt==trts[t])
    
    dataset_trt_no0<-subset(dataset_trt, dataset_trt$ab1 > 0 & dataset_trt$ab2 > 0)
    
    
    
    dataset_trt_no0$pairs<-paste(dataset_trt_no0$taxa1,dataset_trt_no0$taxa2)
    
    
    
    for(r in length(rhos)){
      
      
      
      if(rhos[r] < 0){temp<-subset(dataset_trt_no0, dataset_trt_no0$rho <= rhos[r])}
      
      if(rhos[r] > 0){temp<-subset(dataset_trt_no0, dataset_trt_no0$rho >= rhos[r])}
      
      if(dim(temp)[1]>1){
        
        
        
        temp.graph<-igraph::simplify(igraph::graph.edgelist(as.matrix(temp[,c(2,3)]),directed=FALSE))
        
        edge_list<-data.frame(igraph::get.edgelist(temp.graph,names=TRUE))
        
        
        
        edge_list$pairs<-paste(edge_list$X1,edge_list$X2)
        
        edge_list_pvals<-merge(edge_list,dataset_trt_no0,by="pairs",sort=FALSE  )
        
        
        
        edge_list_pvals$rho_cut<-rhos[r]
        
        edge_list_pvals$trt<-trts[t]
        
        
        
        edge_list_pvals$qval<-p.adjust(edge_list_pvals$p.value,method="fdr",n=length(edge_list_pvals$p.value))
        
        as.matrix(names(edge_list_pvals))
        
        
        
        final.results<-rbind(final.results,edge_list_pvals[,-c(2:3)])	}
      
    }
    
  }
  
  return(final.results)
  
}




comm_stats<-function(dataset, rhos){
  
  
  
  final.results<-data.frame()
  
  
  
  trts<-as.vector(unique(dataset$trt))
  
  
  
  for(t in 1:length(trts)){
    
    
    
    dataset_trt<-subset(dataset, trt==trts[t])
    
    
    
    for(r in 1:length(rhos)){
      
      
      
      if(rhos[r] < 0){temp<-subset(dataset_trt, rho <= rhos[r])}
      
      if(rhos[r] > 0){temp<-subset(dataset_trt, rho >= rhos[r])}
      
      if(dim(temp)[1]>1){
        
        
        
        temp.graph<-igraph::simplify(igraph::graph.edgelist(as.matrix(temp[,c(2,3)]),directed=FALSE))
        
        temp_comm<-igraph::edge.betweenness.community(temp.graph, directed=FALSE)
        
        
        
        member_data<-cbind(row.names(as.matrix(membership(temp_comm))),as.matrix(membership(temp_comm)))
        
        row.names(member_data)<-NULL
        
        
        
        
        
        rho_cut<-rep(rhos[r],dim(member_data)[1])
        
        trt<-rep(trts[t],dim(member_data)[1])
        
        stats<-cbind(trt,rho_cut,member_data)
        
        colnames(stats)<-c("trt","rho_cut","taxon","module")
        
        
        
        
        
        final.results<-rbind(final.results,stats)	}
      
    }
    
  }
  
  return(final.results)
  
}


construct_network <- function(edge_lists, comm_results, plotNetStats=F, plotNetwork=F,plotBetweennessEeigenvalue=F, scale.vertex.size=1, scale.edge.width=1, ...){
  
  
  
  colours <- microbiomeseq_cols()
  
  
  
  rho_cuts <- as.numeric(levels(comm_results$rho_cut)) 
  
  
  
  positive_rho_cuts <- rho_cuts[rho_cuts>0]
  
  
  
  g1 <- NULL
  
  
  
  for (i in unique(edge_lists$trt)){
    
    
    
    for (j in positive_rho_cuts){
      
      
      
      g1<-igraph::graph.edgelist(as.matrix(subset(edge_lists,trt==i)[,3:4]),directed=FALSE)
      
      
      
      E(g1)$color<-ifelse(subset(edge_lists,trt==i)[,"rho"]<0,"red","blue")
      
      
      
      E(g1)$width <- as.numeric(subset(edge_lists,trt==i)[,"rho"])*scale.edge.width
      
      
      
      V(g1)$size<-igraph::degree(g1)*scale.vertex.size
      
      
      
      #Give colours to subcommunities
      
      V(g1)$color<-as.character(sapply(V(g1)$name,function(x) ifelse(x %in% as.character(subset(comm_results,trt==i & rho_cut==j)[,"taxon"]),
                                                                     
                                                                     colours[as.numeric(as.character(subset(comm_results,trt==i & rho_cut==j & taxon==x)[,"module"]))],"white")))
      
      
      
      V(g1)$module<-as.character(sapply(V(g1)$name,function(x) if(x %in% as.character(subset(comm_results,trt==i & rho_cut==j)[,"taxon"])){
        
        return(subset(comm_results,trt==i & rho_cut==j & taxon==x)[,"module"])}))
      
      
      
      #plot network
      
      if(plotNetwork){
        
        pdf(paste("network_",paste(unique(edge_lists$trt),collapse="_vs_"),"_",i,"_",j,".pdf",sep=""),width=30,height=30)
        
        par(mai=c(0,0,1,0))       #this specifies the size of the margins. the default settings leave too much free space on all sides (if no axes are printed)
        
        plot(g1, layout=layout.fruchterman.reingold, main=i, vertex.label.dist=0.05,  vertex.frame.color='blue', 	
             
             vertex.label.color='black', vertex.label.font=2,vertex.label=V(g1)$name, ...)
        
        dev.off()
        
        
        
      }
      
      
      
      #plot betweenness - Eigen value plots
      
      if(plotBetweennessEeigenvalue){
        
        p <- plot_betweenness_eigenvalue(g1)
        
        pdf(paste("Betweeness_Eigenvalue_",paste(unique(edge_lists$trt),collapse="_vs_"),"_",i,"_",j,"_",".pdf",sep=""))
        
        print(p)
        
        dev.off()
        
      }
      
      
      
    }
    
  }
  
  out <- list(edgelists=edge_lists, comm_results=comm_results, graph=g1)
  
  return(out)
  
}







#produce a plot of betweenness Vs eigen value centrality.

plot_betweenness_eigenvalue <- function(g1){
  
  
  
  df<-data.frame(Betweenness=igraph::betweenness(g1),Eigenvalue=igraph::evcent(g1)$vector,Subcommunity=V(g1)$module, Taxa=V(g1)$name,Degree=igraph::degree(g1),Closeness=igraph::closeness(g1))
  
  df_community_mean<-aggregate(df[,c("Eigenvalue","Betweenness")],by=list(Subcommunity=df$Subcommunity),FUN="median")
  
  df_community_summary<-aggregate(df[,c("Eigenvalue","Betweenness")],by=list(Subcommunity=df$Subcommunity),FUN="summary")
  
  colnames(df_community_summary$Eigenvalue)<-paste("Eigenvalue",colnames(df_community_summary$Eigenvalue))
  
  colnames(df_community_summary$Betweenness)<-paste("Betweenness",colnames(df_community_summary$Betweenness))
  
  df_community_mean<-data.frame(df_community_mean,df_community_summary$Eigenvalue,df_community_summary$Betweenness)
  
  df_community_mean<-df_community_mean[df_community_mean$Subcommunity!="white",]
  
  colnames(df_community_mean)<-gsub("_$","",gsub("\\.","_",colnames(df_community_mean)))
  
  
  
  p<-ggplot(df)
  
  p<-p+geom_point(aes(x=Eigenvalue, y=Betweenness,size=Degree,fill=Subcommunity), colour="blue", pch=21,alpha=0.2)
  
  p<-p+scale_fill_manual(values=unique(get.vertex.attribute(g1, "color")))
  
  p<-p+theme_bw()
  
  p<-p+xlab("Eigenvalue Centrality")+ylab("Betweenness Centrality")
  
  p<-p+geom_point(aes(x=Eigenvalue,y=Betweenness,fill=Subcommunity),data=df_community_mean,colour="white",pch=23,size=7)
  
  
  
  
  if(dim(df_community_mean)[1]>0){
    
    p<-p+geom_errorbar(aes(x=Eigenvalue,y=Betweenness,ymin =Betweenness_1st_Qu,ymax = Betweenness_3rd_Qu,colour=Subcommunity),width=0.01,data=df_community_mean)
    
    p<-p+geom_errorbarh(aes(x=Eigenvalue,y=Betweenness,xmin =Eigenvalue_1st_Qu,xmax = Eigenvalue_3rd_Qu,colour=Subcommunity),height=0.01,data=df_community_mean)
    
    p<-p+scale_colour_manual(values=levels(df_community_mean$Subcommunity),guide = FALSE)
    
  }
  
  return(p)
  
}


differential_abundance <- function(physeq, grouping_column,pvalue.threshold=0.05,lfc.threshold=0,
                                   
                                   filename = "NB_significant",output_norm=NULL){
  
  
  
  
  abund_table <- otu_table(physeq)
  
  meta_table <-data.frame(sample_data(physeq))
  
  
  
  
  #==get count data to be used in deseq ==========#
  
  countData = round(as(abund_table, "matrix"), digits = 0)+1
  
  if(!taxa_are_rows(physeq)){
    
    countData = t(countData)
    
  }
  
  #==add a dummy column corresponding to grouping variable ==#
  
  meta_table$Groups <- meta_table[,grouping_column]
  
  #== create deseq compatible matrix and then test differential abundance in the groups=#
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData, meta_table, as.formula(~Groups))
  
  data_deseq_test = DESeq2::DESeq(dds)
  
  #==Extract results of the test =======#
  
  res = DESeq2::results(data_deseq_test, cooksCutoff = FALSE)
  
  res_tax = cbind(as.data.frame(res), as.matrix(countData[rownames(res), ]), OTU = rownames(res))
  
  res_tax_sig = subset(res_tax, padj < pvalue.threshold & lfc.threshold < abs(log2FoldChange))
  
  res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
  
  res_tax$Significant[is.na(res_tax$Significant)] <- "No"
  
  
  
  
  #==normalising abundance data for output ==#
  
  if(!is.null(output_norm)){
    
    physeq <- normalise_data(physeq,output_norm)
    
    data <- as.data.frame(otu_table(physeq))
    
  }
  
  else{
    
    data <- abund_table
    
  }
  
  
  
  
  #==get importance measure of significant features using random forest
  
  subset.data<-data.frame(data[,as.character(res_tax[rownames(res_tax_sig),"OTU"])])
  
  rownames(res_tax_sig) <- colnames(subset.data) #enforce rownames of res_tax_sig to be the same as colnames of subset data for easy indexing during plotting.
  
  
  
  
  rf_res <- randomforest_res(subset.data, meta_table$Groups)
  
  df_accuracy <- rf_res$importance
  
  
  
  
  #organise data for output
  
  df <- NULL
  
  for(i in df_accuracy$Sample){
    
    rank <- (subset(df_accuracy, df_accuracy$Sample==i))$rank
    
    tmp<-data.frame(subset.data[,i],meta_table$Groups, rep(rank), rep(paste(i," padj = ", sprintf("%.5g",res_tax_sig[i,"padj"]), sep=""), dim(data)[1]))
    
    colnames(tmp)<-c("Value","Groups","Rank","Taxa")
    
    if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)}
    
    df <- na.omit(df)
    
  }
  
  
  
  
  data_to_write=NULL
  
  if(!is.null(filename)){
    
    data_to_write<-res_tax_sig[,c("baseMean","log2FoldChange","pvalue","padj")]
    
    data_to_write$Upregulated<-levels(meta_table[,grouping_column])[as.numeric(data_to_write$log2FoldChange>0)+1]
    
    write.csv(data_to_write,paste(filename,paste(levels(meta_table$Groups),collapse="_vs_"),".csv",sep=""))
    
  }
  
  
  
  
  out_put <- list("SignFeaturesTable"=res_tax,"plotdata"=df,"importance"=df_accuracy)
  
  
  
  
  return(out_put)
  
  
  
  
}


taxa.env.correlation <- function(physeq, grouping_column, method="pearson", pvalue.threshold=0.05,
                                 
                                 padjust.method="BH", adjustment=1, num.taxa=50, select.variables=NULL){
  
  
  
  
  method<- match.arg(method,c("pearson", "kendall", "spearman"),several.ok = F)
  
  #enforce orientation
  
  if(taxa_are_rows(physeq)){
    
    physeq <- t(physeq)
    
  }
  
  abund_table <- as.data.frame(otu_table(physeq))
  
  #abund_table <- t(otu_table(physeq))
  
  meta_table <- data.frame(sample_data(physeq))
  
  #get grouping information
  
  groups<-meta_table[,grouping_column]
  
  #select variables to show (exclude) in (from) correlation plot.
  
  if(!is.null(select.variables)){
    
    meta_table <- subset(meta_table,select=select.variables)
    
  }
  
  #pick the numerical environmental variables since correlation function only accepts numerical variable
  
  mt_env <- meta_table[,sapply(meta_table,is.numeric)]
  
  
  
  
  #Now get a filtered abundance table based on selected variables
  
  abund_table_filt<-abund_table[rownames(mt_env),]
  
  
  
  
  # #pick top most  num.taxa taxa
  
  # select.top.taxa <- top.taxa(abund_table, num.taxa)
  
  # abund_table_filt <- select.top.taxa$abund_table
  
  abund_table_filt<-abund_table_filt[,order(colSums(abund_table_filt),decreasing=TRUE)]
  
  #Extract list of top N Taxa
  
  taxa_list<-colnames(abund_table_filt)[1:num.taxa]
  
  #remove "__Unknown__" and add it to others
  
  taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
  
  abund_table_filt<-data.frame(abund_table_filt[,colnames(abund_table_filt) %in% taxa_list])
  
  
  
  
  #Now calculate the correlation between individual Taxa and the environmental data
  
  df <- tables.correlate(abund_table_filt, mt_env, groups, method)
  
  colnames(df)<-c("Taxa","Env","Correlation","Pvalue","Type")
  
  df$Pvalue<-as.numeric(as.character(df$Pvalue))
  
  df$Correlation<-as.numeric(as.character(df$Correlation))
  
  # add column for adjusted p-values
  
  df$AdjPvalue<-rep(0,dim(df)[1])
  
  
  
  
  # correct pvalues for multiple testing
  
  df <- p.adjust.cor(df, adjustment, padjust.method)
  
  #Now we generate the labels for signifant values
  
  df$Significance<-cut(df$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  
  #We ignore NAs
  
  df<-df[complete.cases(df),]
  
  return(df)
  
}




tables.correlate<-function(table1, table2, groups=NULL, method){
  
  df<-NULL
  
  for(i in colnames(table1)){
    
    for(j in colnames(table2)){
      
      
      
      
      if(!is.null(groups)){
        
        for(k in unique(groups)){
          
          a<-table1[groups==k,i,drop=F]
          
          b<-table2[groups==k,j,drop=F]
          
          tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,k)
          
          
          
          
          if(is.null(df)){df<-tmp} else{df<-rbind(df,tmp)}
          
        }
        
      }
      
      else{
        
        
        
        
        a<-table1[,i,drop=F]
        
        b<-table2[,j,drop=F]
        
        tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value)
        
        
        
        
        if(is.null(df)){df<-tmp} else{df<-rbind(df,tmp)}
        
        
        
        
      }
      
      
      
      
    }
    
  }
  
  
  
  
  df<-data.frame(row.names=NULL,df)
  
  return(df)
  
}







plot_taxa_env <- function(df){
  
  p <-ggplot2::ggplot(aes(x=Type, y=Taxa, fill=Correlation), data=df)
  
  p <- p + ggplot2::geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C")
  
  p<-p+ ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  
  p<-p+ ggplot2::geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL)
  
  p<-p+ ggplot2::facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x")
  
  p<-p+ ggplot2::xlab("Groups")
  
  p<-p+ ggplot2::theme(strip.background = element_rect(fill = "white"))
  
  return(p)
  
}




#You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):

# 1 -> donot adjust

# 2 -> adjust Env + Type (column on the correlation plot)

# 3 -> adjust Taxa + Type (row on the correlation plot for each type)

# 4 -> adjust Taxa (row on the correlation plot)

# 5 -> adjust Env (panel on the correlation plot)

#adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")




# df is a data frame

p.adjust.cor <- function(df,adjustment=1,padjust.method="BH"){
  
  if(adjustment==1){
    
    df$AdjPvalue<-df$Pvalue
    
  } else if (adjustment==2){
    
    for(i in unique(df$Env)){
      
      for(j in unique(df$Type)){
        
        sel<-df$Env==i & df$Type==j
        
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
        
      }
      
    }
    
  } else if (adjustment==3){
    
    for(i in unique(df$Taxa)){
      
      for(j in unique(df$Type)){
        
        sel<-df$Taxa==i & df$Type==j
        
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
        
      }
      
    }
    
  } else if (adjustment==4){
    
    for(i in unique(df$Taxa)){
      
      sel<-df$Taxa==i
      
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
      
    }
    
  } else if (adjustment==5){
    
    for(i in unique(df$Env)){
      
      sel<-df$Env==i
      
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
      
    }
    
  }
  
  return(df)
  
}


formatPvalues <- function(pvalue) {
  
  ra<-""
  
  if(pvalue <= 0.1) ra<-"."
  
  if(pvalue <= 0.05) ra<-"*"
  
  if(pvalue <= 0.01) ra<-"**"
  
  if(pvalue <= 0.001) ra<-"***"
  
  return(ra)
  
}

#these are dependance functions for fso function as presented b

#Reference: http://www.nku.edu/~boycer/fso/

sim.binary <- function (df, method = NULL, diag=FALSE, upper=FALSE){
  
  df <- as.matrix(df)
  
  a <- df %*% t(df)
  
  b <- df %*% (1 - t(df))
  
  c <- (1 - df) %*% t(df)
  
  d <- ncol(df) - a - b - c
  
  #1=Jaccard, 2=Baroni-Urbani & Buser, 3=Kulczynski, 4=Ochiai, 5=Sorensen
  
  if (method == 1) {
    
    sim <- a/(a + b + c)
    
  }
  
  else if (method == 2) {
    
    sim <- (a + sqrt(a*d))/(a + b + c + sqrt(a*d))
    
  }
  
  else if (method == 3) {
    
    sim <- 0.5* (a/(a + b) + a/(a + c))
    
  }
  
  else if (method == 4) {
    
    sim <- a/sqrt((a + b) * (a + c))
    
  }
  
  else if (method == 5) {
    
    sim <- 2 * a/(2 * a + b + c)
    
  }
  
  sim2 <- sim[row(sim) > col(sim)]
  
  class(sim2) <- "dist"
  
  attr(sim2, "Labels") <- dimnames(df)[[1]]
  
  attr(sim2, "Diag") <- diag
  
  attr(sim2, "Upper") <- upper
  
  attr(sim2, "Size") <- nrow(df)
  
  attr(sim2, "call") <- match.call()
  
  
  
  return(sim2)
  
}




sim.abund <- function (df, method = NULL, diag = FALSE, upper = FALSE) 
  
{
  
  METHODS <- c("Baroni-Urbani & Buser", "Horn", "Yule (Modified)")
  
  if (!inherits(df, "data.frame")) 
    
    stop("df is not a data.frame")
  
  if (any(df < 0)) 
    
    stop("non negative value expected in df")
  
  if (is.null(method)) {
    
    cat("1 = Baroni-Urbani & Buser\n")
    
    cat("2 = Horn\n")
    
    cat("3 = Yule\n")
    
    cat("Select an integer (1-3): ")
    
    method <- as.integer(readLines(n = 1))
    
  }
  
  df <- as.matrix(df)
  
  sites <- nrow(df)
  
  species <- ncol(df)
  
  sim <- array(0, c(as.integer(sites),as.integer(sites)))
  
  spmax <- apply(df,2,max)
  
  
  
  if (method == 1) {
    
    #compute similarities (Baroni-Urbani & Buser)
    
    for (x in 1:sites) {
      
      for (y in 1:sites) {
        
        h1 <- 0
        
        h2 <- 0
        
        h3 <- 0
        
        for (i in 1:species) {
          
          h1 <- h1 + min(df[x,i],df[y,i])
          
          h2 <- h2 + max(df[x,i],df[y,i])
          
          h3 <- h3 + spmax[i] - max(df[x,i],df[y,i])
          
        }
        
        numer <- h1 + sqrt(h1*h3)
        
        denom <- h2 + sqrt(h1*h3)
        
        sim[x,y] <- ifelse(identical(denom,0), 0, numer/denom)
        
      }
      
    }
    
  }
  
  
  
  else if (method == 2) {
    
    #compute similarities (Horn)
    
    for (x in 1:sites) {
      
      for (y in 1:sites) {
        
        h1 <- 0
        
        h2 <- 0
        
        h3 <- 0
        
        for (i in 1:species) {
          
          if((df[x,i] + df[y,i]) > 0) h1 <- h1 + (df[x,i] + df[y,i]) * log10(df[x,i] + df[y,i])
          
          if(df[x,i] > 0) h2 <- h2 + df[x,i] * log10(df[x,i])
          
          if(df[y,i] > 0) h3 <- h3 + df[y,i] * log10(df[y,i])
          
        }
        
        x.sum <- sum(df[x,])
        
        y.sum <- sum(df[y,])
        
        xy.sum <- x.sum + y.sum
        
        if (identical(xy.sum, 0)) (sim[x,y] <- 0) else (sim[x,y] <- (h1 - h2 - h3)/(xy.sum * log10(xy.sum)-x.sum * log10(x.sum) - y.sum * log10(y.sum)))
        
        if (sim[x,y] < 1.0e-10) sim[x,y] <- 0
        
      }
      
    }
    
    
    
  }
  
  
  
  else if (method == 3) {
    
    #compute similarities (Yule)
    
    for (x in 1:sites) {
      
      for (y in 1:sites) {
        
        h1 <- 0
        
        h2 <- 0
        
        h3 <- 0
        
        h4 <- 0
        
        for (i in 1:species) {
          
          h1 <- h1 + min(df[x,i], df[y,i])
          
          h2 <- h2 + max(df[x,i], df[y,i]) - df[y,i]
          
          h3 <- h3 + max(df[x,i], df[y,i]) - df[x,i]
          
          h4 <- h4 + spmax[i] - max(df[x,i], df[y,i])
          
        }
        
        numer <- sqrt(h1*h4)
        
        denom <- sqrt(h1*h4) + sqrt(h2*h3)
        
        sim[x,y] <- ifelse(identical(denom,0), 0, numer/denom)
        
      }
      
    }
    
  }  
  
  
  
  else stop("Non convenient method")
  
  sim >- t(sim)
  
  sim2 <- sim[row(sim) > col(sim)]
  
  attr(sim2, "Size") <- sites
  
  attr(sim2, "Labels") <- dimnames(df)[[1]]
  
  attr(sim2, "Diag") <- diag
  
  attr(sim2, "Upper") <- upper
  
  attr(sim2, "method") <- METHODS[method]
  
  attr(sim2, "call") <- match.call()
  
  class(sim2) <- "dist"
  
  return(sim2)
  
}




st.acr <- function(sim, diag=FALSE, upper = FALSE)
  
{
  
  dis <- 1 - sim
  
  #use "shortest" or "extended"
  
  edis <- as.matrix(stepacross(dis, path = "shortest", toolong = 1))
  
  sim <- 1 - edis
  
  amax <- max(sim)
  
  amin <- min(sim)
  
  sim <- (sim-amin)/(amax-amin)
  
  sim2 <- sim[row(sim) > col(sim)]
  
  attr(sim2, "Size") <- nrow(sim)
  
  attr(sim2, "Labels") <- dimnames(sim)[[1]]
  
  attr(sim2, "Diag") <- diag
  
  attr(sim2, "Upper") <- upper
  
  attr(sim2, "call") <- match.call()
  
  class(sim2) <- "dist"
  
  return(sim2)
  
}

#/===dependencies===#


generateFSO<-function(physeq,grouping_column,method=1,indices=NULL,filename=NULL, type=1,step_across=F){
  
  data <- as.data.frame(otu_table(physeq))
  
  meta_table <- data.frame(sample_data(physeq))
  
  #pick only the numeric variables of metadata
  
  Env <- meta_table[,sapply(meta_table,is.numeric)]
  
  #create dummy column for grouping information
  
  meta_table$Groups <- meta_table[,grouping_column]
  
  #get indices of variables to investigate
  
  if(is.null(indices)){
    
    indices<-seq(1:dim(Env)[2])
    
  }
  
  sim <- sim.abund(data,method=method)
  
  dis.ho <- 1 - sim
  
  df<-NULL
  
  if(type==1){
    
    for(i in names(Env)[indices]){
      
      param.fso<-fso(Env[,i],dis.ho,permute=1000)
      
      tmp<-data.frame(mu=param.fso$mu,param=param.fso$data, Groups=as.factor(meta_table$Groups),
                      
                      label=rep(paste(i,"(",round(param.fso$r,2)," ",formatPvalues(param.fso$p),")",sep=""),dim(data)[1]))
      
      if(is.null(df)){df<-tmp} else {df<-rbind(df,tmp)}
      
    }
    
  }
  
  else if(type==2){
    
    if(step_across){
      
      sim <- st.acr(sim)
      
    }
    
    whole.mfso<-mfso(as.formula(paste("~",paste(lapply(indices,function(x) paste("Env[,",x,"]",sep="")),collapse="+"))),dis=dis.ho,data=Env,scaling=2,permute=1000)
    
    whole.mfso$var<-names(Env)[indices]
    
    names(whole.mfso$data)<-names(Env)[indices]
    
    #print(whole.mfso)
    
    for (i in 1:(ncol(whole.mfso$mu) - 1)) {
      
      for (j in (i + 1):ncol(whole.mfso$mu)) {
        
        cat(paste("Processing",names(whole.mfso$data)[i]," and ",names(whole.mfso$data)[j],"\n"))
        
        tmp<-data.frame(x=whole.mfso$mu[, i],y=whole.mfso$mu[, j], Groups=as.factor(meta_table$Groups),label=rep(paste("x=mu(",names(whole.mfso$data)[i],")",", y=mu(",names(whole.mfso$data)[j],")",sep=""),dim(data)[1]))
        
        if(is.null(df)){df<-tmp} else {df<-rbind(df,tmp)}
        
      }
      
    }
    
  }
  
  p <- ggplot(df, aes(param, mu))+geom_point(aes(colour = Groups)) +geom_smooth(method="lm", size=1, se=T) +theme_bw()+facet_wrap( ~ label , scales="free", ncol=3)
  
  p<-p+theme(strip.background = element_rect(fill = "white"))
  
  return(p)
  
}


kruskal_abundance <- function(physeq, grouping_column,pvalue.threshold=0.05)
  
{
  
  abund_table <- otu_table(physeq)
  
  meta_table <-data.frame(sample_data(physeq))
  
  meta_table$Groups <- meta_table[,grouping_column]
  
  
  
  
  kruskal.wallis.table <- data.frame()
  
  data <- as.data.frame(abund_table)
  
  for (i in 1:dim(data)[2]){
    
    ks.test <- kruskal.test(data[,i], g=meta_table$Groups)
    
    # Store the result in the data frame
    
    kruskal.wallis.table <- rbind(kruskal.wallis.table,data.frame(id=names(data)[i],p.value=ks.test$p.value))
    
  }
  
  kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]
  
  kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)
  
  kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value, decreasing=FALSE), ]
  
  kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
  
  kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor
  
  rownames(kruskal.wallis.table) <- kruskal.wallis.table$id
  
  
  
  
  #==================significant feature selection =================================#
  
  last.significant.element <- max(which(kruskal.wallis.table$q.value <= pvalue.threshold))
  
  selected <- 1:last.significant.element
  
  sig_res <-kruskal.wallis.table$id[selected]
  
  
  
  
  #==random forest classifier ==#
  
  subset.data<-data.frame(data[,as.character(kruskal.wallis.table[rownames(kruskal.wallis.table),"id"])])
  
  kruskal.wallis.table$id <- colnames(subset.data) #enforce that ids and colnames of subset data remain the same for easy indexing later on
  
  subset.data <- subset.data[,sig_res]
  
  rf_res <- randomforest_res(subset.data, meta_table$Groups)
  
  df_accuracy <- rf_res$importance
  
  
  
  
  df <- NULL
  
  for(i in df_accuracy$Sample){
    
    rank <- (subset(df_accuracy, df_accuracy$Sample==i))$rank
    
    tmp<-data.frame(subset.data[,i],meta_table$Groups, rep(rank), rep(paste(i," p.adj = ",sprintf("%.10g",kruskal.wallis.table[kruskal.wallis.table$id==i,"q.value"]),sep=""),dim(data)[1]))
    
    colnames(tmp)<-c("Value","Groups","Rank","Taxa")
    
    if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)}
    
    df <- na.omit(df)
    
  }
  
  
  
  
  out <- list("SignfeaturesTable"=kruskal.wallis.table, "plotdata"=df, "importance"=df_accuracy)
  
  
  
  
  return(out)
  
  
  
  
}



identify.num.variables <- function(physeq){
  
  meta_table <- sample_data(physeq)
  
  num.variables <- meta_table[,sapply(meta_table,is.numeric)]
  
  names.num.variables <- names(num.variables)
  
  return(names.num.variables)
  
}




get.num.variables <- function(physeq){
  
  meta_table <- sample_data(physeq)
  
  num.variables <- meta_table[,identify.num.variables(physeq)]
  
  not.num.variables <- meta_table[,!colnames(meta_table)%in%colnames(num.variables)]
  
  out <- list(num.variables=num.variables, notnum.variables=not.num.variables)
  
  return(out)
  
}




select.vars <- function(df1, df2, select.variables=NULL){
  
  if(!is.null(select.variables)){
    
    df1 <- df1[,colnames(df1)%in%select.variables]
    
  }
  
  out<-cbind(df1,df2)
  
  return(out)
  
}





module.roles <- function(comm_graph){
  
  
  
  
  td <- network_degree(comm_graph)
  
  
  
  
  wmd <- within_module_degree(comm_graph)
  
  
  
  
  z <- zscore(wmd)
  
  
  
  
  amd <- among_module_connectivity(comm_graph)
  
  
  
  
  pc <- participation_coeffiecient(amd, td)
  
  
  
  
  zp <- data.frame(z,pc)
  
  
  
  
  nod.roles <- assign_module_roles(zp)
  
  
  
  
  return(nod.roles)
  
  
  
  
}




# find total degree for each of the features in the graph




network_degree <- function(comm_graph){
  
  
  
  
  ki_total <-NULL
  
  
  
  
  net_degree <- degree(comm_graph)
  
  
  
  
  for(i in 1:length(V(comm_graph))){
    
    
    
    
    ki <- net_degree[i]
    
    
    
    
    tmp <- data.frame(taxa=names(ki), total_links=ki)
    
    
    
    
    if(is.null(ki_total)){ki_total<-tmp} else{ki_total <- rbind(ki_total, tmp)}
    
    
    
    
  }
  
  
  
  
  return(ki_total)
  
  
  
  
}




#compute within-module degree for each of the features




within_module_degree <- function(comm_graph){
  
  
  
  
  mods <- get.vertex.attribute(comm_graph, "module")
  
  
  
  
  vs <- as.list(V(comm_graph))
  
  
  
  
  modvs <- data.frame("taxon"= names(vs), "mod"=mods)
  
  
  
  
  sg1 <- decompose.graph(comm_graph,mode="strong")
  
  
  
  
  df <- data.frame()
  
  
  
  
  for(mod in unique(modvs$mod)){
    
    
    
    
    mod_nodes <- subset(modvs$taxon,modvs$mod==mod)
    
    
    
    
    neighverts <- unique(unlist(sapply(sg1,FUN=function(s){if(any(V(s)$name %in% mod_nodes)) V(s)$name else NULL})))
    
    
    
    
    g3 <- induced.subgraph(graph=comm_graph,vids=neighverts)
    
    
    
    
    mod_degree <- degree(g3)
    
    
    
    
    for(i in mod_nodes){
      
      
      
      
      ki <- mod_degree[which(names(mod_degree)==i)]
      
      
      
      
      tmp <- data.frame(module=mod, taxa=names(ki), mod_links=ki)
      
      
      
      
      df <- rbind(df,tmp)
      
      
      
      
    }
    
    
    
    
  }
  
  
  
  
  return(df)
  
  
  
  
}

#calculate the degree (links) of each node to nodes in other modules.




among_module_connectivity <- function(comm_graph){
  
  
  
  
  mods <- get.vertex.attribute(comm_graph, "module")
  
  
  
  
  vs <- as.list(V(comm_graph))
  
  
  
  
  modvs <- data.frame("taxa"= names(vs), "mod"=mods)
  
  
  
  
  df <- data.frame()
  
  
  
  
  for(i in modvs$taxa){
    
    
    
    
    for(j in modvs$taxa){
      
      
      
      
      if(are_adjacent(graph=comm_graph, v1=i , v2=j)){
        
        
        
        
        mod <- subset(modvs$mod, modvs$taxa==j)
        
        
        
        
        tmp <- data.frame(taxa=i, taxa2=j, deg=1, mod_links=mod)
        
        
        
        
        df <- rbind(df, tmp)
        
        
        
        
      }
      
      
      
      
    }
    
    
    
    
  }
  
  
  
  
  out <- aggregate(list(mod_links=df$deg), by=list(taxa=df$taxa, module=df$mod), FUN=sum)
  
  
  
  
  return(out)
  
  
  
  
}




#compute within-module degree z-score which

#measures how well-connected a node is to other nodes in the module.




zscore <- function(mod.degree){
  
  
  
  
  ksi_bar <- aggregate(mod_links ~ module, data=mod.degree, FUN = mean)
  
  
  
  
  ksi_sigma <- aggregate(mod_links ~ module, data=mod.degree, FUN = sd)
  
  
  
  
  z <- NULL
  
  
  
  
  for(i in 1:dim(mod.degree)[1]){
    
    
    
    
    mod_mean <- ksi_bar$mod_links[which(ksi_bar$module == mod.degree$module[i])]
    
    
    
    
    mod_sig <- ksi_sigma$mod_links[which(ksi_bar$module == mod.degree$module[i])]
    
    
    
    
    z[i] <- (mod.degree$mod_links[i] - mod_mean)/mod_sig
    
    
    
    
  }
  
  
  
  
  z <- data.frame(row.names=rownames(mod.degree), z, module=mod.degree$module)
  
  
  
  
  return(z)
  
  
  
  
}







#The participation coefficient of a node measures how well a  node is distributed

# in the entire network. It is close to 1 if its links are uniformly

#distributed among all the modules and 0 if all its links are within its own module.




participation_coeffiecient <- function(mod.degree, total.degree){
  
  
  
  
  p <- NULL
  
  
  
  
  for(i in total.degree$taxa){
    
    
    
    
    ki <- subset(total.degree$total_links, total.degree$taxa==i)
    
    
    
    
    taxa.mod.degree <- subset(mod.degree$mod_links, mod.degree$taxa==i)
    
    
    
    
    p[i] <- 1 - (sum((taxa.mod.degree)**2)/ki**2)
    
    
    
    
  }
  
  
  
  
  p <- as.data.frame(p)
  
  
  
  
  return(p)
  
  
  
  
}







assign_module_roles <- function(zp){
  
  
  
  
  zp <- na.omit(zp)
  
  
  
  
  zp$roles <- rep(0, dim(zp)[1])
  
  
  
  
  outdf <- NULL
  
  
  
  
  for(i in 1:dim(zp)[1]){
    
    
    
    
    df <- zp[i, ]
    
    
    
    
    if(df$z < 2.5){ #non hubs
      
      
      
      
      if(df$p < 0.05){
        
        
        
        
        df$roles <- "ultra peripheral"
        
        
        
        
      }
      
      else if(df$p < 0.620){
        
        
        
        
        df$roles <- "peripheral"
        
        
        
        
      }
      
      else if(df$p < 0.80){
        
        
        
        
        df$roles <- "non hub connector"
        
        
        
        
      }
      
      else{
        
        
        
        
        df$roles <- "non hub kinless"
        
        
        
        
      }
      
      
      
      
    }
    
    else { # module hubs
      
      
      
      
      if(df$p < 0.3){
        
        
        
        
        df$roles <- "provincial hub"
        
        
        
        
      }
      
      else if(df$p < 0.75){
        
        
        
        
        df$roles <- "connector hub"
        
        
        
        
      }
      
      else {
        
        
        
        
        df$roles <- "kinless hub"
        
        
        
        
      }
      
      
      
      
    }
    
    
    
    
    if(is.null(outdf)){outdf <- df}else{outdf <- rbind(outdf, df)}
    
    
    
    
  }
  
  
  
  
  return(outdf)
  
  
  
  
}




plot_roles <- function(node.roles, roles.colors=NULL){
  
  
  
  
  x1<- c(0, 0.05, 0.62, 0.8, 0, 0.30, 0.75)
  
  x2<- c(0.05, 0.62, 0.80, 1,  0.30, 0.75, 1)
  
  y1<- c(-Inf,-Inf, -Inf, -Inf,  2.5, 2.5, 2.5)
  
  y2 <- c(2.5,2.5, 2.5, 2.5, Inf, Inf, Inf)
  
  
  
  
  lab <- c("ultra peripheral","peripheral" ,"non-hub connector","non-hub kinless","provincial"," hub connector","hub kinless")
  
  
  
  
  if(is.null(roles.colors)){roles.colors <- c("#E6E6FA", "#DCDCDC", "#F5FFFA", "#FAEBD7", "#EEE8AA", "#E0FFFF", "#F5F5DC")}
  
  
  
  
  p <- ggplot() + geom_rect(data=NULL, mapping=aes(xmin=x1, xmax=x2, ymin=y1,ymax=y2, fill=lab))
  
  
  
  
  p <- p + guides(fill=guide_legend(title="Topological roles"))
  
  
  
  
  p  <- p + scale_fill_manual(values = roles.colors)
  
  
  
  
  p <- p + geom_point(data=node.roles, aes(x=p, y=z,color=module)) + theme_bw()
  
  
  
  
  p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("Participation Coefficient")+ylab(" Within-module connectivity z-score")
  
  
  
  
  return(p)
  
}



module_env_correlation <- function(co_occur_res, select.variables=NULL, ...){
  
  
  
  
  comdat <- co_occur_res$comm_data
  
  
  
  
  comm_cor <- co_occur_res$net$comm_results
  
  
  
  
  comm_graph <- co_occur_res$net$graph
  
  
  
  
  meta_table <- co_occur_res$meta_table[sapply(co_occur_res$meta_table, is.numeric)]
  
  
  
  
  net.comm <- data.frame(comdat[, colnames(comdat)%in%comm_cor$taxon])
  
  
  
  
  net.comm <- net.comm[rownames(net.comm)%in%rownames(meta_table),]
  
  
  
  
  if(!is.null(select.variables)){ meta_table<- meta_table[, select.variables]}
  
  
  
  
  comm.taxa.cor <- tables.correlate(net.comm, meta_table, method= "pearson")
  
  
  
  
  colnames(comm.taxa.cor) <- c("Taxa", "Env", "Correlation", "Pvalue")
  
  
  
  
  comm.taxa.cor$AdjPvalue<-rep(0,dim(comm.taxa.cor)[1])
  
  
  
  
  # correct pvalues for multiple testing
  
  comm.taxa.cor$Pvalue <- as.numeric(as.character(comm.taxa.cor$Pvalue))
  
  
  
  
  comm.taxa.cor <- p.adjust.cor(comm.taxa.cor)
  
  
  
  
  comm_cor$Taxa <- comm_cor$taxon
  
  
  
  
  comm.cor.merge <- merge.data.frame(comm_cor, comm.taxa.cor, by="Taxa")
  
  
  
  
  graph.btn <- betweenness(comm_graph)
  
  
  
  
  graph.btn <- data.frame(row.names = names(graph.btn), Taxa = names(graph.btn) ,betweenness=graph.btn)
  
  
  
  
  graph.btn <- graph.btn[rownames(graph.btn)%in%colnames(net.comm),]
  
  
  
  
  comm.cor.merge  <- merge.data.frame(comm.cor.merge, graph.btn, by="Taxa")
  
  
  
  
  df <- NULL
  
  
  
  
  for(i in levels(comm.cor.merge$module)){
    
    
    
    
    modi <- subset(comm.cor.merge, module==i)
    
    
    
    
    if(dim(modi)[1]>2){
      
      
      
      
      modi <- modi[which(modi$betweenness == max(modi$betweenness)),]
      
      tax <- modi$taxon
      
      pvalue <- as.numeric(as.character(modi$AdjPvalue))
      
      corr <- as.numeric(as.character(modi$Correlation))
      
      env <- modi$Env
      
      trt <- modi$trt
      
      
      
      
      tmp <- data.frame("Taxa"=paste("mod",i,"-", tax) , "Env"=env, "Type"=trt, "tax"=tax, "AdjPvalue"=pvalue, "Correlation"=corr)
      
      
      
      
      if(is.null(df)){df<-tmp} else {df <- rbind(df, tmp)}
      
      
      
      
    }
    
    
    
    
  }
  
  
  
  
  df$Significance<-cut(df$AdjPvalue,breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  
  
  
  
  df <- na.omit(df)
  
  
  
  
  return(df)
  
  
  
  
}

#1. edgernorm

edgeRnorm = function(physeq, ...){
  
  
  
  abund_table <- otu_table(physeq)
  
  
  
  # Enforce orientation.
  
  if(!taxa_are_rows(physeq) ){
    
    abund_table <- t(abund_table)
    
  }
  
  x = as(abund_table, "matrix")
  
  # See if adding a single observation, 1, 
  
  # everywhere (so not zeros) prevents errors
  
  # without needing to borrow and modify 
  
  # calcNormFactors (and its dependent functions)
  
  # It did. This fixed all problems. 
  
  # Can the 1 be reduced to something smaller and still work?
  
  x = x + 1
  
  # Now turn into a DGEList
  
  y = edgeR::DGEList(counts=x, remove.zeros=TRUE)
  
  # Perform edgeR-encoded normalization, using the specified method (...)
  
  z = edgeR::calcNormFactors(y, ...)
  
  # A check that we didn't divide by zero inside `calcNormFactors`
  
  if( !all(is.finite(z$samples$norm.factors)) ){
    
    stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
    
  }
  
  otu_table(physeq) <- otu_table(z$counts, taxa_are_rows=TRUE)
  
  return(physeq)
  
}

#2. variance stabilisation

deseq_varstab = function(physeq, sampleConditions=rep("A", nsamples(physeq)), ...){
  
  
  
  abund_table <- otu_table(physeq)
  
  
  
  # Enforce orientation.
  
  if(!taxa_are_rows(physeq) ){
    
    abund_table <- t(abund_table)
    
  }
  
  
  
  x = as(abund_table, "matrix")
  
  # The same tweak as for edgeR to avoid NaN problems
  
  # that cause the workflow to stall/crash.
  
  x = x + 1
  
  
  
  cds = DESeqDataSetFromMatrix(x, DataFrame(sampleConditions), ~ 1)
  
  # First estimate library size factors
  
  cds = estimateSizeFactors(cds)
  
  # Variance estimation, passing along additional options
  
  cds = estimateDispersions(cds, ...)
  
  #vsmat = varianceStabilizingTransformation(cds)
  
  vsmat <- varianceStabilizingTransformation(cds)
  
  otu_table(physeq) <- otu_table(assay(vsmat), taxa_are_rows=TRUE)
  
  return(physeq)
  
}

#3.proportion

# Normalize total sequences represented 




# Scale by dividing each variable by its standard deviation.

#physeq = transform_sample_counts(physeq, function(x) x/sd(x))

# Center by subtracting the median

#physeq = transform_sample_counts(physeq, function(x) (x-median(x)))

proportion = function(physeq){
  
  normf = function(x, tot=max(sample_sums(physeq))){ tot*x/sum(x) }
  
  physeq = transform_sample_counts(physeq, normf)
  
  return(physeq)
  
}




#4.random sampling 

randomsubsample = function(physeq, smalltrim=0.15, replace=TRUE,meta=F){
  
  # Set the minimum value as the smallest library quantile, n`smalltrim` 
  
  samplemin = sort(sample_sums(physeq))[-(1:floor(smalltrim*nsamples(physeq)))][1]
  
  physeqr = rarefy_even_depth(physeq, samplemin, rngseed=TRUE,replace=replace, trimOTUs=TRUE)
  
  return(physeqr)
  
}




#5.relative transformation

relative <- function(physeq,norm.meta=F,select.variables=NULL){
  
  if(norm.meta){
    
    get.vars <- get.num.variables(physeq)
    
    
    
    norm.variables <- get.vars$num.variables/rowSums(get.vars$num.variables)
    
    
    
    meta_table <- select.vars(norm.variables, get.vars$notnum.variables, select.variables)
    
    
    
    sample_data(physeq) <- meta_table
    
  }
  
  else if(!norm.meta){
    
    abund_table <- otu_table(physeq)
    
    otu_table(physeq) <- abund_table/rowSums(abund_table)
    
  }
  
  return(physeq)
  
}

#6. log relative transformation

log_relative <- function(physeq, norm.meta=F, select.variables=NULL){
  
  if(norm.meta){
    
    get.vars <- get.num.variables(physeq)
    
    norm.variables <- log(get.vars$num.variables/rowSums(get.vars$num.variables))
    
    meta_table <- select.vars(norm.variables, get.vars$notnum.variables, select.variables)
    
    sample_data(physeq) <- meta_table
    
  }
  
  else if(!norm.meta){
    
    abund_table <- otu_table(physeq)
    
    otu_table(physeq) <- log((abund_table+1)/(rowSums(abund_table)+dim(abund_table)[2]))
    
  }
  
  return(physeq)
  
}




#7.

scale.meta <- function(physeq,type="scale",select.variables=NULL){
  
  
  
  get.vars <- get.num.variables(physeq)
  
  num.variables <- get.vars$num.variables
  
  
  
  type <- match.arg(type, c("scale","log","sqrt"))
  
  
  
  if(type=="scale"){
    
    norm.variables <- data.frame(scale(num.variables))
    
  }
  
  else if(type=="log"){
    
    norm.variables <- log2(num.variables)
    
  }
  
  else if(type=="sqrt"){
    
    norm.variables <- sqrt(num.variables)
    
  }
  
  
  
  sample_data(physeq) <- select.vars(norm.variables, get.vars$notnum.variables, select.variables)
  
  
  
  return(physeq)
  
}




normalise_data <- function(physeq, norm.method, ...){
  
  norm.method = match.arg(norm.method,c("edgernorm","varstab","randomsubsample",
                                        
                                        "proportion","relative","log-relative","scale"))
  
  switch(norm.method,
         
         "randomsubsample"=randomsubsample(physeq),
         
         "proportion"=proportion(physeq),
         
         "varstab"=deseq_varstab(physeq, ...),
         
         "edgernorm"=edgeRnorm(physeq, ...),
         
         "log-relative"=log_relative(physeq, ...),
         
         "relative"=relative(physeq, ...),
         
         "scale"=scale.meta(physeq, ...)
         
  )
  
}


ordination <- function(physeq,which_distance="bray",method,grouping_column,pvalue.cutoff=0.05){
  
  
  
  
  meta_table <- data.frame(sample_data(physeq))
  
  meta_table$Groups <- meta_table[,grouping_column]
  
  
  
  
  sol<-NULL
  
  if(method=="PCoA"){
    
    sol<-cmdscale(phyloseq::distance(physeq,which_distance),eig=T)
    
  }
  
  else{
    
    sol<-phyloseq::ordinate(physeq,method,distance=which_distance)
    
  }
  
  
  
  
  dist<-phyloseq::distance(physeq,which_distance)
  
  adonis_res <- vegan::adonis(dist ~ Groups, data=meta_table)
  
  
  
  
  betadisper_pw <- beta_disper(physeq,grouping_column,pvalue.cutoff,which_distance)
  
  betadisper_pw <- betadisper_pw$betadisper_res
  
  
  
  
  out<- list("solution"=sol,"betadispersion"=betadisper_pw,"adonis_res"=adonis_res, "groups"=meta_table$Groups)
  
  return(out)
  
}




beta_disper <- function(physeq,grouping_column,pvalue.cutoff=0.05,which_distance){
  
  
  
  
  meta_table <- data.frame(sample_data(physeq))
  
  meta_table$Groups <- meta_table[,grouping_column]
  
  # compute beta dispersion
  
  mod<- vegan::betadisper(phyloseq::distance(physeq,method=which_distance),meta_table$Groups,type="centroid")
  
  # compute pairwise beta dispersion for all levels in the grouping variable
  
  pmod <- vegan::permutest(mod, permutations = 99, pairwise = TRUE)
  
  p.values <- pmod$pairwise$observed
  
  # extract significantly dispersed pairs and assign significant labels
  
  p.values <- p.values[!is.na(p.values <= pvalue.cutoff)]
  
  signi_label <- paste(cut(p.values,breaks=c(-Inf,0.001,0.01,0.05, Inf), label=c("***", "**", "*", ".")))
  
  groups_compared <- names(p.values)
  
  
  
  
  betadisper_res <- data.frame(groups_compared, p.values, signi_label)
  
  out <- list("betadisper_res"=betadisper_res, "pmod"=pmod)
  
  return(out)
  
}



perform_anova <- function(df,meta_table,grouping_column,pValueCutoff){
  
  
  
  
  dt<-data.table::data.table(data.frame(df,.group.=meta_table[,grouping_column]))
  
  #specifying a p-value cutoff for the ggplot2 strips
  
  pval<-dt[, list(pvalue = sprintf("%.2g",
                                   
                                   tryCatch(summary(aov(value ~ .group.))[[1]][["Pr(>F)"]][1],error=function(e) NULL))),
           
           by=list(measure)]
  
  #Filter out pvals that we are not significant
  
  pval<-pval[!pval$pvalue=="",]
  
  pval<-pval[as.numeric(pval$pvalue)<=pValueCutoff,]
  
  
  
  
  #using sapply to generate significances for pval$pvalue using the cut function.
  
  pval$pvalue<-sapply(as.numeric(pval$pvalue),function(x){as.character(cut(x,breaks=c(-Inf, 0.001, 0.01, pValueCutoff, Inf),label=c("***", "**", "*", "")))})
  
  
  
  
  #Update df$measure to change the measure names if the grouping_column has more than three classes
  
  if(length(unique(as.character(meta_table[,grouping_column])))>2){
    
    df$measure<-as.character(df$measure)
    
    if(dim(pval)[1]>0){
      
      for(i in seq(1:dim(pval)[1])){
        
        df[df$measure==as.character(pval[i,measure]),"measure"]=paste(as.character(pval[i,measure]),as.character(pval[i,pvalue]))
        
      }
      
    }
    
    df$measure<-as.factor(df$measure)
    
  }
  
  #Get all possible pairwise combination of values in the grouping_column
  
  s<-combn(unique(as.character(df[,grouping_column])),2)
  
  
  
  
  #df_pw will store the pair-wise p-values
  
  df_pw<-NULL
  
  for(k in unique(as.character(df$measure))){
    
    #We need to calculate the coordinate to draw pair-wise significance lines
    
    #for this we calculate bas as the maximum value
    
    bas<-max(df[(df$measure==k),"value"])
    
    #Calculate increments as 10% of the maximum values
    
    inc<-0.1*bas
    
    #Give an initial increment
    
    bas<-bas+inc
    
    for(l in 1:dim(s)[2]){
      
      #Do a pair-wise anova
      
      tmp<-c(k,s[1,l],s[2,l],bas,paste(sprintf("%.2g",tryCatch(summary(aov(as.formula(paste("value ~",grouping_column)),data=df[(df$measure==k) & (df[,grouping_column]==s[1,l] | df[,grouping_column]==s[2,l]),] ))[[1]][["Pr(>F)"]][1],error=function(e) NULL)),"",sep=""))
      
      #Ignore if anova fails
      
      if(!is.na(as.numeric(tmp[length(tmp)]))){
        
        #Only retain those pairs where the p-values are significant
        
        if(as.numeric(tmp[length(tmp)])<pValueCutoff){
          
          if(is.null(df_pw)){df_pw<-tmp}else{df_pw<-rbind(df_pw,tmp)}
          
          #Generate the next position
          
          bas<-bas+inc
          
        }
        
      }
      
    }
    
  }
  
  if(!is.null(df_pw)){
    
    df_pw<-data.frame(row.names=NULL,df_pw)
    
    names(df_pw)<-c("measure","from","to","y","p")
    
  }
  
  out <- list("df_pw"=df_pw, "df"=df)
  
  return(out)
  
}



plot.ordination <- function(ordination.res, method, pvalue.cutoff=0.05, show.pvalues=T, N=5,
                            
                            extra_marginspace=0.35){
  
  sol <- ordination.res$solution
  
  adn_res <- ordination.res$adonis_res
  
  betadisper_res <- ordination.res$betadispersion
  
  groups <- ordination.res$groups
  
  #This function is applied to each level of NMDS (group) and it uses also function
  
  #cov.wt to calculate covariance matrix.
  
  veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100){
    
    theta <- (0:npoints) * 2 * pi/npoints
    
    Circle <- cbind(cos(theta), sin(theta))
    
    t(center + scale * t(Circle %*% chol(cov)))
    
  }
  
  
  
  
  ord_res<-data.frame(x=sol$points[,1],y=sol$points[,2],Groups=groups)
  
  ord_res.mean=aggregate(ord_res[,1:2],list(group=ord_res$Groups),mean)
  
  
  
  
  plot.new()
  
  ord<-ordiellipse(sol, groups,display = "sites", kind ="se", conf = 0.95, label = T)
  
  dev.off()
  
  #Generate ellipse points
  
  df_ell <- data.frame()
  
  #include a condition to check positive definitenessas required by
  
  #cholskey decomposition to avoid errors in vegancovellipse
  
  for(g in levels(ord_res$Groups)){
    
    if(g!="" && (g %in% names(ord)) && all(eigen(ord[[g]]$cov)$values>0)){
      
      df_ell <- rbind(df_ell, cbind(as.data.frame(with(ord_res[ord_res$Groups==g,],
                                                       
                                                       veganCovEllipse(ord[[g]]$cov ,ord[[g]]$center,ord[[g]]$scale))),Groups=g))
      
    }
    
  }
  
  colnames(df_ell)<-c("x","y","Groups")
  
  
  
  
  #coloring function
  
  gg_color_hue<-function(n){
    
    hues=seq(15,375,length=n+1)
    
    hcl(h=hues,l=65,c=100)[1:n]
    
  }
  
  
  
  
  cols=gg_color_hue(length(unique(ord_res$Groups)))
  
  
  
  
  p<-ggplot2::ggplot(data=ord_res,aes(x,y,colour=Groups))
  
  p<-p + ggplot2::geom_point(alpha=0.5,size = 2)
  
  p<-p+ggplot2::theme_bw()
  
  p<-p+ ggplot2::annotate("text",x=ord_res.mean$x,y=ord_res.mean$y,label=ord_res.mean$group,size=6,colour=cols,family="Courier",fontface="bold",alpha=0.8,vjust=0.3)
  
  p<-p+ ggplot2::geom_path(data=df_ell, aes(x=x, y=y), size=1, linetype=1,alpha=0.3)
  
  
  
  
  #axis labels
  
  if(method=="NMDS"){
    
    #annotate plot with stress value from NMDS results
    
    stress.value <- sol$stress
    
    stress.label <- paste("STRESS=",round(stress.value,4))
    
    p <- p + ggplot2::annotation_custom(grob = textGrob(label = stress.label, hjust = 0, gp = gpar(cex = 1.5,fontsize=8)),
                                        
                                        ymin = max(ord_res$y), ymax = max(ord_res$y),
                                        
                                        xmin = extra_marginspace+max(ord_res$x),xmax = extra_marginspace+max(ord_res$x))
    
    p<-p+xlab("NMDS1")+ylab("NMDS2")
    
  }
  
  else if(method=="PCoA"){
    
    p<-p+xlab(paste("Dim1 (",sprintf("%.4g",sol$eig[1]),"%)",sep=""))+ylab(paste("Dim2 (",sprintf("%.4g",sol$eig[2]),"%)",sep=""))
    
  }
  
  
  
  
  #add the adonis results on the plot using custom annoatation
  
  #this only happens if adonis results turn out significant
  
  gt<-NULL #ggtable table to accomodate ordination plot and adonis results
  
  if(!is.null(adn_res)){
    
    adn_pvalue<-adn_res[[1]][["Pr(>F)"]][1]
    
    adn_rsquared<-round(adn_res[[1]][["R2"]][1],3)
    
    #use the bquote function to format adonis results to be annotated on the ordination plot.
    
    signi_label <- paste(cut(adn_pvalue,breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ".")))
    
    adn_res_format <- bquote(atop(atop("PERMANOVA",R^2==~.(adn_rsquared)), atop("p-value="~.(adn_pvalue)~.(signi_label), phantom())))
    
    if(adn_pvalue<=pvalue.cutoff){
      
      gt <- plot_adonis_res(p,ord_res, adn_res_format,extra_marginspace)
      
    }
    
  }
  
  
  
  
  # add a table of beta dispersion results
  
  anova_table<-NULL
  
  if(!is.null(betadisper_res)){
    
    anova_table <- plot_betadisper(betadisper_res, show.pvalues,pvalue.cutoff, N)
    
    th <- sum(anova_table$heights)
    
  }
  
  
  
  
  out<-p
  
  if(!is.null(gt) && !is.null(anova_table)){
    
    out <-gridExtra::grid.arrange(gt,anova_table,heights = unit.c(unit(1, "null"), th))
    
  }
  
  else if(is.null(gt) && !is.null(anova_table)){
    
    out <- gridExtra::grid.arrange(p,anova_table,heights = unit.c(unit(1, "null"), th))
    
  }
  
  return(out)
  
}




# ==========function to annotate adonis (PERMANOVA) results onto an ordination plot

plot_adonis_res <- function(p, ord_res, adn_res,extra_marginspace){
  
  p<-p+theme(legend.position = "none") #get rid of the legend
  
  p<-p+theme(plot.margin = unit(c(1,8,1,1), "lines")) #allow extra space in margins
  
  #annotate plot with adonis results
  
  p <- p + annotation_custom(grob = textGrob(label = adn_res, hjust = 0, gp = gpar(cex = 1.5,fontsize=12)),
                             
                             ymin = median(ord_res$y), ymax = median(ord_res$y), xmin =extra_marginspace+ max(ord_res$x), xmax = extra_marginspace+max(ord_res$x))
  
  gt <- ggplot_gtable(ggplot_build(p))
  
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  
  return(gt)
  
}




#=====================A function to plot betadispersion results onto an ordination plot

plot_betadisper <- function(betadisper_res, show.pvalues=T, pvalue.cutoff, N){
  
  #order the significantly dispersed groups in incresing size of p-value
  
  if(dim(betadisper_res)[1]>=2){
    
    betadisper_res<-betadisper_res[order(betadisper_res$p.value),]
    
  }
  
  
  
  
  colnames(betadisper_res)<-c("groups","p_value","label")
  
  #display significantly dispersed groups,select the number of groups to display
  
  #and whether or not to show p-values select groups to display on plot
  
  anova_table_display <- subset(betadisper_res, p_value<=pvalue.cutoff)
  
  print(anova_table_display)
  
  if(dim(anova_table_display)[1]<N){
    
    N <- dim(anova_table_display)[1]
    
  }
  
  if(show.pvalues){
    
    anova_table_display$p_label <- paste(rep("p-value="),sprintf("%.e",anova_table_display$p_value)," ",anova_table_display$label, sep="")
    
    anova_table_display<-anova_table_display[1:N, c("groups","p_label")]
    
  }
  
  else{
    
    anova_table_display<-anova_table_display[1:N,c("groups","label")]
    
  }
  
  #theme_minimal in tableGrob provides a white background
  
  anova_table <- gridExtra::tableGrob(anova_table_display,rows = NULL,cols = NULL,theme =ttheme_minimal()) #set rows/cols to null to get rid of the row numbering
  
  title <- grid::textGrob("BETA-DISPERSION",gp=gpar(fontsize=12))
  
  anova_table <- gtable::gtable_add_rows( anova_table, heights = grobHeight(title) + unit(5,"mm"), pos = 0)
  
  anova_table <- gtable::gtable_add_grob( anova_table,  title, 1, 1, 1, ncol(anova_table))
  
  ##calculate height of table and pass it to glob to avoid white space between plot and table
  
  return(anova_table)
  
}




plot_ordisurf <- function(sol, meta_table, env.variable, grouping_column){
  
  
  
  
  groups <- meta_table[,grouping_column] #get grouping information from meta data
  
  
  
  
  df=data.frame(x=sol$point[,1],y=sol$point[,2],Groups=groups)
  
  #Add a dummy variable corrresponding to the selected variable
  
  meta_table$var <- meta_table[,env.variable]
  
  
  
  
  #fit a surface for a selected variable onto ordination stats
  
  ordi<- vegan::ordisurf(sol,meta_table$var ,plot = FALSE, bs="ds")
  
  ordi.grid <- ordi$grid #extracts the ordisurf object
  
  #str(ordi.grid) #it's a list though - cannot be plotted as is
  
  ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
  
  ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
  
  ordi.mite.na <- data.frame(na.omit(ordi.mite)) #gets rid of the nas
  
  
  
  
  #make the plot
  
  p<-ggplot2::ggplot()+stat_contour(data = ordi.mite.na, aes(x = x, y = y, z = z, colour = ..level..),positon="identity") #can change the binwidth depending on how many contours you want
  
  p<-p+ ggplot2::geom_point(data=df,aes(x,y,fill=Groups),pch=21,size=3)
  
  p<-p+ ggplot2::scale_colour_continuous(high = "darkgreen", low = "darkolivegreen1") #here we set the high and low of the colour scale.  Can delete to go back to the standard blue, or specify others
  
  p<-p+ ggplot2::labs(colour = paste(env.variable)) #another way to set the labels, in this case, for the colour legend
  
  p<-p+ ggplot2::theme_bw()
  
  return(p)
  
}




plot_anova_diversity <- function(physeq, method, grouping_column,pValueCutoff=0.05, outfile="anova_diversity.csv")
  
{
  
  #enforce orientation
  
  if(taxa_are_rows(physeq)){
    
    physeq <- t(physeq)
    
  }
  
  abund_table <- otu_table(physeq)
  
  meta_table <- sample_data(physeq)
  
  
  
  
  #get diversity measure using selected methods
  
  div.df <- alpha_div(physeq,method)
  
  
  
  
  #=add grouping information to alpha diversity measures
  
  df<-data.frame(div.df,(meta_table[,grouping_column])[as.character(div.df$sample),])
  
  
  
  
  #perform anova of diversity measure between groups
  
  anova_res <- perform_anova(df,meta_table,grouping_column,pValueCutoff)
  
  df_pw <- anova_res$df_pw #get pairwise p-values
  
  write.csv(df_pw, file = outfile)
  
  #Draw the boxplots
  
  p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)
  
  p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
  
  p<-p+theme_bw()
  
  p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  p<-p+facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Observed Values")+xlab("Samples")
  
  p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("Groups")
  
  
  
  
  #This loop will generate the lines and signficances
  
  if(!is.null(df_pw)){ #this only happens when we have significant pairwise anova results
    
    for(i in 1:dim(df_pw)[1]){
      
      p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))), measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
      
      p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
      
    }
    
  }
  
  return(p)
  
  
  plot_anova_env <- function(physeq, grouping_column, pValueCutoff=0.05,select.variables=NULL){
    
    
    
    
    #get meta data from phyloseq object
    
    meta_table <- as.data.frame(sample_data(physeq))
    
    #pick numerical variables of environmental data
    
    env_table <- meta_table[,sapply(meta_table,is.numeric)]
    
    df<- reshape2::melt(as.matrix(env_table))
    
    names(df)<-c("sample","measure","value")
    
    #Incorporate categorical data in df
    
    df<-data.frame(df,(meta_table[, grouping_column])[as.character(df$sample),])
    
    
    
    
    #do anova of environmental variables between groups
    
    anova_res <- perform_anova(df,meta_table,grouping_column,pValueCutoff)
    
    df_pw <- anova_res$df_pw #get pairwise p-values
    
    df <- anova_res$df #get updated environmental measure information
    
    
    
    
    #pick selected variables
    
    if(!is.null(select.variables)){
      
      df <- df[which(df$measure%in%select.variables),]
      
      df_pw<-df_pw[which(df_pw$measure%in%select.variables),]
      
    }
    
    #Draw the boxplots
    
    p<-ggplot2::ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)
    
    p<-p+ggplot2::geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
    
    p<-p+ggplot2::theme_bw()+geom_point(size=1,alpha=0.2)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    p<-p+ggplot2::facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Observed Values")+xlab("Groups")
    
    p<-p+ggplot2::theme(strip.background = element_rect(fill = "white"))
    
    #This loop will generate the lines and signficances
    
    for(i in 1:dim(df_pw)[1]){
      
      p<-p+ ggplot2::geom_path(inherit.aes=F,aes(x,y),data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))), measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
      
      p<-p+ ggplot2::geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
      
    }
    
    return(p)
    
  }
  
  
  plot_cca <- function(physeq, grouping_column, pvalueCutoff=0.01,norm_method=NULL, env.variables=NULL,
                       
                       num.env.variables=NULL, exclude.variables=NULL, draw_species=F){
    
    #extract abundance and meta data from supplied phyloseq object
    
    abund_table <- otu_table(physeq)
    
    meta_table <- data.frame(sample_data(physeq))
    
    
    
    
    #Use adonis to find significant environmental variables
    
    abund_table.adonis <- vegan::adonis(abund_table ~ ., data=meta_table)
    
    #pick significant features
    
    bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<=pvalueCutoff]
    
    #throw out na in case any exist
    
    bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
    
    #pick out the selected variables if at all they are part of the significant ones
    
    if(!is.null(env.variables)&&(env.variables%in%bestEnvVariables)){
      
      bestEnvVariables <- env.variables
      
    }
    
    #provide number of variables to display
    
    if(!is.null(num.env.variables)){
      
      if(num.env.variables>length(bestEnvVariables)){
        
        stop(cat(paste("Choose a number less than",length(bestEnvVariables))))
        
      }else{
        
        bestEnvVariables <- bestEnvVariables[1:num.env.variables]
        
      }
      
    }
    
    #exclude selected variables from appearing on the plot
    
    if(!is.null(exclude.variables)&&(exclude.variables%in%bestEnvVariables)){
      
      bestEnvVariables <- bestEnvVariables[!(bestEnvVariables%in%exclude.variables)]
      
    }
    
    #We are now going to use only those environmental variables in cca that were found significant
    
    eval(parse(text=paste("sol <- cca(abund_table ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=meta_table)",sep="")))
    
    
    
    
    scrs<- vegan::scores(sol,display=c("sp","wa","lc","bp","cn"))
    
    #Extract site data first
    
    df_sites<-data.frame(scrs$sites,meta_table[,grouping_column])
    
    colnames(df_sites)<-c("x","y","Groups")
    
    
    
    
    #Draw sites
    
    p<-ggplot2::ggplot()
    
    p<-p+ggplot2::geom_point(data=df_sites,aes(x,y,colour=Groups))
    
    #Draw biplots
    
    multiplier <- vegan:::ordiArrowMul(scrs$biplot)
    
    df_arrows<- scrs$biplot*multiplier
    
    colnames(df_arrows)<-c("x","y")
    
    df_arrows=as.data.frame(df_arrows)
    
    p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.5)
    
    p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.5)
    
    # Draw species
    
    df_species<- as.data.frame(scrs$species)
    
    colnames(df_species)<-c("x","y")
    
    # Either choose text or points
    
    #p<-p+geom_text(data=df_species,aes(x,y,label=rownames(df_species)))
    
    if(draw_species){
      
      p<-p+geom_point(data=df_species,aes(x,y,shape="Species"))+scale_shape_manual("",values=2)
      
    }
    
    p<-p+theme_bw()+xlab("CCA1")+ylab("CCA2")
    
    return(p)
    
  }
  
  
  
  
  plot_signif <- function(df=NULL,top.taxa=20,...){
    
    
    
    
    #==plot the significant features information and random classifier results
    
    p<-NULL
    
    df <- na.omit(df)
    
    #pick the top significant taxa, this is to avoid overwhelming clutter.
    
    df <- df[which(as.numeric(df$Rank)%in%(1:top.taxa)),]
    
    max.rank <- max(as.numeric(df$Rank))
    
    df$bar_height <- max.rank-as.numeric(df$Rank)+1
    
    if(max.rank >=20){
      
      df$bar_height <- 1 + max.rank-round(as.numeric(df$Rank)/(top.taxa/10))
      
    }
    
    df$rank_label <-NULL
    
    for (i in 1:dim(df)[1]){
      
      rank <- as.numeric(df$Rank[i])
      
      df$rank_label[i] <- paste(paste(rep("\u25ac", df$bar_height[i]),collapse=""),rank) #"\u25aa" "\u25ae" "\u25ac"
      
    }
    
    
    
    
    if(!is.null(df)){
      
      p<-ggplot(df,aes(Groups,Value,colour=Groups))
      
      p<-p+geom_boxplot(outlier.size=NA)+geom_jitter(position = position_jitter(height = 0, width=0))+theme_bw()
      
      p<-p+ facet_grid( ~rank_label+Taxa, scales="free_x")
      
      p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+theme(strip.text = element_text(size = 10, colour = "black", angle = 90,vjust=0)) #vjust=0 aligns the labels to the bottom in this case
      
      p<-p+theme(strip.background =element_rect(fill="white"))+theme(plot.margin = unit(c(1, 1, 0, 1), "lines"))#set background colour for facet lbels
      
    }
    
    
    
    
    return(p)
    
  }
  
  
  
  
  #==generate stand alone plot for random forest results =============#
  
  plot_MDA <- function(df_accuracy, top.taxa=20){
    
    mda_plot <-NULL
    
    if(!is.null(df_accuracy)){
      
      df_accuracy <- df_accuracy[which(df_accuracy$rank%in%c(1:top.taxa)),]
      
      mda_plot <- ggplot(data = df_accuracy,aes(x=Sample,y=Value)) + theme_bw()
      
      mda_plot <- mda_plot+geom_bar(stat = "identity",fill="darkblue",width = 0.5)
      
      mda_plot <- mda_plot + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
      
      mda_plot <- mda_plot + xlab("Taxa description") + ylab("Mean Decrease in Accuracy")
      
    }
    
    return(mda_plot)
    
  }
  
  
  
  
  plot_MA <- function(res_tax=NULL, label=F){
    
    p1 <- NULL
    
    p2<-NULL
    
    if(!is.null(res_tax)){
      
      tax.display = NULL
      
      tax.aggregate = "OTU"
      
      p1 <- ggplot(data = res_tax, aes(x = baseMean, y = log2FoldChange, color = Significant)) + geom_point(size = 2) + scale_x_log10()
      
      p1 <- p1+scale_color_manual(values=c("black", "red"))+labs(x="Mean abundance",y="Log2 fold change")+theme_bw()
      
      if(label == T){
        
        if (!is.null(tax.display)){
          
          rlab <- data.frame(res_tax, Display = apply(res_tax[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
          
        }
        
        else {
          
          rlab <- data.frame(res_tax, Display = res_tax[,tax.aggregate])
          
        }
        
        p1 <- p1 + geom_text(data = subset(rlab, Significant == "Yes"), aes(label = Display), size = 4, vjust = 1)
        
      }
      
      
      
      
      #== produce a bar version of ma plot =========================#
      
      res_tax$cols<- ifelse(res_tax$log2FoldChange>=0, "Up regulated", "Down regulated")
      
      p2<-ggplot(data=res_tax,aes(x=rownames(res_tax),y=log2FoldChange,fill=cols))+geom_bar(stat="identity",width=0.5)+theme_bw()
      
      p2<-p2 + theme(axis.text.x = element_text(angle = 90,hjust = 1))+ xlab("Taxa Description") +ylab("Log2 Fold Change")+
        
        geom_text(aes(label=round(as.numeric(baseMean),1)), vjust=0, angle=90)
      
      p2<-p2+scale_fill_manual(values = c("Up regulated" = "darkblue", "Down regulated" = "red"))+theme(legend.title=element_blank())
      
    }
    
    
    
    
    out <- list("maplot"=p1,"lfcplot"=p2)
    
    return(out)
    
  }
  
  
  
  
  #== plot of multiple testing corrections ====
  
  plot_corrections <- function(corrections_table, pvalue.Cutoff){
    
    
    
    
    if(!is.null(corrections_table)){
      
      plot(corrections_table$p.value, corrections_table$E.value,main='Multitesting corrections',
           
           xlab='Nominal p-value',ylab='Multitesting-corrected statistics',log='xy',col='blue',panel.first=grid(col='#BBBBBB',lty='solid'))
      
      lines(corrections_table$p.value,corrections_table$FWER,pch=20,col='darkgreen', type='p')
      
      lines(corrections_table$p.value,corrections_table$q.value,pch='+',col='darkred', type='p')
      
      abline(h=pvalue.Cutoff, col='red', lwd=2)
      
      legend('topleft', legend=c('E-value', 'p-value', 'q-value'), col=c('blue', 'darkgreen','darkred'), lwd=2,bg='white',bty='o')
      
    }
    
    
    
    
  }
  
  
  plot_taxa <- function(physeq,grouping_column,method="hellinger",number.taxa=21,filename=NULL){
    
    
    
    
    #==extract components of the phyloseq object
    
    abund_table <- otu_table(physeq)
    
    meta_table <- data.frame(sample_data(physeq))
    
    #Enforce orientation of the phyloseq object
    
    if(taxa_are_rows(physeq) ){
      
      abund_table <- t(abund_table)
      
    }
    
    
    
    
    #===Calculate beta diversity and extract measure for local contribution to beta diversity
    
    beta_div<-beta.div(abund_table,method=method,sqrt.D=F,samp=T,nperm=999)
    
    df_LCBD<-data.frame(Sample=names(beta_div$LCBD),LCBD=beta_div$LCBD,p.LCBD=beta_div$p.LCBD)
    
    
    
    
    #=== add grouping information to the LCBD results
    
    df_LCBD<-data.frame(df_LCBD,Groups=meta_table[rownames(df_LCBD),grouping_column])
    
    
    
    
    if(!is.null(filename)){
      
      write.csv(df_LCBD,paste(filename,"_LCBD",".csv",sep=""))
      
    }
    
    
    
    
    select.top.taxa <- top.taxa(abund_table, number.taxa)
    
    new_x <- select.top.taxa$abund_table
    
    number.taxa <- select.top.taxa$number.taxa
    
    
    
    
    #arrange data for plotting in a format compatible to ggplot
    
    df<-NULL
    
    for (i in 1:dim(new_x)[2]){
      
      tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Groups=meta_table[,grouping_column])
      
      if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
      
    }
    
    df<-data.frame(df,df_LCBD[as.character(df$Sample),c("LCBD","p.LCBD","Groups")])
    
    #==plot the data
    
    colours <- microbiomeseq_cols()
    
    p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Groups, drop=TRUE,scale="free",space="free_x")
    
    p<- p+ guides(fill=guide_legend(ncol=1))+scale_fill_manual(values=colours[1:(number.taxa+1)])+theme_bw()+xlab("Samples")
    
    p<-p+ scale_y_continuous(expand = c(0.02,0))+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    
    p<-p+geom_point(aes(Sample,-0.02,size=LCBD))+theme(strip.background = element_rect(fill = "white"))
    
    return(p)
    
  }
  
  
  
  
  top.taxa <- function(abund_table, number.taxa){
    
    #==== sort the abundance table by total abundance of each taxa  in decreasing order
    
    abund_table<-abund_table[,order(colSums(abund_table),decreasing=TRUE)]
    
    #Extract list of top number.taxa Taxa
    
    taxa_list<-colnames(abund_table)[1:number.taxa]
    
    #remove "__Unknown__" and add it to others
    
    taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
    
    number.taxa<-length(taxa_list)
    
    #Generate a new table with everything added to Others
    
    new_x<-data.frame(abund_table[,colnames(abund_table) %in% taxa_list],Others=rowSums(abund_table[,!colnames(abund_table) %in% taxa_list]))
    
    
    
    
    out <- list("abund_table"=new_x, "number.taxa"=number.taxa)
    
    return(out)
    
  }
  
  
  
  
  microbiomeseq_cols <- function(){
    
    colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00",
                 
                 "#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000",
                 
                 "#FFFF00",grey.colors(1000));
    
    return(colours)
    
  }
  
  
  randomforest_res <- function(data, groups){
    
    IDs_map<-data.frame(row.names=colnames(data),"taxa"=colnames(data))
    
    val<-randomForest::randomForest(groups~ ., data=data, importance=T, proximity=T,ntree=1500,keep.forest=F)
    
    imp<- randomForest::importance(val)
    
    df_accuracy<-data.frame(row.names=NULL,Sample=rownames(imp),Value=abs(as.numeric(imp[,"MeanDecreaseAccuracy"])),Index=rep("Mean Decrease Accuracy",dim(imp)[1]))
    
    
    
    
    #Rearrange the features in terms of importance for ggplot2 by changing factor levels
    
    df_accuracy$Sample<-IDs_map[as.character(df_accuracy$Sample),"taxa"]
    
    df_accuracy_order<-as.character(IDs_map[rownames(imp),"taxa"][order(abs(as.numeric(imp[,"MeanDecreaseAccuracy"])),decreasing=T)])
    
    df_accuracy$Sample<-factor(as.character(df_accuracy$Sample),levels=df_accuracy_order)
    
    df_accuracy$rank <- base::rank(df_accuracy$Value, ties.method = "min")
    
    df_accuracy$rank <- max(df_accuracy$rank)-df_accuracy$rank+1
    
    
    
    
    out<-list("importance"=df_accuracy)
    
    return(out)
    
  }
  
  
  taxa_level <- function(physeq,which_level){
    
    #enforce orientation
    
    if(taxa_are_rows(physeq)){
      
      physeq <- t(physeq)
      
    }
    
    OTU <- otu_table(physeq)
    
    SAM <- sample_data(physeq)
    
    OTU_taxonomy <- tax_table(physeq)
    
    new_abund_table<-NULL
    
    if(which_level=="Otus"){
      
      OTU_tree <- phy_tree(physeq)
      
      new_abund_table<-OTU
      
    } else {
      
      list<-na.omit(unique(OTU_taxonomy[,which_level]))
      
      new_abund_table<-NULL
      
      for(i in list){
        
        rt <- na.omit(rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i])
        
        tmp<-data.frame(rowSums(OTU[,rt]))
        
        if(i==""){colnames(tmp)<-c("__Unknowns__")} else {colnames(tmp)<-paste("",i,sep="")}
        
        if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
        
      }
      
    }
    
    OTU<-as.data.frame(as(new_abund_table,"matrix"))
    
    #Convert the data to phyloseq format
    
    OTU = otu_table(as.matrix(OTU), taxa_are_rows = FALSE)
    
    TAX = tax_table(as.matrix(OTU_taxonomy))
    
    SAM = sample_data(SAM)
    
    #reconstruct the phyloseq object
    
    physeq<-NULL
    
    if(which_level=="Otus"){
      
      physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,midpoint(OTU_tree))
      
    } else {
      
      physeq<-merge_phyloseq(phyloseq(OTU),SAM)
      
    }
    
    return(physeq)
    
  }
  
  
  
  
  
  