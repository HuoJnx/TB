## description
#### Here is an advance spearman correlation network constructor.
#### Inputs: ① otu_table ② p threshold ③ top (Only for debug)
#### Processing procedures are as follows:
#### 1. Use WGCNA to get the spearman correlation matrix including coefficient and p-value.
#### 2. Use igraph to convert the coefficient and p matrix to long format.
#### 3. Use BH method to adjust the p-value, then mask the coefficient > p threshold as 0.
#### 4. Save the new adjacent matrix, node_res, edge_res and gml format network.

## self function
pretty_print=function(tt){
    print(tt)
    print("")
}




## get arguments
myargs = commandArgs(trailingOnly=TRUE)

#### get otu path
otu_path=myargs[1]
pretty_print(sprintf("The path of the otu table is %s.",otu_path))

#### get p_thresh
p_thresh=as.double(myargs[2])
if(is.na(p_thresh)){
    p_thresh=0.05
    pretty_print("Can't get the customized threshold of p-value, use 0.05.")
    pretty_print("Additionally, the p-value will be correted by 'BH' method in p.adjust().")
}
pretty_print(sprintf("The threshold of the p-value is %s.",p_thresh))
#### get the top
top=as.integer(myargs[3])
if(is.na(top)){
    pretty_print("Use the whole otu table to construct the network.")
}else{
    pretty_print(sprintf("Construct the table using top %s OTUs.",top))
}





## import package
library(WGCNA)
library(igraph)




## get otu_table
otu_table=read.csv(otu_path,sep="\t",row.names = 1)

if(!is.na(top)){
    otu_table=otu_table[0:top,]
}





## get cor and p res
pretty_print("Get spearman matrix from corAndPvalue function in WGCNA package. It will take a while.")
res=corAndPvalue(t(otu_table),method="spearman",alternative="two.side")
res_p=res$p
res_cor=res$cor

#### when convert from adjacent to graph, if value==0, the edge will be dropped, so should assign a very small number
res_p[res_p==0]=1e-10
res_cor[res_cor==0]=1e-10

#### save the original cor
pretty_print("Save the original corr matrix before adjustment.")
write.table(res_cor,"original_cor.tsv",sep="\t",col.names = NA)





## adjust the p-value and then get the index needed to mask.
#### the g_p is just for converting the adjacent p matrix to data.frmae format
g_p=graph_from_adjacency_matrix(as.matrix(res_p),mode="undirected",weight=TRUE,diag=FALSE,add.colnames=NA)
E(g_p)$weight=p.adjust(E(g_p)$weight,"BH")

#### get the mask_index and del_index
mask_index=(E(g_p)$weight> 0.05)
del_index=which(mask_index)

#### save the adjust p matrix
pretty_print("Save the adjust-p (BH) matrix.")
new_p=get.adjacency(g_p,attr="weight")
new_p=data.frame(as.matrix(new_p))
rownames(new_p)=rownames(res_p)
colnames(new_p)=colnames(res_p)
write.table(new_p,"adjust_p.tsv",sep="\t",col.names = NA)

## get the corr graph
g_cor=graph_from_adjacency_matrix(as.matrix(res_cor),mode="undirected",weight=TRUE,diag=FALSE,add.colnames=NA)

#### mask the p
E(g_cor)$weight[mask_index]=0





## rebuild the new adjacent matrix
pretty_print("Save the adjaceny_weight matrix after p adjustment.")
new_adj=get.adjacency(g_cor,attr="weight")
new_adj=data.frame(as.matrix(new_adj))
rownames(new_adj)=rownames(res_cor)
colnames(new_adj)=colnames(res_cor)
write.table(new_adj,"adjacency_weigth.tsv",sep="\t",col.names = NA)







## save others

#### delete the edges with weight==0, becasue we don't want to save redundant info.
g_cor=delete.edges(g_cor, del_index)

#### save the gml
pretty_print("Save gml.")
write.graph(g_cor, "weight_res.gml", format = "gml")

#### save node and edge table
df_cors=get.data.frame(g_cor,what="both")


#### save the node table
pretty_print("Save node_res.")
df_cor_v=df_cors$vertices
df_cor_v$id=rownames(df_cor_v)
df_cor_v$label=colnames(res_cor)
write.table(df_cor_v,"node_res.tsv",sep="\t",row.names=F)

#### save the edge table
df_cor_e=df_cors$edges
colnames(df_cor_e)=c("source","target","weight")
write.table(df_cor_e,"edge_res.tsv",sep="\t",row.names=F)