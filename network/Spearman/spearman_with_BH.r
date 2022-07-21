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
    cat(tt,"\n\n")
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
library(tidyverse)
library(WGCNA)
library(igraph)




## get otu_table
otu_table=read.csv(otu_path,sep="\t",header=T,comment.char="~",row.names = 1)

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





## adjust the p-value and then get the index needed to mask.

#### get graph
g_p=graph_from_adjacency_matrix(as.matrix(res_p),mode="undirected",weight=T,diag=F,add.colnames="label")


g_cor=graph_from_adjacency_matrix(as.matrix(res_cor),mode="undirected",weight=T,diag=F,add.colnames="label")


#### check the length of g_p and g_cor
if(length(E(g_p))==length(E(g_cor))){
    pretty_print("The length of p-value equals to corr, continue.")
}else{
    stop("The lenght of p-value NOT equals to corr, stop.")
}

#### the g_p is just for converting the adjacent p matrix to data.frmae format

E(g_p)$weight=p.adjust(E(g_p)$weight,"BH")

#### get the mask_index and del_index

mask_index=(E(g_p)$weight> 0.05)|(is.na(E(g_p)$weight))
del_index=which(mask_index)




## polish the corr graph


#### mask the p
E(g_cor)$weight[mask_index]=0


#### delete the edges with weight==0, becasue we don't want to save redundant info.
g_cor=delete.edges(g_cor, del_index)

#### delete the nodes without edges
del_index2=which(degree(g_cor)==0)
g_cor=delete_vertices(g_cor,del_index2)







## save the gml
pretty_print("Save gml.")
write.graph(g_cor, "weight_res.gml", format = "gml")