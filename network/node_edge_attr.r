## description
#### Get attributes for nodes including degree, closeness, betweenness_centrality and eigenvector_centrality
#### Arguementss:
#### 1. Path of network
#### 2. Path of df_taxa
#### 3. annotation rank
#### 4. output directory






## self function
pretty_print=function(tt){
    cat(tt,"\n\n")
}

options(warn=-1)



## hello prompt
pretty_print("Hello, here is a R script for exploring node and edge attributes.")





## get arguments
myargs = commandArgs(trailingOnly=TRUE)

#### get network path (gml)
net_path=myargs[1]
pretty_print(sprintf("The path of the network is %s",net_path))

#### get taxa table
taxa_path=myargs[2]
pretty_print(sprintf("The path of the df_taxa is %s",taxa_path))

#### get anotation rank
annotation_rank=myargs[3]
if(is.na(annotation_rank)){
    annotation_rank="species"
    pretty_print(sprintf("Not input the annotation rank, use '%s' as default.",annotation_rank))
}
pretty_print(sprintf("The annotation rank is '%s'",annotation_rank))

#### get the output directory
out_dir=myargs[4]
dir.create(out_dir,showWarnings=F)
pretty_print(sprintf("The output directory is %s.",out_dir))




## library
library(tidyverse,warn.conflicts=F)
library(igraph,warn.conflicts=F)
library(GGally,warn.conflicts=F)



## read data

#### read taxonomy data
pretty_print("Read taxonomy data.")
df_taxa=read_tsv(taxa_path)
df_taxa=df_taxa%>%mutate(`#OTU ID`=as.character(`#OTU ID`))

#### read net
pretty_print("Read network file in gml format.")
net=read_graph(file=net_path,format="gml")




## polish 

#### Use the abs weight
E(net)$corr = E(net)$weight
E(net)$weight = abs(E(net)$weight)

#### get data.frame
pretty_print("Make df_node.")
label=vertex_attr(net,"label")
if(is.null(label)){
    label=vertex_attr(net,"name")
}
df_node=tibble(label=as.character(vertex_attr(net,"label")),id=as.vector(V(net)))

#### merge
df_node=left_join(df_node,df_taxa,by=c("label"="#OTU ID"))


## save the new graph with rank annotation
#### annotation
V(net)$superkingdom=df_node$superkingdom
V(net)$kingdom=df_node$kingdom
V(net)$phylum=df_node$phylum
V(net)$class=df_node$class
V(net)$order=df_node$order
V(net)$family=df_node$family
V(net)$genus=df_node$genus
V(net)$species=df_node$species
#### save
pretty_print("Save new graph with rank annotation.")
netname_vector=strsplit(basename(net_path), "\\.")[[1]]

new_netname=paste0(netname_vector[1],"_with_taxa",".",netname_vector[2])
graph_path=file.path(out_dir,new_netname)
write.graph(net,graph_path,format="gml")

#### only keep the id, label and needed rank
df_node=df_node%>%select(id,label,all_of(annotation_rank))
df_rank=df_node




## node attribute
#### degree
df_node["degree"]=degree(net)
df_node["neighbor"]=graph.knn(net, V(net),weights = NA )$knn
df_node["weight_degree"]=strength(net)
df_node["weight_neighbor"]=graph.knn(net, V(net))$knn

#### closeness
pretty_print("calculate closeness centrality.")
df_node["closeness_centrality"] = closeness(net)

#### betweenness_centrality
pretty_print("calculate betweenness centrality.")
df_node["betweenness_centrality"]=betweenness(net)

#### eigenvector_centrality
pretty_print("calculate eigenvector centrality.")
df_node["eigenvector_centrality"]=evcent(net)$vector


#### save results

###### df_node
pretty_print("Save df_node.")
df_name="df_node_summary.tsv"
df_path=file.path(out_dir,df_name)
write_tsv(df_node,file=df_path)

###### summary
pretty_print("Save summary for node attributes.")

df_stat=df_node%>%select(-id,-label,-all_of(annotation_rank))
for (e in colnames(df_stat)){
    summ = summary(df_stat[[e]])
    file_name=sprintf("Node_summary_%s.txt",e)
    file_path=file.path(out_dir,file_name)
    capture.output(summ,file=file_path)
}

###### pair plot
######## scatter
pretty_print("Save plots for node attributes.")
fig=df_stat%>%ggpairs()+theme(axis.text.x = element_text(angle = 90))
fig_name="Node_pairplot.pdf"
fig_name=file.path(out_dir,fig_name)
ggsave(fig_name,fig,width = 16,height = 16)

######## density
func=function(){
    fig=df_stat%>%ggpairs(lower = list(continuous = "density"))+theme(axis.text.x = element_text(angle = 90))
}

fig = try(func(),silent = T)
if ("try-error" %in% class(fig)){
    pretty_print("Error in density plot, skip.")
}else{
    fig_name="Node_pairplot_density.pdf"
    fig_path=file.path(out_dir,fig_name)
    ggsave(fig_path,fig,width = 16,height = 16)
}



## edge attribute
#### weigth & betweenness_centrality
pretty_print("Make df_edge.")
df_edge = as_tibble(as_edgelist(net))
colnames(df_edge)=c("source","target")
df_edge["weight"]=E(net)$weight
df_edge["betweenness_centrality"]=edge.betweenness(net)

#### merge taxa information
###### merge for source node
df_edge=left_join(df_edge,df_rank,by=c("source"="id"))
new_name1=paste0("source","_",annotation_rank)
rename_vector=setNames(annotation_rank,new_name1)
df_edge=df_edge%>%rename(rename_vector,source_label=label)

###### merge for target node
df_edge=left_join(df_edge,df_rank,by=c("target"="id"))
new_name2=paste0("target","_",annotation_rank)
rename_vector=setNames(annotation_rank,new_name2)
df_edge=df_edge%>%rename(rename_vector,target_label=label)

#### reoder the df_edge
col_order=c("source","target","source_label","target_label",new_name1,new_name2,"weight","betweenness_centrality")
df_edge=df_edge%>%dplyr::select(all_of(col_order))

#### save results

###### df_edge
pretty_print("Save df_edge.")
df_name="df_edge_summary.tsv"
df_path=file.path(out_dir,df_name)
write_tsv(df_edge,file=df_path)

###### summary
pretty_print("Save summary for edge attributes.")

df_stat=df_edge%>%select(weight,betweenness_centrality)
for (e in colnames(df_stat)){
    summ = summary(df_stat[[e]])
    file_name=sprintf("Edge_summary_%s.txt",e)
    file_path=file.path(out_dir,file_name)
    capture.output(summ,file=file_path)
}

###### pair plot
pretty_print("Save plots for edge attributes.")
######## sactter
fig=df_stat%>%ggpairs()
fig_name="Edge_pairplot.pdf"
fig_name=file.path(out_dir,fig_name)
ggsave(fig_name,fig,width = 8,height = 8)

######## density
func=function(){
    fig=df_stat%>%ggpairs(lower = list(continuous = "density"))
}

fig = try(func(),silent = T)
if ("try-error" %in% class(fig)){
    pretty_print("Error in density plot, skip.")
}else{
    fig_name="Edge_pairplot_density.pdf"
    fig_path=file.path(out_dir,fig_name)
    ggsave(fig_path,fig,width = 8,height = 8)
}


## finish prompt
pretty_print("All finished.")