## Description
#### 1. It masks the adjacent matrix with p>0.05 as 0.
#### 2. It outputs a gml network file, edge_res and node_res for further analysis like cytoscape. 


## import functions

#### pandas
import pandas as pd
from pandas import DataFrame as DF, Series

#### rpy2
import rpy2.robjects as ro

#### dataframe conversion
from rpy2.robjects import Formula
from rpy2.robjects import NA_Logical as NA
from rpy2.rinterface import NULL as NULL
from rpy2.robjects import r, rl
from rpy2.robjects.packages import importr

from rpy2.robjects.packages import STAP


#### dataframe conversion
from rpy2.robjects import pandas2ri

pandas2ri.activate()
from rpy2.robjects.pandas2ri import py2rpy as p2r

#### self
def pretty_print(t):
    t=str(t)
    print(t+"\n")

    
    
    
    
    
## parse input
import sys

weight_path=sys.argv[1]
p_path=sys.argv[2]
pretty_print('''
The path of the weight_df is {},

the path of the p_df is {},
'''.format(weight_path,p_path))

try:
    p_thresh=sys.argv[3]
    p_thresh=float(p_thresh)
except Exception as e:
    p_thresh=0.05
    pretty_print(e)
    pretty_print("No P threshold input, use default threshold=={}".format(p_thresh))
pretty_print("P threshold = {}".format(p_thresh))








## read data
weight_df=pd.read_csv(weight_path,sep="\t",index_col=0)
p_df=pd.read_csv(p_path,sep="\t",index_col=0)




## get the masked adjacent matrix
#### mask the p>threshold
pretty_print("Mask the corr with p>{}".format(p_thresh))
masked_weight_df=weight_df.mask(p_df>p_thresh,other=0)
#### save the masked matrix
masked_weight_df.index.name=""
masked_weight_df.to_csv("significant_adjacency_weight.tsv",sep="\t")





## get the gml format and save
pretty_print("Save weight_res.gml from masked adjacent matrix.")
igraph=importr("igraph")
string = """
adj2gml=function(adj_df,file_name="weight_res.gml"){
    library(igraph)
    g = graph_from_adjacency_matrix(as.matrix(adj_df), mode = 'undirected', weighted = TRUE, diag = FALSE)
    write.graph(g,"weight_res.gml",format="gml")
}
"""
adj2gml = STAP(string, "_")
adj2gml.adj2gml(masked_weight_df)




## get the node_res
pretty_print("Save node_res.tsv")
node_label=masked_weight_df.index
node_id=[e for e in range(len(node_label))]
node_res=DF({"id":node_id,"label":node_label})
node_res.to_csv("node_res.tsv",sep="\t",index=False)





## get the edge_res (melt from adjacent matrix)
#### melt
pretty_print("Save edge_res.tsv")
edge_res=masked_weight_df.copy()
edge_res.index=node_res["id"]
edge_res.columns=node_res["id"]
edge_res.index.name="source"
edge_res.columns.name=""
edge_res=edge_res.reset_index()
edge_res=edge_res.rename(columns={"index":"source"})
edge_res=edge_res.melt(id_vars=["source"],var_name=["target"],value_name="weight")

#### filter the existed one
edge_res=edge_res.query("weight!=0")

#### sort by source and target
edge_res=edge_res.sort_values(["source","target"])

#### save
edge_res.to_csv("edge_res.tsv",sep="\t",index=False)