## description

'''
Procedure:
1. Build a phyloseq object from otu_path, meta_path and taxa_path.
2. Run SPIEC-EASI using phyloseq object.
3. Output the node_res, edge_res, and gml format for further analysis like cytoscape.


Arguments:

1. otu_path: The path of otu table
2. meta_path: The path of meta data
3. taxa_path: The path of taxa information
4. ncores: The cores we used to do the jobs
5. top (optional): Select the first N OTUs from the original otu table to do the jobs. It's for testing.
'''

## basic function

def pretty_print(t):
    print(t+"\n")
    
import sys
otu_path=sys.argv[1]
meta_path=sys.argv[2]
taxa_path=sys.argv[3]
ncores=int(sys.argv[4])
pretty_print('''
The otu_path is {},
the meta_path is {},
the taxa_path is {}.'''.format(otu_path,meta_path,taxa_path))
pretty_print("Use {} cores to do the jobs.".format(ncores))

try:
    top=sys.argv[5]
    top=int(top)
    pretty_print("Top = {}".format(top))
except Exception as e:
    pretty_print("Error: {}, so top = None".format(e))
    top=None
    


import pandas as pd
from pandas import DataFrame as DF, Series


## rpy2

#### rpy2 util

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
from rpy2.robjects.conversion import rpy2py as r2p

#### some r function
rprint = ro.globalenv.find("print")
from rpy2.robjects.packages import data as rdata




def pair_otu_taxa(df_otu, df_taxa,prompt=True):
    """
    The phyloseq require the df_otu and df_taxa should have the same index, so should pair them.
    
    Before using the function, should confirm that all the taxaid in df_otu is includes in df_taxa, otherwise will give unexpected outputs.
    
    The function will first use the index in df_otu to reindex the df_taxa, so you should confrim that the index in df_otu is taxid.
    
    The function will then replace the np.nan and pd.NA by None, which is accepted by p2r function in rpy2.
    """
    # ------------------------------ prompt ------------------------------
    if prompt:
        pretty_print(
            "The original length of otu_table is {}, that of taxa_table is {}, match the index of taxa_table based on otu_table automaticlly.".format(
                len(df_otu), len(df_taxa)
            )
        )
        
    # ------------------------------ match index ------------------------------
    df_taxa = df_taxa.reindex(index=df_otu.index)
    
    # ------------------------------ replace the np.nan and pd.NA by None (p2r only accept None, not np.nan or pd.NA) ------------------------------
    df_taxa=df_taxa.mask(df_taxa.isna(),other=None)
    return df_taxa


def build_phy(df, dfg, dft=NA, taxa_are_rows=True, prompt=True):
    string = """
    build_phy=function(matrix,index,dfg, dft=NA, taxa_are_rows=T){
        ##rebuid df
        df = data.frame(matrix)
        rownames(df) = index
        colnames(df) = rownames(dfg)

        ## make phy obj

        #### include otu info
        otu_obj=otu_table(df,taxa_are_rows=taxa_are_rows)

        #### include meta info 【will create a new columns call 'sample_name', for used in case】
        dfg["sample_name"]=rownames(dfg)
        meta_obj=sample_data(dfg)

        #### include taxa info
        if (class(dft)!=class(NA)){
            matrix_t=as.matrix(dft)#without this step will get error
            tax_obj=tax_table(matrix_t)
            phy_obj=phyloseq(otu_obj,meta_obj,tax_obj)
        }else{
            phy_obj=phyloseq(otu_obj,meta_obj)
        }

        return(phy_obj)
    }
"""
    my_phy_sub = STAP(string, "_")
    if type(dft) != type(NA):
        dft = pair_otu_taxa(df, dft,prompt=prompt)
    phy_obj = my_phy_sub.build_phy(df.values, df.index, dfg, dft, taxa_are_rows)
    return phy_obj

def read_data(otu_path,meta_path,taxa_path,sep1="\t",sep2=",",sep3=","):
    df_otu=pd.read_csv(otu_path,index_col=0,sep=sep1)
    df_meta=pd.read_csv(meta_path,index_col=0,sep=sep2)
    df_taxa=pd.read_csv(taxa_path,index_col=0,sep=sep3)
    return(df_otu,df_meta,df_taxa)

def spiec_easi(phy_obj,method="mb",lambda_min_ratio=1e-3,nlambda=30,sel_criterion="bstars",rep_num=50,seed=1000,ncores=40):
    string = """
    spiec_easi=function(phy_obj,method="mb",lambda.min.ratio=1e-3,nlambda=30,sel.criterion="bstars",rep.num=50,seed=1000,ncores=40,out_dir="."){

        ## run spiec
        pargs = list(rep.num=rep.num, seed=seed, ncores=ncores)
        spiec_res=spiec.easi(phy_obj,method=method,lambda.min.ratio=lambda.min.ratio,nlambda=nlambda,sel.criterion=sel.criterion,pulsar.select=TRUE,pulsar.params=pargs)
        taxid=rownames(otu_table(phy_obj))

        ## graph
        graph_res=adj2igraph(getRefit(spiec_res), vertex.attr = list(label = taxid ))

        ## node
        print("Get node table.")
        node_res= data.frame(id = as.vector(V(graph_res)), label = V(graph_res)$label)
        write.table(node_res,"node_res.tsv",sep="\t", row.names=F)
        ## edge
        print("Get edge table.")
        sebeta = symBeta(getOptBeta(spiec_res), mode = "maxabs")
        edge_res = summary(sebeta)
        #### give header
        edge_res = data.frame(edge_res)
        names(edge_res) = c("source", "target", "weight")
        #### sort
        edge_res = edge_res[order(edge_res$source, edge_res$target),]
        #### save
        write.table(edge_res,"edge_res.tsv",sep="\t", row.names=F)
        ## adjacency_weight matrix
        print("Get adjacency weigth matrix.")
        E(graph_res)$weight = edge_res$weight
        adjacency_weight = as.data.frame(as.matrix(get.adjacency(graph_res, attr = "weight")))
        rownames(adjacency_weight) = taxid
        colnames(adjacency_weight) = taxid
        write.table(adjacency_weight,"adjacency_weight.tsv",sep="\t",col.names = NA)
        ## save graph
        print("Save graph.")
        write.graph(graph_res, "weight_res.gml", format = "gml")
        return(list(node = node_res,edge = edge_res, adjacency_weight = adjacency_weight))
    }
"""
    my_spiec_easi_sub = STAP(string, "my_spiec_easi_sub")
    spiec_res=my_spiec_easi_sub.spiec_easi(phy_obj,method=method,lambda_min_ratio=lambda_min_ratio,nlambda=nlambda,sel_criterion=sel_criterion,rep_num=rep_num,seed=seed,ncores=ncores)
    
phyloseq = importr("phyloseq")
spieceasi = importr("SpiecEasi")
igraph = importr("igraph")
Matrix = importr("Matrix")


df_otu, df_meta, df_taxa = read_data(otu_path, meta_path, taxa_path)

if type(top) != type(None):
    pretty_print("Use the top {} taxa to do the SPIEC-EASI.".format(top))
    df_otu = df_otu.iloc[0:top]
else:
    pretty_print("Use all the taxa to do the SPIEC-EASI.")

my_phy = build_phy(df_otu, df_meta, df_taxa)

spiec_easi(my_phy, ncores=ncores)