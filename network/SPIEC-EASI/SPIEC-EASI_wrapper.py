## description

'''
otu_path: The path of otu table
meta_path: The path of meta data
ncores: The cores we used to do the jobs
top (optional): Select the first N OTUs from the original otu table to do the jobs. It's for testing.
'''

## basic function

def pretty_print(t):
    print(t+"\n")
    
import sys
otu_path=sys.argv[1]
meta_path=sys.argv[2]
ncores=int(sys.argv[3])
pretty_print('''
The otu_path is {},
the meta_path is {}.'''.format(otu_path,meta_path))
pretty_print("Use {} cores to do the jobs.".format(ncores))

try:
    top=sys.argv[4]
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


## help print
def pretty_print(t):
    print(t+"\n")
    
## help pair the index of df_meta and df_otu
def pair_otu_meta(df_otu,df_meta,how="left"):
    """
    phyloseq not require the df_otu and df_taxa to have the same samples.
    
    But sometimes we will extract data from the built phyloseq object.

    So for safety reason, we'd better pair them first.
    """
    if how=="left":
        pretty_print("Use the index of df_otu.")
        df_meta=df_meta.reindex(index=df_otu.columns)
    elif how=="right":
        pretty_print("Use the index of df_meta.")
        df_otu=df_otu.reindex(columns=df_meta.index)
    return(df_otu,df_meta)

## help pair the index of df_otu and df_taxa
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
            """The original length of otu_table is {},
            that of taxa_table is {},
            match the index of taxa_table based on otu_table automaticlly.""".format(len(df_otu), len(df_taxa))
        )
        
    # ------------------------------ match index ------------------------------
    df_taxa = df_taxa.reindex(index=df_otu.index)
    
    # ------------------------------ replace the np.nan and pd.NA by None (p2r only accept None, not np.nan or pd.NA) ------------------------------
    #### df_taxa = df_taxa.replace(to_replace=[np.nan,pd.NA], value=None)
    df_taxa=df_taxa.mask(df_taxa.isna(),other=None)
    return df_taxa



def build_phy(df, df_meta, dft=NA,taxa_are_rows=True, how_to_pair_meta="left", prompt=True):
    ## R wrapper for building a phyloseq object from df_otu, df_meta and df_taxa
    string = """
    build_phy=function(matrix,index,df_meta, dft=NA, taxa_are_rows=T){
        ## rebuid df
        df = data.frame(matrix)
        rownames(df) = index
        colnames(df) = rownames(df_meta)
        
        ## make phy obj
        
        #### include otu info
        otu_obj=otu_table(df,taxa_are_rows=taxa_are_rows)
        
        #### include meta info 【will create a new columns call 'sample_name', for used in case】
        df_meta["sample_name"]=rownames(df_meta)
        meta_obj=sample_data(df_meta)
        
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

    ## pair the index of df_meta and df_otu
    df, df_meta=pair_otu_meta(df,df_meta, how=how_to_pair_meta)
        
    ## pair the index of df_otu and df_taxa 
    if type(dft) != type(NA):
        dft = pair_otu_taxa(df, dft,prompt=prompt)

    phy_obj = my_phy_sub.build_phy(df.values, df.index, df_meta, dft, taxa_are_rows)
    return phy_obj

otu_sep="\t"
meta_sep=","
def read_otu(otu_path,sep=otu_sep):
    df_otu=pd.read_csv(otu_path,index_col=0,sep=sep)
    return(df_otu)

def read_meta(meta_path,sep=meta_sep):
    df_meta=pd.read_csv(meta_path,index_col=0,sep=sep)
    return(df_meta)

def spiec_easi(phy_obj,method="mb",lambda_min_ratio=1e-3,nlambda=30,sel_criterion="bstars",rep_num=50,seed=1000,ncores=40):
    string = """
    spiec_easi=function(phy_obj,method="mb",lambda.min.ratio=1e-3,nlambda=30,sel.criterion="bstars",rep.num=50,seed=1000,ncores=40,out_dir="."){
        ## import package
        library(igraph)
        library(Matrix)
        library(SpiecEasi)
        
        ## run spiec and get output
        pargs = list(rep.num=rep.num, seed=seed, ncores=ncores)
        spiec_res=spiec.easi(phy_obj,method=method,lambda.min.ratio=lambda.min.ratio,nlambda=nlambda,sel.criterion=sel.criterion,pulsar.select=TRUE,pulsar.params=pargs)
        
        
        bm=symBeta(getOptBeta(spiec_res), mode="maxabs")
        
        
        ## graph
        #### get graph
        graph_res <- adj2igraph(bm,
                  diag=F,
                  vertex.attr = list(label=taxa_names(phy_obj)))
        
        ## save graph
        print("Save graph.")
        write.graph(graph_res, "weight_res.gml", format = "gml")
    }
"""
    my_spiec_easi_sub = STAP(string, "my_spiec_easi_sub")
    spiec_res=my_spiec_easi_sub.spiec_easi(phy_obj,method=method,lambda_min_ratio=lambda_min_ratio,nlambda=nlambda,sel_criterion=sel_criterion,rep_num=rep_num,seed=seed,ncores=ncores)
    
phyloseq = importr("phyloseq")



df_otu=read_otu(otu_path)
df_meta=read_meta(meta_path)


if type(top) != type(None):
    pretty_print("Use the top {} taxa to do the SPIEC-EASI.".format(top))
    df_otu = df_otu.iloc[0:top]
else:
    pretty_print("Use all the taxa to do the SPIEC-EASI.")

my_phy = build_phy(df_otu, df_meta)

spiec_easi(my_phy, ncores=ncores)

pretty_print("All finished.")