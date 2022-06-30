## help print
def pretty_print(t):
    print(t+"\n")
    
## help pair the index of otu_table and taxa_table
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
## R wrapper for building a phyloseq object from otu_df, meta_df and taxa_df

def build_phy(df, dfg, dft=NA, taxa_are_rows=True, prompt=True):
    string = """
    build_phy=function(matrix,index,dfg, dft=NA, taxa_are_rows=T){
        ## rebuid df
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