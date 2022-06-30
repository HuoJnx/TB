def sort_by_row(df,sort_index="median",ascending=False):
    # ------------------- description -----------------------
    """
    In recent practice, we often need to sort the rows by the summary value of the rows, including mean, median and so on.
    However, nowadays pandas have no this kind of function. So I make myself.
    
    It supports the input dataframe with columns in 'object' dtype, it will just ignore it.
    
    Arguments:
    1. df: the input dataframe.
    
    2. index: the index you used to sort the df, the supported index are "mean" and "median".
    
    3. ascending: For sort.
    """

    # ------------------- get describe and sort ----------------------
    if sort_index=="median":
        df_describe=df.median(axis=1,numeric_only=True)
    elif sort_index=="mean":
        df_describe=df.mean(axis=1,numeric_only=True)
    else:
        raise Exception("Should be 'median' or 'mean'.")
    df_describe=df_describe.sort_values(ascending=ascending)
    
    # ------------------- use the index in df_describe to sort the original df ----------------------
    index_after_sort=df_describe.index
    df=df.reindex(index=index_after_sort)
    return(df)

def convert_rank(df_otu,df_taxa,df_taxa_id,select_rank,sort_index="median"):
    # ------------------------------- prepare -------------------------------
    
    ## get the name and id for the select rank
    select_df_taxa_name=df_taxa.reindex(columns=[select_rank])
    select_df_taxa_name.columns=["taxa_name"]
    select_df_taxa_id=df_taxa_id.reindex(columns=[select_rank])
    select_df_taxa_id.columns=["#OTU ID"]
    # ------------------------------- process the otu table -------------------------------
    
    ## merge the df_otu with taxa name and id info
    df_merge=df_otu.merge(select_df_taxa_name,left_index=True,right_index=True,how="left")
    df_merge=df_merge.merge(select_df_taxa_id,left_index=True,right_index=True,how="left")
    
    ## deal with the unclassified
    df_merge["taxa_name"]=df_merge["taxa_name"].fillna("unclassified")
    df_merge["#OTU ID"]=df_merge["#OTU ID"].fillna(0)

    ## because the df_otu need to be sort by median, so the taxa info should be splited into another dataframe first
    df_taxa_info=df_merge[["taxa_name","#OTU ID"]]
    df_taxa_info=df_taxa_info.drop_duplicates()
    df_merge=df_merge.drop(columns=["#OTU ID"])
    ## convert the original taxa to wanted taxa (just groupby and sum the abs abundance)
    df_coverted=df_merge.groupby("taxa_name",as_index=True).sum()
    
    ## sort the df_otu by descending
    df_coverted=sort_by_row(df_coverted,sort_index=sort_index, ascending=ascending)
    
    ## merget again with the taxa info
    df_coverted=df_coverted.merge(df_taxa_info,left_on="taxa_name",right_on="taxa_name",how="left")
    
    # ------------------------------- process the corresponding df_taxa -------------------------------

    ## drop the col after the select col
    index_select_rank=df_taxa.columns.get_loc(select_rank)
    df_taxa_coverted=df_taxa.iloc[:,0:index_select_rank+1]
    
    ## drop duplicates
    df_taxa_coverted=df_taxa_coverted.drop_duplicates(subset=[select_rank])
    
    ## match the index of otu table to df_taxa
    df_taxa_coverted=df_coverted.reindex(columns=["taxa_name","#OTU ID"]).merge(df_taxa_coverted,left_on="taxa_name",right_on=select_rank,how="left")
    
    ## final polishment for df_taxa
    df_taxa_coverted=df_taxa_coverted.set_index("#OTU ID")
    df_taxa_coverted=df_taxa_coverted.drop(columns=["taxa_name"])
    
    # ------------------------------- process the otu table -------------------------------
    ## final polishment for df_otu
    df_coverted=df_coverted.set_index("#OTU ID")
    df_coverted=df_coverted.drop(columns="taxa_name")
    
    # ------------------------------ get the relative abundance --------------------------------------

    ## can get a rel otu table too
    df_coverted_rel=df_coverted/df_coverted.sum()
    
    # ------------------------------ output -----------------------------------
    ## return
    pretty_print("Successfully convert to '{}' rank.".format(select_rank))
    return(df_coverted,df_coverted_rel,df_taxa_coverted)