def convert_rank(df_otu,df_taxa,select_rank,sort_index="median"):
    # ------------------------------- prepare -------------------------------
    ## Notice: Sould confirm that all the taxaid in df_otu is includes in df_taxa, otherwise will give unexpected outputs.
    df_taxa=pair_otu_taxa(df_otu,df_taxa,prompt=False)
    
    # ------------------------------- process the otu table -------------------------------
    
    ## merge the df_otu with df_taxa
    select_df_taxa=df_taxa.reindex(columns=[select_rank])
    df_merge=df_otu.merge(select_df_taxa,left_index=True,right_index=True)
    
    ## convert the None as "NaN_symbol", otherwise will be drop when using groupby
    df_merge[select_rank]=df_merge[select_rank].fillna("NaN_symbol")

    ## convert the original taxa to wanted taxa (just groupby and sum the abs abundance)
    df_coverted=df_merge.groupby(select_rank,as_index=True).sum()
    
    ## sort the df_otu by descending
    df_coverted,sorted_index=sort_by_row(df_coverted,sort_index=sort_index, ascending=False, refresh_index=True, drop_non_numeric=True)
    
    # ------------------------------- process the corresponding df_taxa -------------------------------

    ## drop the col after the select col
    index_select_rank=df_taxa.columns.get_loc(select_rank)
    df_taxa_coverted=df_taxa.iloc[:,0:index_select_rank+1]
    
    ## drop duplicates (If normal, will not have duplicates, but just in case.)
    df_taxa_coverted=df_taxa_coverted.drop_duplicates(subset=[select_rank])
    
    ## match the index of otu table to df_taxa
    df_taxa_coverted.index=df_taxa_coverted[select_rank]
    df_taxa_coverted=sort_by_index_from_other(df_taxa_coverted,sorted_index,axis="index",refresh_index=True)

    ## replace the NaN_symbol, np.nan【caused by df_otu having some taxid not included in df_taxa】 to None
    df_taxa_coverted=df_taxa_coverted.replace(["NaN_symbol",np.nan],None)
    
    # ------------------------------ get the relative abundance --------------------------------------

    ## can get a rel otu table too
    df_coverted_rel=df_coverted/df_coverted.sum()
    
    # ------------------------------ output -----------------------------------
    ## return
    pretty_print("Successfully convert to '{}' rank.".format(select_rank))
    return(df_coverted,df_coverted_rel,df_taxa_coverted)


def select_top(df_otu,df_taxa,sort_index="median",top=10,out_of_top="Other"):
    # ------------------------------- description -------------------------------
    """
    Sometimes we don't need so many taxa, we many want to mark and group the out of top taxa
    
    Arguments:
    1. df_otu: The absolute otu table
    2. df_taxa: The taxa table
    3. sort_index: When select top, which index you want to use? 'median' or 'mean'
    4. top: define the top X
    5. out_of_top: The label of the out of top taxa, default: 'Other'
    """
    # ------------------------------- prepare -------------------------------
    ## check whether the length of df > top
    if top>=len(df_otu):
        raise Exception("The length of df ({}) <= the top({}).".format(len(df_otu),top))
    else:
        print("The length of df ({}) > the top({}), select the first {} and others will be labeled as '{}'".format(len(df_otu),top, top, out_of_top))
    
    ## Notice: Sould confirm that all the taxaid in df_otu is includes in df_taxa, otherwise will give unexpected outputs.
    df_taxa=pair_otu_taxa(df_otu,df_taxa,prompt=False)
    
    
    # ------------------------------- process the otu table -------------------------------
    ## sort the df_otu by descending (Without this, can't to top step)
    df_coverted,sorted_index=sort_by_row(df_otu,sort_index=sort_index, ascending=False, refresh_index=True, drop_non_numeric=True)
    
    ## split the df into df_top (df_coverted)  and df_out_of_top
    
    df_top=df_coverted.loc[0:top]
    df_out_of_top=df_coverted.loc[top:]
    
    ## calculate the sum of out of top, and then concat back to the df_coverted
    df_out_of_top=df_out_of_top.sum()
    df_coverted=df_top
    df_coverted.loc[len(df_coverted)]=df_out_of_top
    
    # ------------------------------- process the corresponding df_taxa -------------------------------
    
    ## match the index of otu table to df_taxa
    df_taxa_coverted=sort_by_index_from_other(df_taxa,sorted_index,axis="index",refresh_index=True)
    df_taxa_coverted=df_taxa_coverted.loc[0:top]
    
    ## add a taxa row as out_of_top
    df_taxa_coverted.loc[len(df_taxa_coverted)]=out_of_top
    
    # ------------------------------ get the relative abundance --------------------------------------
    
    ## can get a rel otu table too
    df_coverted_rel=df_coverted/df_coverted.sum()
    
    # ------------------------------ output -----------------------------------
    pretty_print("Successfully select the top {}.".format(top))
    return(df_coverted,df_coverted_rel,df_taxa_coverted)

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
    
    # ------------------- use the index in df_describe to sort the original df
    index_after_sort=df_describe.index
    df=df.reindex(index=index_after_sort)
    return(df)



