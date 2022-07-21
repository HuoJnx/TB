def select_top(df_otu,df_taxa,sort_index="median",top=10,out_of_top="Others"):
    
    # ------------------------------- prepare -------------------------------
    ## check whether the length of df > top
    if top>len(df_otu):
        print ("The length of df ({}) < the top({}), return original.".format(len(df_otu),top))
        return(df_otu,df_taxa)
    else:
        print("The length of df ({}) >= the top({}), select the first {} and others will be labeled as '{}'".format(len(df_otu),top, top, out_of_top))
    
    
    # ------------------------------- process the otu table -------------------------------
    
    ## sort the df_otu by descending (Without this, can't to top step)
    df_coverted=sort_by_row(df_otu,sort_index=sort_index, ascending=False)
    
    ## split the df into df_top (df_coverted)  and df_out_of_top
    
    df_top=df_coverted.iloc[0:top]
    df_out_of_top=df_coverted.iloc[top:]
    
    ## calculate the sum of out of top, and then concat back to the df_coverted
    df_out_of_top=df_out_of_top.sum()
    df_coverted=df_top
    OTU_ID_others=99999999
    df_coverted.loc[OTU_ID_others]=df_out_of_top
    
    # ------------------------------- process the corresponding df_taxa -------------------------------
    df_taxa_coverted=df_taxa.reindex(index=df_coverted.index)
    df_taxa_coverted.loc[OTU_ID_others]=out_of_top
    
    # ------------------------------ output -----------------------------------
    pretty_print("Successfully select the top {}.".format(top))
    return(df_coverted,df_taxa_coverted)