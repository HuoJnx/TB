def plot_alpha_diversity(df_otu,df_meta,df_taxa,df_taxa_id,rank,group_col="group",stat_m="kruskal",  pair_m="dunn", adj_m="BH",label_samples= True):
    # -------------------------------- description --------------------------------
    """
    Plot the alpha diversity, including the 'simpson' and 'shannon'.
    
    Arguments:
    1. df, original df_otu, absolute value.
    2. df_meta, dataframe which index is sample names and a col indicating the group info of the samples
    3. rank, the rank you want to plot.
    4. group_col, the col name in df_meta indicating the group info.
    5. stat_m, the test used for general comparesion.
    6. pair_m, the test used for pair_wise comparesion.
    7. adj_m, the adjustment method for pair_wise comparesion.
    8. label_samples, whether to plot the sample names in the jitter plot.
    """
    # -------------------------------- prepare ----------------------------------
    
    ## index list
    index_list=["simpson","shannon"]
    fig_list=[]
    df_list=[]
    ## convert_rank
    df_converted,_,df_taxa_converted,_ = convert_rank(df_otu, df_taxa, df_taxa_id,"order")
    
    ## get relative abundance
    dfr=df_converted/df_converted.sum()
    ## get alpha diversity for each sample【a DF with only one col】
    for index in index_list:

        res = dfr.apply(func=vegan.diversity, axis=0, index=index).T
        
        #### merge the alpha diversity for each sample with the group information
        res_df = res.merge(df_meta, left_index=True, right_index=True)
        
        index = index.capitalize()
        res_df.columns = [index, group_col]
        
        #### description results

        df_descr=description(res_df,group_col,index)
        
        #### plot box and jitter
        fig = (
            gp.ggplot(res_df)
            + gp.aes_string(x=group_col, y=index)
            + gp.geom_jitter()
            + gp.geom_boxplot(outlier_shape=NA, fill=NA)
        )
        #### plot hypothesis test
        res_list = plot_hypothesis(fig, res_df, group_col, index, stat_m=stat_m, pair_m=pair_m, adj_m=adj_m)
        fig = get_rlist(res_list,"fig")
        
        #### plot labels of samples
        if label_samples:
            fig = fig + ggrepel.geom_text_repel(gp.aes(label=df_meta.index))
        
        #### plot theme
        fig = mytheme(fig)
        fig_list.append(fig)
        df_list.append(df_descr)
    return(fig_list,df_list)
vegan = importr("vegan")