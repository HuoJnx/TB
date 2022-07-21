string="""
extract_phy=function(phy_obj){
    ## extract df_otu and transport
    df_otu=as.data.frame(as.matrix(otu_table(phy_obj)))
    df_otu_t=t(df_otu)
    
    ## extract df_meta
    df_meta=as.data.frame(as.matrix(sample_data(phy_obj)))
    group_col=colnames(df_meta)[1]
    df_group=df_meta[group_col]
    
    return(list(df_otu_t=df_otu_t,df_group=df_group,group_col=group_col))
}
adonis_wrap=function(phy_obj,method="bray",by=NULL){

    ## extract data from phyloseq object
    data_list=extract_phy(phy_obj)
    df_otu_t=data_list[[1]]
    df_group=data_list[[2]]
    group_col=data_list[[3]]
    
    ## build formula
    form=formula(sprintf("df_otu_t ~ %s",group_col))
    
    ## run
    res = adonis2(form, data = df_group, method = method,by=by)
    return(res)
}
adnois_pair=function(phy_obj,method="bray",p.adjust.m="BH"){

    ## extract data from phyloseq object
    data_list=extract_phy(phy_obj)
    df_otu_t=data_list[[1]]
    df_group=data_list[[2]]
    group_col=data_list[[3]]
    df_group=df_group[[group_col]]
    
    ## run
    res=pairwise.adonis(df_otu_t, df_group, p.adjust.m=p.adjust.m, sim.method=method)
    return(res)
}

"""
my_adonis=STAP(string,"my_adonis")

def plot_beta_diverstiy_common(df,dfg,df_taxa,df_taxa_id,
                               rank,method="MDS",dis="bray",
                               sample_color="group",sample_label="sample_name",
                               ellipse=True,clr_shift=0):
    # ---------------------------------------- description ---------------------------------------------
    """
    Plot the heatmap and the necessary diagnosis plots.
    
    Arguments:
    1. df: df_otu
    2. df_taxa: df_taxa
    3. rank: The rank you want to plot.
    4. method: The method you want to use for dimension reduction.
    5. dis: The distance method you use. Specifically, you can set dis as 'clr', then it will use the clr + euclidean distance.
    6. sample_color: The col names to define the group of each sample in dfg.
    7. sample_label: The col names to define the name of each sample in dfg. If you use build_phy to build the phyloseq object, the col name will be 'sample_name'
    8. clr_shift: The default shift in clr transform of microbiome package is small,
        so the clr results will not be stable when some samples are at abnormally low abundance.
        In this case you can manually set a larger clr_shift before the clr transformation, to 
        ignore the abnormal samples temporarily. 【The best way is not to use the shift, just drop the abnormal samples.】
    
    Outputs:
    1. fig: Rel taxa ordinate

    """
    # ---------------------------------------------- prepare --------------------------------------------
    
    # convert rank
    dfa,_,dft,_=convert_rank(df,df_taxa,df_taxa_id, rank)
    dfr=dfa/dfa.sum()
    
    # build abs and rel
    if dis=="clr":
        pretty_print("""
        Use the clr + PCA mode:
        1. Plotting: ignore the method argument and use 'PCoA' instead.
        2. Adonis: ignore the dis argument and use 'euclidean' instead.
        """)
        my_phy=build_phy(dfa,dfg,dft)
        my_phy=microbiome.transform(my_phy,transform="shift",shift=clr_shift)
        my_phy=microbiome.transform(my_phy,transform="clr")
        real_dis="euclidean"
        method="PCoA"# it is because clr + PCA = clr + euclidean + PCoA (PCoA=MDS)
    else:
        my_phy=build_phy(dfr,dfg,dft)
        real_dis=dis
    
    # get ordinate results
    ord_res=phyloseq.ordinate(my_phy,method,real_dis)
    
    # ---------------------------------------------- plot ----------------------------------------------
    
    ## plot main
    fig=phyloseq.plot_ordination(my_phy,ord_res,type="sample",color=sample_color)
    fig=fig+ggrepel.geom_text_repel(gp.aes_string(label=sample_label),size=5,show_legend=False)
    fig=mytheme(fig,"theme_bw")
    
    ## plot labels
    if dis=="clr":
        fig=fig+gp.labs(title="clr + PCA, rank = {}".format(rank))
    else:
        fig=fig+gp.labs(title="{} + {}, rank = {}".format(dis,method, rank))
    
    ## plot ellipse
    if ellipse:
        stat_ellipse=r["stat_ellipse"]
        fig=fig+stat_ellipse(level=0.95)
    
    ## adonis general table
    df_adonis=r2p(my_adonis.adonis_wrap(my_phy,by=NULL,method=real_dis))
    df_adonis.index.name=rank
    
    ## adnois pair-wais table
    df_adonis_pair=r2p(my_adonis.adnois_pair(my_phy,method=real_dis))
    df_adonis_pair.index.name=rank
    
    return(fig,df_adonis,df_adonis_pair)

phyloseq=importr("phyloseq")
microbiome=importr("microbiome")
pairwiseAdonis=importr("pairwiseAdonis")