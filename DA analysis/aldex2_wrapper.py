def my_aldex(
    phy_obj,
    denom="iqlr",
    effect=True,
    use="welch",
    alpha=0.05,
    sort="effect",
    ascending=True,
    only_sig=True,
    plot_limit=50,
    color="compare",
    highlight="<",
    xy_lab=["Taxa","Effect"],
):
    # -----------------------------necessary-------------------------------
    string = """
    my_aldex = function(phy_obj, denom="iqlr", effect=T ) {
    
        ## get inputs from phy_obj for aldex2
        df_otu=as.data.frame(as.matrix(otu_table(phy_obj)))
        
        df_meta=as.data.frame(as.matrix(sample_data(phy_obj)))
        df_meta_group_col_name=names(df_meta)[1]
        df_group=df_meta[[1]]
        print(sprintf("Use %s in meta table as factor.",df_meta_group_col_name))
        
        ## run ALDEX2
        res=aldex(df_otu,df_group,denom=denom,effect=effect)
        
        ## return
        return(list(res=res,df_otu=df_otu,df_meta=df_meta))
    }

    """

    ## run aldex
    my_aldex_sub = STAP(string, "_")
    res_list = my_aldex_sub.my_aldex(phy_obj, denom=denom, effect=effect)
    
    ##parse results
    res=get_rlist(res_list, "res")
    df_otu=get_rlist(res_list, "df_otu")
    df_meta=get_rlist(res_list,"df_meta")

    # -----------------------------data polishment-relative-------------------

    ## drop some redundant
    df_final = res.drop(columns=["rab.all", "overlap"])

    ## rename the necessary cols
    if use == "wilcox":
        print("Use wilcox, drop Welch results")
        df_final = df_final.rename(columns={"wi.ep": "p_val", "wi.eBH": "q_val"})
    elif use == "welch":
        print("Use Welch, drop Wilcox results")
        df_final = df_final.rename(columns={"we.ep": "p_val", "we.eBH": "q_val"})

    ## add diff_abn
    df_final["diff_abn"] = df_final["q_val"].apply(
        lambda x: "True" if x < alpha else "False"
    )

    ## filter non-sig?
    if only_sig:
        df_final = df_final.query("diff_abn == 'True'")

    ## calculate the compare col
    fac_list = df_meta.iloc[:, 0].unique()
    val1 = "{} < {}".format(fac_list[0], fac_list[1])
    val2 = "{} > {}".format(fac_list[0], fac_list[1])
    df_final["compare"] = df_final["effect"].apply(lambda x: val1 if x > 0 else val2)

    ## sort by the col you want
    df_final = df_final.sort_values(by=sort, ascending=ascending)
    
    ## change the order
    first_batch = ["p_val", "q_val", "diff_abn", "effect", "compare"]
    other_batch = df_final.columns.difference(first_batch).to_list()
    first_batch.extend(other_batch)
    new_order = first_batch
    df_final = df_final.reindex(columns=new_order)
    
    ## reset index
    df_final=df_final.reset_index()
    
    # --------------------------------plot----------------------------------------
    
    ## limit the length of df for plot
    if len(df_final) > plot_limit:
        print("Plot limited is {}, will cut some results.".format(plot_limit))
        df_plot = df_final.copy()
        df_plot["abs_W"]=df_plot["W"].abs()
        df_plot=df_plot.sort_values(by="abs_W",ascending=False)
        df_plot = df_plot.iloc[0:plot_limit, :]
    else:
        df_plot = df_final
        
    ## manually add the color for each value
    palette = ["#00AFBB", "#FC4E07"]# The first is blue and second is orange
    df_final["color"] = df_final[color].map(
        lambda x: palette[1] if highlight in x else palette[0]# Highlight in orange and the other in blue
    )
    
    ## some version of gp.scale_colour_manual will raise error if not drop duplicates
    color_map = df_final.reindex(columns=[color, "color"]).drop_duplicates()
    
    ## get the sort value for plot
    if ascending:
        sorting = "ascending"
    else:
        sorting = "descending"
        
    ## plot
    if len(df_plot)==0:
        print("No significant, so no plot.")
        return(None,None,None)
    fig = (
        ggpubr.ggdotchart(
            df_plot,
            x="index",
            y="effect",
            color=color,
            add="segments",
            add_params=r.list(color=color, size=1),
            rotate=1,  # Rotate vertically
            y_text_col=1,  # Color y text by groups
            ggtheme=ggpubr.theme_pubr(),  # ggplot2 theme
            shape=17,
            xlab=xy_lab[0],
            ylab=xy_lab[1],
            sorting=sorting,
        )
        + ggpubr.theme_cleveland()
        + gp.scale_color_manual(
            breaks=color_map[color].values, values=color_map["color"].values
        )
    )
    
    # --------------------------------df_check----------------------------------------
    
    df_check = df_otu.loc[df_final["index"], dfg.index]
    
    return df_final, df_check, fig