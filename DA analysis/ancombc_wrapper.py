def my_ancombc(
    phy_obj,
    p_adj_method="BH",
    sort="W",
    color="compare",
    ascending=True,
    only_sig=True,
    alpha=0.05,
    conserve=False,
    xy_lab=["Taxa", "W-statistics"],
    highlight="<",
    plot_limit=50,
):
    ### -----------------------------explaination-------------------------------
    """
    Notice: 
    
    1. As default, the W <0 means the first condition has higher abundance, and W>0 means the
        second condition has higher abundance.
    2. df are taxa in rows, samples in columns.
    
    Arguments:
    
    sort: use which col to sort, supported cols are index, beta, se, W, p_val, q_val, diff_abn, compare.
    color: use which col to define color, supported cols are index, beta, se, W, p_val, q_val, diff_abn,
            compare.
    only_sig: If True, drop non-sig rows
    alpha: significance threshold
    conserve: Whether to use the conserve mode
    xy_lab: x, y label in plotting
    highlight: The palette has two color, one is blue another is orange, and orange is defined as
                highlight color. So you can define the highlight keyword, if check the keyword, the
                conditin will be in orange.
    plot_limit: limit the rows to plot
    """
    ###
    # -----------------------------necessary-------------------------------
    ## R-wrapper
    string = """
    my_ancombc = function(phy_obj, p_adj_method = "holm",
        alpha = 0.05, conserve = T) {

        ## get the colname of the first column of meta table
        df_meta=sample_data(phy_obj)
        df_meta_group_col_name=names(df_meta)[1]
        print(sprintf("Use %s in meta table as factor.",df_meta_group_col_name))
        
        ## run ANCOMBC
        out = ancombc(phyloseq = phy_obj, formula = df_meta_group_col_name,
            p_adj_method = p_adj_method, alpha = alpha, conserve = conserve)
        
        ## get original df_otu, df_meta
        df_otu=as.data.frame(as.matrix(otu_table(phy_obj)))
        df_meta=as.data.frame(as.matrix(sample_data(phy_obj)))
        
        ## return
        return(list(res=out$res,df_otu=df_otu,df_meta=df_meta))
    }

    """

    ## run ancombc
    my_ancombc_sub = STAP(string, "_")
    if conserve:
        pretty_print("Use conserve mode.")
    else:
        pretty_print("Non-conserve mode.")
    res_list = my_ancombc_sub.my_ancombc(
        phy_obj, p_adj_method, alpha=alpha, conserve=conserve
    )
    
    ##parse results
    res=get_rlist(res_list, "res")
    df_otu=get_rlist(res_list, "df_otu")
    df_meta=get_rlist(res_list,"df_meta")
    
    ## merge
    df_list = []
    for k, v in res.items():
        v=r2p(v)
        original_name = v.columns[0]
        v.columns = [k]
        df_list.append(v)
    df_final = pd.concat(df_list, axis=1)
    

    # -----------------------------data polishment-relative-------------------

    ## filter the significant
    if only_sig:
        df_final = df_final.query("diff_abn==1")
        pretty_print("Keep significant, threshold is {}.".format(alpha))
    else:
        pretty_print("Keep all.")

    ## sort
    try:
        df_final = df_final.sort_values(by=sort, ascending=ascending)
        pretty_print("Sorted by {}.  ".format(sort))
    except:
        pretty_print("Not sorted.")

    # -----------------------------plotting-relative-------------------

    ## reset index for plot
    df_final = df_final.reset_index()
    df_final.index.name = original_name

    ## change the label of significant indicator
    df_final["diff_abn"] = df_final["diff_abn"].map({0: "False", 1: "True"})

    ## add the compare col
    fac_list = df_meta.iloc[:, 0].unique()
    val1 = "{} < {}".format(fac_list[0], fac_list[1])
    val2 = "{} > {}".format(fac_list[0], fac_list[1])
    df_final["compare"] = df_final["W"].apply(lambda x: val1 if x > 0 else val2)

    ## manually add the color for each value
    palette = ["#00AFBB", "#FC4E07"]# The first is blue and second is orange
    df_final["color"] = df_final[color].map(
        lambda x: palette[1] if highlight in x else palette[0]# Highlight in orange and the other in blue
    )
    color_map = df_final.reindex(columns=[color, "color"]).drop_duplicates()

    # -----------------------------plotting----------------------------

    ## limit the max df length
    if len(df_final) > plot_limit:
        pretty_print("Plot limited is {}, will cut some results.".format(plot_limit))
        df_plot = df_final.copy()
        df_plot["abs_W"]=df_plot["W"].abs()
        df_plot=df_plot.sort_values(by="abs_W",ascending=False)
        df_plot = df_plot.iloc[0:plot_limit, :]
    else:
        df_plot = df_final

    ## get the sort value for plot
    if ascending:
        sorting = "ascending"
    else:
        sorting = "descending"

    ## plot
    if len(df_plot)==0:
        pretty_print("No significant, so no plot.")
        return None,None,None
    fig = (
        ggpubr.ggdotchart(
            df_plot,
            x="index",
            y="W",
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

    # -----------------------------check df----------------------------

    df_check = df_otu.loc[df_final["index"], :]
    #df_check=0
    return df_final, df_check, fig