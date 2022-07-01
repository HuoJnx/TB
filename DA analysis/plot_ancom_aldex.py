def plot_ancom_aldex(
    df_final,
    col_to_color,
    col_for_color,
    sorting="none"
    x_y=["index","W"],
    xy_lab=["Taxa","W-statistics"],
):
    
    # plot config
    ## some version of gp.scale_colour_manual will raise error if not drop duplicates
    color_map=df_final[[col_to_color,col_for_color]].drop_duplicates()
    
    # ---------------------------plot----------------------------
    
    fig = (
        ggpubr.ggdotchart(
            df_final,
            x=x_y[0],
            y=x_y[1],
            color=col_to_color,
            sorting=sorting, #default is "descending". However, something bug in here, now 'none' equalled to 'ascending'
            add="segments",
            add_params=r.list(color=col_to_color, size=1),
            rotate=1,  # Rotate vertically
            y_text_col=0,  # Whether to coloring text in y. However, something bug in here, now it won't coloring text no matter you set what.
            ggtheme=ggpubr.theme_pubr(),  # ggplot2 theme
            shape=17,
            xlab=xy_lab[0],
            ylab=xy_lab[1],
        )
        + ggpubr.theme_cleveland()
        + gp.scale_colour_manual(breaks=color_map[col_to_color].values, values=color_map[col_for_color].values)
        )
    return fig
