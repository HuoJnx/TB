## rpy2

#### rpy2 util
%load_ext rpy2.ipython
import rpy2.robjects as ro

#### dataframe conversion
from rpy2.robjects import Formula
from rpy2.robjects import NA_Logical as NA
from rpy2.rinterface import NULL as NULL
from rpy2.robjects import r, rl
from rpy2.robjects.packages import importr

from rpy2.robjects.packages import STAP


def get_rlist(r_list, string, prompt=False):
    index=r_list.names.index(string)
    if prompt:
        print("The index of the results is {}".format(index))
    res=r_list[r_list.names.index(string)]
    if type(res)==type(ro.vectors.DataFrame(DF())):
        print("'{}' is data.frame, convert to pandas dataframe.".format(string))
        res=r2p(res)
    return(res)



#### dataframe conversion
from rpy2.robjects import pandas2ri

pandas2ri.activate()
from rpy2.robjects.pandas2ri import py2rpy as p2r
from rpy2.robjects.conversion import rpy2py as r2p

#### some r function
rprint = ro.globalenv.find("print")
from rpy2.robjects.packages import data as rdata

#### ggplot related
import rpy2.robjects.lib.ggplot2 as gp
from rpy2.ipython.ggplot import image_png


ggpubr = importr("ggpubr")
ggrepel = importr("ggrepel")
rstatix = importr("rstatix")
stringr = importr("stringr")
gridExtra = importr("gridExtra")
ggsave = r["ggsave"]
get_legend = r["get_legend"]


def ggsave_wrap_test(fig, path="test.png", size=[5, 5]):
    ggsave(path, fig, width=size[0], height=size[1])


def ggsave_wrap_pro(fig, fig_dir, fig_name, fig_fmt, size=[5, 5],prompt=False):
    os.makedirs(fig_dir, exist_ok=True)
    fig_path = os.path.join(fig_dir, fig_name) + "." + fig_fmt
    ggsave(fig_path, fig, width=size[0], height=size[1])
    if prompt:
        print("Successfully save figure to {}.".format(fig_path))


def plot(fig, size=[400, 400]):
    return image_png(fig, width=size[0], height=size[1])


def mytheme(
    fig,
    base="theme_bw",
    only_font=False,
    axis_size=15,
    axis_title_size=25,
    legend_title_size=25,
    legend_text_size=15,
    x_angle=0,
):

    # ---------------------------- set the basic theme ---------------------------

    if not only_font:
        print("Set basic theme to {}, them set font.".format(base))
        basic_list = [
            "theme_bw",
            "theme_classic",
            "theme_dark",
            "theme_grey",
            "theme_light",
            "theme_linedraw",
            "theme_minimal",
            "theme_void",
        ]
        if base in basic_list:
            base_theme = r[base]
            fig = fig + base_theme()
        else:
            raise Exception(
                "Only support the follow basic theme: {}".format(", ".basic_list)
            )
    else:
        print("Skip setting the basic theme, only set the font")

    # ----------------------------- set the font of axis and legend ----------------------------

    fig = fig + gp.theme(
        ## general
        axis_text=gp.element_text(size=axis_size),
        axis_title=gp.element_text(size=axis_title_size),
        legend_title=gp.element_text(size=legend_title_size),
        legend_text=gp.element_text(size=legend_text_size),
        ## specific
        axis_text_x=gp.element_text(angle=x_angle),
    )

    return fig


def my_box(fig):
    fig = fig + gp.geom_boxplot(outlier_shape=NA) + gp.geom_jitter()
    return fig


def extract_legend(fig, drop=False):
    l = get_legend(fig)
    if drop:
        fig = fig + gp.theme(legend_position="none")
    return (fig, l)


    
# plot statistics

string = """
plot_hypothesis_sub = function(fig, data, group_col, value_col,
    stat_m = "kruskal", pair_m = "wilcox", adj_m = "BH") {

    # ------------------------------------ conduct hypothesis test ------------------------------------
    if (stat_m == "kruskal") {
        form = formula(paste(value_col, "~", group_col))
        whole_res = data %>%
            kruskal_test(form)
    }

    if (is.na(pair_m)) {
        pair_res = NA
    } else if (pair_m == "wilcox") {
        pair_res = data %>%
            wilcox_test(form, p.adjust.method = adj_m)
    } else if (pair_m == "dunn") {
        pair_res = data %>%
            dunn_test(form, p.adjust.method = adj_m)
    } else {
        stop("Only support 'wilcox', 'dunn' or NA.")
    }

    # ------------------------------------ get info for plot ------------------------------------
    if (!is.na(pair_m)) {
        pair_res = pair_res %>%
            add_xy_position(x = group_col)
    }


    # ------------------------------------ plot ------------------------------------
    fig = fig + labs(subtitle = get_test_label(whole_res, detailed = T))
    if (!is.na(pair_m)){
        fig = fig + labs(subtitle = get_test_label(whole_res, detailed = T))
        fig = fig + stat_pvalue_manual(pair_res, hide.ns = T)
    }
    
    pair_res=pair_res%>%dplyr::select(all_of(c("group1","group2","statistic","p","p.adj","p.adj.signif")))
    return(list(general=whole_res,pair=pair_res,fig=fig))
}

try_get=function(fig){
    #print(class(fig))
    print(fig$mapping)
}
"""
my_stat = STAP(string, "_")
def plot_hypothesis(fig, data, group_col, value_col, stat_m = "kruskal", pair_m = "wilcox", adj_m = "BH"):
    # ------------------------------------ description -----------------------------------
    """
    Add general and pair_wise hypothesis test for plot.
    
    Arguments:
    1. fig: The fig you want to add info.
    2. data: The dataframe you plot the fig.
    3. group_col: The col in the data infering the group info for samples.
    4. value_col: The col in the data infering the value you used for general and pair_wise comparesion,.
    5. stat_m: The statistics method for doing the general comparesion, now accept 'kruskal', will support 'anova' in the future.
    6. pair_m: THe statistics method for doing the pair_wise comparesion, now only accept 'wilcox' or 'dunn'.
    7. adj_m: The adjust method for p-value results adjustment.
    
    Outputs:
    A fig with hypothesis test results.
    """
    # ------------------------------------ prompt ------------------------------------
    pretty_print("The used group col is '{}', value col is '{}'.".format(group_col,value_col))
    pretty_print("Method for general test is '{}'.".format(stat_m))
    if type(pair_m)==type(NA):
        pretty_print("Don't need pair_wise comparesion.")
    else:
        pretty_print("Method for pair_wise comparesion is '{}', that for p-value adjustment is '{}'.".format(pair_m, adj_m))
        
    # ------------------------------------- plot --------------------------------------
    fig = my_stat.plot_hypothesis_sub(fig, data, group_col, value_col, stat_m = stat_m, pair_m = pair_m, adj_m = adj_m)
    return(fig)


