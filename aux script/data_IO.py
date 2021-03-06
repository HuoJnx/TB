#-------------------------- read --------------------------
otu_sep="\t"
meta_sep=","
taxa_sep="\t"

def read_otu(otu_path,sep=otu_sep):
    df_otu=pd.read_csv(otu_path,index_col=0,sep=sep)
    return(df_otu)

def read_meta(meta_path,sep=meta_sep):
    df_meta=pd.read_csv(meta_path,index_col=0,sep=sep)
    return(df_meta)

def read_taxa(taxa_path,sep=taxa_sep):
    df_taxa=pd.read_csv(taxa_path,index_col=0,sep=sep)
    return(df_taxa)

def read_taxa_id(taxa_path,sep=taxa_sep):

    def toint(e):
        try:
            e=int(e)
        except:
            e=pd.NA
        return(e)

    df_taxa_id=pd.read_csv(taxa_path,index_col=0,sep=sep)
    df_taxa_id=df_taxa_id.applymap(toint)

    return(df_taxa_id)

def read_phy_data(otu_path,meta_path,taxa_path):
    df_otu=read_otu(otu_path)
    df_meta=read_meta(meta_path)
    df_taxa=read_taxa(taxa_path)
    return(df_otu,df_meta,df_taxa)

def read_all_data(otu_path,meta_path,taxa_path,taxa_id_path):
    df_otu=read_otu(otu_path)
    df_meta=read_meta(meta_path)
    df_taxa=read_taxa(taxa_path)
    df_taxa_id=read_taxa_id(taxa_id_path)
    return(df_otu,df_meta,df_taxa,df_taxa_id)

#-------------------------- write --------------------------


def write_otu(df_otu,otu_path,sep=otu_sep):
    df_otu.write(otu_path,sep=otu_sep)

def write_meta(df_meta,meta_path,sep=meta_sep):
    df_meta.write(meta_path,sep=meta_sep)

def write_taxa(df_taxa,taxa_path,sep=taxa_sep):
    df_taxa.write(taxa_path,sep=taxa_sep)

#-------------------------- use --------------------------
#otu_path="otu_tables (new_version, only include the one with full information from order to species)/OTU_ABS.tsv"
#meta_path="group_info (maunal, drop DT8, HT2).csv"
#taxa_path="taxa_tables/Whole_taxa_lineage_taxa_name.tsv"
#taxa_id_path="taxa_tables/Whole_taxa_lineage_taxa_id.tsv"
#df_otu,df_meta,df_taxa=read_phy_data(otu_path,meta_path,taxa_path)
#df_otu,df_meta,df_taxa,df_taxa_id=read_all_data(otu_path,meta_path,taxa_path,taxa_id_path)