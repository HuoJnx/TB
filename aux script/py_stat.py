def description(data,group_col,value_col,round=3):
    ## select group, and select value
    df_group_select=data.groupby(group_col)[value_col]
    ## get basic description
    base_describe=df_group_select.describe()
    ## get sem
    sem=df_group_select.sem()
    sem.name="sem"
    ## merge
    final_describe=base_describe.merge(sem,left_index=True,right_index=True)
    ## round
    final_describe=final_describe.round(round)
    return(final_describe)

def df_save_pro(df,df_dir,df_name,df_fmt,header=True,index=True,prompt=True):
    os.makedirs(df_dir, exist_ok=True)
    df_path=os.path.join(df_dir,df_name) + "." + df_fmt
    if df_fmt=="tsv":
        df.to_csv(df_path,sep="\t",header=header,index=index)
    elif df_fmt=="csv":
        df.to_csv(df_path,sep=",",header=header,index=index)
    if prompt:
        print("Successfully save dataframe to {}.".format(df_path))