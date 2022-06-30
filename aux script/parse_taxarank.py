## description
"""
For parsing the results from taxaranks, but only for that from running with argument '-t' like the bellow:

taxaranks -i Whole_taxa_id.csv -o Whole_taxa_lineage -t

""" 
## prepare

#### normal import
import pandas as pd
import sys,os

#### self-function
def pretty_print(t):
    print(t + "\n")

#### get argv

infile=sys.argv[1]

#### get file_name and dir_path for input
file_name=os.path.basename(infile)
file_name_without_suffix=file_name.split(".")[0]
dir_path=os.path.dirname(infile)

#### set output name
file_name_taxaname=file_name_without_suffix+"_taxa_name"+".tsv"
file_name_taxaid=file_name_without_suffix+"_taxa_id"+".tsv"

#### prompt
pretty_print("""
The input file is {},
output the taxa_name as {},
taxa_id as {}.""".format(infile,file_name_taxaname,file_name_taxaid))





## read df (The low_memory=False is for muting the mixed types waringing.)
df=pd.read_csv(infile,sep="\t",low_memory=False)

## drop the 'taxa_searched' column, and set the 'user_taxa' column as '#OTU ID'(biom format)
df=df.drop(columns="taxa_searched")
df=df.set_index("user_taxa")
df.index.name="#OTU ID"

## seperate to 2 parts


#### df_name

###### extract the columns not containing "taxid"
df_name=df.loc[:,~df.columns.str.contains("taxid")]
###### save
df_name.to_csv(file_name_taxaname,sep="\t")

#### df_id

def toint(e):
    try:
        e=int(e)
    except:
        e=pd.NA
    return(e)

###### extract the columns contains taxid
df_taxid=df.loc[:,df.columns.str.contains("taxid")]
###### to int
df_taxid=df_taxid.applymap(toint)
###### drop the "taxid" in columns
df_taxid.columns=df_taxid.columns.str.removesuffix("_taxid")
###### save
df_taxid.to_csv(file_name_taxaid,sep="\t")

#### prompt

pretty_print("Finished.")