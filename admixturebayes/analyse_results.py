import pandas as pd
import os

def generate_summary_csv(list_of_summaries, filename= 'summaries.csv', reference_tree=None):
    if reference_tree is not None:
        df=pd.DataFrame.from_items([('summary', [summ.name for summ in list_of_summaries])]+
                                   [('output', [summ.output for summ in list_of_summaries])]+
                                   [('value', [summ.summary_of_phylogeny(reference_tree) for summ in list_of_summaries])])
    else:    
        df=pd.DataFrame.from_items([('summary', [summ.name for summ in list_of_summaries])]+
                                   [('output', [summ.output for summ in list_of_summaries])])
    df.to_csv(filename)

def save_to_csv(list_of_tuples, summaries, filename='results.csv', origin_layer=(1,1)):
    df=pd.DataFrame.from_items([('iteration',list_of_tuples[0])]+[(summ_object.name,summ_col) for summ_object,summ_col in zip(summaries,list_of_tuples[1:])])
    if origin_layer is not None: 
        df['origin']=origin_layer[0]
        df['layer']=origin_layer[1]
    df.to_csv(filename)
    
def get_permut_filename(filename):
    parts=filename.split('.')
    return '.'.join(parts[:-1])+'permuts.'+parts[-1]
    
def save_permuts_to_csv(list_of_permuts, filename='results-permuts.csv'):
    with open(filename, 'w') as f:
        f.write('flip_number,'+','.join(map(str,list(range(len(list_of_permuts[0])))))+'\n')
        for n,permut in enumerate(list_of_permuts):
            f.write(str(n)+','+','.join(map(str,permut))+'\n')
    
def save_pandas_dataframe_to_csv(pd_dataframe, filename='results.csv'):
    pd_dataframe.to_csv(filename)