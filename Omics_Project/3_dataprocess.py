import pandas as pd

def sum_signal_matrix(data:list, df:pd.DataFrame, samples:list, label:str): # function for concatenating the average value of gene signal for a given probeset based on gene symbol identifier
    count = 0
    avg_df = pd.DataFrame()
    for gene in data:
        table = df.loc[df[label] == gene]
        av_sig = dict()
        for sample in samples:
            mean_val = table[sample].median()
            av_sig[sample] = mean_val
        table = pd.DataFrame(av_sig,index=[1])
        table[label] = gene
        count += 1
        print(f'{count}: parsing information for gene {gene}')
        avg_df = pd.concat([avg_df,table],ignore_index=True)
    return avg_df

def group_by_genes(df,data:list,label:str): # funcition for grouping df rows by gene symbol
    fc_avarage = pd.DataFrame()
    for gene in data:
        table = df[df[label] == gene]
        fc_avarage = pd.concat([fc_avarage,table],ignore_index=True)
    return fc_avarage

def median_signal_matrix(data:list, df:pd.DataFrame, samples:list, label:str): # function for concatenating the average value of gene signal for a given probeset based on gene symbol identifier
    count = 0
    avg_df = pd.DataFrame()
    for gene in data:
        table = df.loc[df[label] == gene].median()
        count += 1
        print(f'{count}: parsing information for gene {gene}')
        avg_df = pd.concat([avg_df,table],ignore_index=True)
    return avg_df

def median_collapse(df:pd.DataFrame, # a DEG data frame with 3 columns: **logFC, adj.P.Val and Gene Symbol**
                    glist:list, # list with unique gene sumbols identified by the affymemtrix analysis
                    label:str,# column of the dataframe that will be used to group and filter
                    verbose:bool): 
    empty = pd.DataFrame()
    for gene in glist:
        table = df[df[label] == gene]
        #print(test)
        median = table.median(numeric_only=True) # calculate the median value of every variable
        #print(median)
        if median['logFC'] in table['logFC'].tolist(): # in case the probeset contains an even amount of probes...
            table = table.loc[table['logFC']==median['logFC']] # select row that contains the median logFC value
            table['Probeset Type'] = 'odd' # add a control column
            empty = pd.concat([empty,table],ignore_index=True) # add this row to an external empty df to save it
            if verbose == True:
                print(f'parsing information for gene {gene}, odd probeset')
        else:# in case the probeset contains an odd amount of probes...
            dic = dict()
            dic['logFC'] = float(median['logFC']) # add the median logFC value to save it
            p_vals = dict()
            for p_val in table['adj.P.Val'].tolist(): 
                diff = p_val - median['adj.P.Val'] # calculate the difference between each p_val and the median calculated one
                p_vals[str(p_val)] = diff # add to ann empty dict the difference as a value and the original p_val as a key (in string)
            diff_list = list()
            for key, value in p_vals.items(): 
                if value > 0:
                    diff_list.append(key) # if the calculate difference is bigger than zero, add the respective key (the original p_val) to a separate list
            for keys in diff_list:
                p_vals.pop(keys) # remove p_vals that greater than the mean
            max_diff = max(list(p_vals.values())) # find the closest p_val to the mean from the resultant ones
            #print(list(p_vals.values()))
            for key, value in p_vals.items():
                if value == max_diff:
                    p_val = key # fetch the desired orginal p_val that is closest to the median one
            dic['adj.P.Val'] = round(float(p_val),15) # add the desired p_val to the previous dict
            dic['Gene Symbol'] = gene # add the Gene Symbol column with the respective gene
            table = pd.DataFrame(dic, index = [0]) # transform dict into a df 
            table['Probeset Type'] = 'even' # add a control column
            empty = pd.concat([empty,table],ignore_index=True) # add this row to an external empty df to save it 
            if verbose == True:
                print(f'parsing information for gene {gene}, even probeset') 
            #print(p_val)
    return empty

def write(file,text): # write result matrix into another file for further processing 
    with open(file, 'w') as file:
        file.write(text)

