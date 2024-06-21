
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

# This collection was created for data analysis regarding mass spectrometry and the functions described are to be used along with 
# with the MassTRIX data sheet <<not comparing isotopes>>. 

# Most of the functions' syntaxes follow the same principle: function(KEGG_cids_list, File_1, File_2), in which KEGG_cids_list is a list with
# all the KEGG_cids from intended compound, perhaps obtained from KEGG_mapper tool for metabolism comparison, File_1 is 
# the name of the file refering to the first organism to be compared (usually the reference one) and File_2 is the name of 
# the second file refering to the second organism to be compared

# For specific signal attribution, make sure the data sheet obtained from MassTRIX has at least 1 ppm or less since two peaks
# may be attributed to the same molecule with the same ionization pattern because of a small observed shift 

# ----------------------------------------------------------------------------------------------------

def plot_cids(x, y, z): #Creates a plot for a list of KEGG_cids showing comparative values of normalized peak height for metabolites in both BY and KO types

# x = list with KEGG cids; y = type 1 file; z = type 2 file
# ---------------------------------------------
    Type1 = pd.read_csv(str(y), sep = '\t') # Reads file 1
    Type2 = pd.read_csv(str(z), sep = '\t') # Reads file 2
    T1allmet_one = Type1.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid') # File 1 df without duplicated KEGG cids
    T1allmet = Type1.set_index('KEGG_cid') # Sets KEGG cids as the df index for file 1
    T2allmet_one  = Type2.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid') # File 2 df without duplicated KEGG cids
    T2allmet = Type2.set_index('KEGG_cid') # Sets KEGG cids as the df index for file 2
    cid_listT1 = [cidT1 for cidT1 in Type1.KEGG_cid]# Creates a list with all identified KEGG cids in type 1

    cid_listT2 = [cidT2 for cidT2 in Type2.KEGG_cid]# Creates a list with all identified KEGG cids in type 2
# ---------------------------------------------

    ord = [] # List with normalized signals. Each metabolite is analyzed in both files and their signals are added in a list, which in turn id added to Ord.
    # If more than 1 metabolite is requested, their signals for both files will appear in tandem in Ord list
    for n in x:
        if n in cid_listT1 and n in cid_listT2: 
            val = [] # List with normalized signals for a coumpound in both files
            #-----------------------------------------------------------------------------------
            # Signal calculation for both files type
            # Check next function for description
            heighta = T2allmet.peak_height[str(n)].sum() 
            heightb = T1allmet.peak_height[str(n)].sum()
            enkBY = T1allmet.peak_height.HMDB01045.sum()
            enkKO = T2allmet.peak_height.HMDB01045.sum()
            qnt_a = heighta/enkKO
            qnt_b = heightb/enkBY
            val.append(qnt_a)
            val.append(qnt_b)
            #-----------------------------------------------------------------------------------
            for c in val: 
                ord.append(c)
        else: # Condition to specify cids that are not present in both files at the same time (or at neither of them)
            if n in cid_listT1:
                print(f'{n} appears only in {y} type')  
            elif n in cid_listT2:
                print(f'{n} appears only in {z} type')
            else:
                print(f'{n} does not appear as a signal')
    abc = [] # List with the name of a selected compound with 'T1' and 'T2' added to it.
    for n in x:
        if n in cid_listT1 and n in cid_listT2:
            name = [] #
            a = T2allmet_one.KEGG_name[str(n)].split(';')[0] # Extract the name from file 1
            b = T1allmet_one.KEGG_name[str(n)].split(';')[0] # Extract the name from file 2
            if a[-1] == ')':
                a = a[:-9] + '\nT2' # Excludes ionization pattern and add 'T2' to the name
                b = b[:-9] + '\nT1' # Excludes ionization pattern and add 'T1' to the name
            else:
                a = a + '\nT2' # Add 'T2' to the name
                b = b + '\nT1' # Add 'T1' to the name    
            name.append(a)
            name.append(b)
            for c in name:
                abc.append(c)
        else:
            continue
    with sns.axes_style('whitegrid'):
        plt.subplots(figsize=(6,6))
        c = sns.barplot(ord,abc)
    return c

def df_cids(x, y, z): #Creates a DataFrame for metabolites in both BY and KO types, comparing the normalized peak values for each metabolite in a list

# x = list with KEGG cids; y = type 1 file; z = type 2 file
# ---------------------------------------------
    Type1 = pd.read_csv(str(y), sep = '\t') # Reads file 1
    Type2 = pd.read_csv(str(z), sep = '\t') # Reads file 2
    T1allmet_one = Type1.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid') # File 1 df without duplicated KEGG cids
    T1allmet = Type1.set_index('KEGG_cid') # Sets KEGG cids as the df index for file 1
    T2allmet_one  = Type2.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid') # File 2 df without duplicated KEGG cids
    T2allmet = Type2.set_index('KEGG_cid') # Sets KEGG cids as the df index for file 2
    
    cid_listT1 = [cidT1 for cidT1 in Type1.KEGG_cid]# Creates a list with all identified KEGG cids in type 1

    cid_listT2 = [cidT2 for cidT2 in Type2.KEGG_cid]# Creates a list with all identified KEGG cids in type 2
# ---------------------------------------------

    ordT2 = [] # List with normalized values for type 2 file
    ordT1 = [] # List with normalized values for type 1 file
    for n in x:
        if n in cid_listT1 and n in cid_listT2:
            val = []
            heighta = T2allmet.peak_height[str(n)].sum() # Calculates peak intesity for the used standard in T2
            heightb = T1allmet.peak_height[str(n)].sum() # Calculates peak intesity for the used standard in T1
            enkT1 = T1allmet.peak_height.HMDB01045.sum() # Calculates peak intesity for the used standard in T1
            enkT2 = T2allmet.peak_height.HMDB01045.sum() # Calculates peak intesity for the used standard in T2
            qnt_a = heighta/enkT2 # Calculates the normalized signal for x in T2
            qnt_b = heightb/enkT1 # Calculates the normalized signal for x in T1
            val.append(qnt_a)
            val.append(qnt_b)
            ordT2.append(val[0])
            ordT1.append(val[1])
        else:
            if n in cid_listT1:
                print(f'{n} appears only in {y} type')  
            elif n in cid_listT2:
                print(f'{n} appears only in {z} type')
            else:
                print(f'{n} does not appear as a signal')
    abcKO = [] # List with compound names for type 2 file
    abcBY = [] # List with compound names for type 1 file
    for n in x:
        if n in cid_listT1 and n in cid_listT2:
            name = []
            a = T2allmet_one.KEGG_name[str(n)].split(';')[0] # Extract the name of the selected compound based on its KEGG cid for type 2 file
            b = T1allmet_one.KEGG_name[str(n)].split(';')[0] # Extract the name of the selected compound based on its KEGG cid for type 1 file
            if a[-1] == ')': # Condition to remove the ionization ion from the name if it appears in it <<Name '[Formula + Ion]' >>, 
                a = a[:-9] + '\nT2'
                b = b[:-9] + '\nT1'
            else:
                a = a + '\nT2'
                b = b + '\nT1'    
            name.append(a)
            name.append(b)
            abcKO.append(name[0][:-3])
            abcBY.append(name[1][:-3])                
        else:
            continue  

    dictKO = dict(zip(abcKO, ordT2)) # Creates a dictionary for type 2 file in which the key is the compound name and the value its normalized signal
    dictBY = dict(zip(abcBY, ordT1)) # Creates a dictionary for type 1 file in which the key is the compound name and the value its normalized signal
    dict_list = [dictKO,dictBY] # Creates a list with the two dictionaries for transformation into a DataFrame
    
    c = pd.DataFrame(dict_list, index=['T2','T1']).transpose()
    return c

def quantify(x, y, z): # Calculate the normalized value for an specific metabolite based on its KEGG_cid

# x = KEGG cid; y = type 1 file; z = type 2 file
# ---------------------------------------------
    Type1 = pd.read_csv(str(y), sep = '\t') # Reads file 1
    Type2 = pd.read_csv(str(z), sep = '\t') # Reads file 2
    T1allmet_one = Type1.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid') # File 1 df without duplicated KEGG cids
    T1allmet = Type1.set_index('KEGG_cid') # Sets KEGG cids as the df index for file 1
    T2allmet_one  = Type2.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid') # File 2 df without duplicated KEGG cids
    T2allmet = Type2.set_index('KEGG_cid') # Sets KEGG cids as the df index for file 2
    cid_listT1 = [] # Creates a list with all identified KEGG cids in type 1
    cid_listT1 = [cidT1 for cidT1 in Type1.KEGG_cid]# Creates a list with all identified KEGG cids in type 1

    cid_listT2 = [cidT2 for cidT2 in Type2.KEGG_cid]# Creates a list with all identified KEGG cids in type 2
# ---------------------------------------------

    enkT1 = T1allmet.peak_height.HMDB01045.sum() # Calculates peak intesity for the used standard in T1
    enkT2 = T2allmet.peak_height.HMDB01045.sum() # Calculates peak intesity for the used standard in T2
    if x in cid_listT1 and x in cid_listT2:
        heighta = T2allmet.peak_height[str(x)].sum() # Calculates peak intesity for x in T1
        heightb = T1allmet.peak_height[str(x)].sum() # Calculates peak intesity for x in T2
        qnt_a = heighta/enkT2 # Calculates the normalized signal for x in T2
        qnt_b = heightb/enkT1 # Calculates the normalized signal for x in T1
        name = T2allmet_one.KEGG_name[str(x)].split(';')[0] # Takes the metabolite name
        if name[-1] == ')': # Removes the ion indication if it is present in the name
            name = name[:-9] 
        else:
            name = name
        return print(f'The normalization for {name} is {qnt_b} for T1 and {qnt_a} for T2')
    else:
        if x in cid_listT1: # Calculates normalization if the coumpound is present only in T1
            b = T1allmet.peak_height[str(x)].sum()
            qnt_b = b/enkT1
            name = T1allmet_one.KEGG_name[str(x)].split(';')[0]
            return print(f'The normalization for {name} is {qnt_b} in T1 only')
        elif x in cid_listT2: # Calculates normalization if the coumpound is present only in T2
            a = T2allmet.peak_height[str(x)].sum()
            qnt_a = a/enkT2
            name = T2allmet_one.KEGG_name[str(x)].split(';')[0]
            return print(f'The normalization for {name} is {qnt_a} in T2 only') 

def cid_color(x, y, z): #Creates a DataFrame with every KEGG cid and associates a color based on the molecule's presence in only one or both of the files analyzed
    
# x = type 1 file; y = type 2 file; z = file to send <<KEGG_cid | color>> pairs
# ---------------------------------------------
    Type1 = pd.read_csv(str(x), sep = '\t') # Reads file 1
    Type2 = pd.read_csv(str(y), sep = '\t') # Reads file 2
    T1allmet_one = Type1.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid') # File 1 df without duplicated KEGG cids
    T1allmet = Type1.set_index('KEGG_cid') # Sets KEGG cids as the df index for file 1
    T2allmet_one  = Type2.drop_duplicates(subset='KEGG_cid',keep='first').set_index('KEGG_cid') # File 2 df without duplicated KEGG cids
    T2allmet = Type2.set_index('KEGG_cid') # Sets KEGG cids as the df index for file 2
    cid_listT1 = [] # Creates a list with all identified KEGG cids in type 1
    cid_listT1 = [cidT1 for cidT1 in Type1.KEGG_cid]# Creates a list with all identified KEGG cids in type 1

    cid_listT2 = [cidT2 for cidT2 in Type2.KEGG_cid]# Creates a list with all identified KEGG cids in type 2
# ---------------------------------------------

    cid_listALL = [] # List with all KEGG cids in both files (not repeated)
    for n in cid_listT2:
        if n in cid_listT1:
            continue
        else:
            if n.startswith('C'):
                cid_listALL.append(n) # Add KEGG cids in T2 if they are not present in T1


    for n in cid_listT1:
        if n.startswith('C'):
            cid_listALL.append(n) # Add all T1 KEGG cids


    colors = [] # List with associated color for each KEGG cid in cid_listALL
    for cid in cid_listALL:
        if cid in cid_listT1 and cid in cid_listT2:
             colors.append('yellow') # If cid is present in both files, the attributed color is yellow
        else:
            if cid in cid_listT1:
                colors.append('blue') # If cid is present only in T1, the attributed color is blue
            elif cid in cid_listT2:
                colors.append('red') # If cid is present only in T1, the attributed color is red
            else:
                continue
    df = pd.DataFrame(colors,index = cid_listALL).reset_index() # both lists are joined in a DataFrame
    return df.to_csv(str(z), header=None, index=None, sep=' ', mode='a') # Sends the KEGG cids and colors to a csv (z) file that 
                                                                         # can be used in KEGG mapper tool