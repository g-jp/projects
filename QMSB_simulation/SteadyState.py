from stimator import read_model
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from collections import Counter
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import math as m
from tqdm import tqdm

def parameter_extraction(model:str):
    parameters = {}
    for line in model.split('\n'):
        if '=' in line and 'init' not in line:
            parm_value = line.split('=')
            parameter = parm_value[0].strip()
            value = float(parm_value[-1].strip())
            parameters[parameter] = value
    return parameters

def script(file,k,value):
    cut=file.split('\n')
    #cut[0] = cut[0] + str(k)
    for n in cut:
        if 'title' in n:
            cut[cut.index(n)] = n + ' ' + str(k) + '=' + str(value)
        if str(k) in n and ',' not in n:
            replace = n.split('=')
            replace[-1]=str(value)
            new = ' = '.join(replace)
            cut[cut.index(n)] = new
            final = '\n'.join(cut)
            return str(final)

def scan_prep(file,reaction,B):
    cut=file.split('\n')
    for n in cut:
        if 'init' in n:
            last = n
        if 'title' in n:
            cut[cut.index(n)] = n + str(reaction)
        if str(reaction) in n and ',' in n:
            add = n.split(', ')
            add[-1] = 'B*' + str(add[-1])
            add_new = ','.join(add)
            cut[cut.index(n)] = add_new
        if 'init' in n:
             cut[cut.index(n)] = ''
    #print(cut)
    final = cut.extend(['B = ' + str(B),'',last])
    final = '\n'.join(cut)
    return final

def sensibility_period(models):
    times = []
    for model in models:
        analysis = SteadyState(model,tf=20000,npoints=1000).response_rate_max(['NFkBnucTOT'])
        time = analysis['time NFkBnucTOT'].tolist()
        period_s = float(time[1]) - float(time[0])
        period_m = period_s/60
        times.append(period_m)
    
    sensibilities = []
    for time in times[1:]:
        sensibility = ((time - times[0])/times[0])/0.2
        sensibilities.append(sensibility)

    param_chang = []
    for model in models[1:]:
        split = model.split('\n')
        for line in split:
            if '(1' in line and 'title' not in line:
                parameter = line.split('=')[0].strip()
                param_chang.append(parameter)
    df = pd.DataFrame(zip(sensibilities,times[1:]),param_chang,columns=['Sensibility 20%','Period'])
    return df

def parameter_extraction(model:str):
        parameters = {}
        for line in model.split('\n'):
            if '=' in line and 'init' not in line and '~' not in line:
                parm_value = line.split('=')
                parameter = parm_value[0].strip()
                value = float(parm_value[-1].strip())
                parameters[parameter] = value
        return parameters

class SteadyState:
    all = []
    def __init__(self,model:str,change={},tf=1000,npoints=1000,outputs=[],scan=False): # A constructor that is used to insert instance variables (i.e., variables that are particular to the created object)
        for line in model.split('\n'):
            if "title" in line:
                title = line.strip('title').strip()
        print(f"Processing {title}...")
        assert type(model) == str, f'Model {model} is not a string!'
        self.simtime = tf
        self.resolution = npoints
        self.parmchange = change
        self.model = model
        self.outputs = outputs
        self.parameters = parameter_extraction(self.model)
        self.model_solve = read_model(model).solve(tf=tf,npoints=npoints,outputs=outputs) #Creates the "model" instance variable that is the solved model 
        if scan == True:
            self.scanparm = change
            self.model_scan = read_model(model).scan(change,tf=tf,npoints=npoints,outputs=outputs) #Creates the "model" instance variable that is the scanned model
        SteadyState.all.append(self)
        

    def __repr__(self):
        title = self.model.split('\n')
        emp_list = [element for element in title if len(element) > 1]
        value = emp_list[0]
        return value
    
    def steady_scan(self, lists, stimuli, sensibility=False):
        sa_values = {} # Empty dictionary for sensibility values 
        ss_values = {} # Empty dictionary for steady-state values
    #_________________________________________________________________________________________________________________________
    # The scan object is a result of applying the scan function to the model
    # The list required is a list containing a string with the variables to be analyzed
    # The stimuli list is a list contaning different percentage values of change for a paramter as to determine the sensibility of the system to that parameter
    # The sensibility term is an option so as to decide if the sensibility analysis will be performed
    #_________________________________________________________________________________________________________________________
        print("Sensibility analysis...")
        for variable in lists: # For a variable in a list containing all variables in strings
            for index in tqdm(range(0, len(self.model_scan))): # For a number assigined to all tested conditions
                ss=self.model_scan[index][variable].tolist()[-1] # Iterates through each test condition and each variable, evaluating it's Steady-State (SS) value
                ss_values.update({variable+' '+str(index):ss})  # Adds this value to a dictionary
        if sensibility == True: # Condition to peform a sensibility analysis
            for variable in lists:
                for index in range(0,len(self.model_scan)):
                    ssn = ss_values[variable+' '+str(index)] #Steady-state in a specific condition
                    ss0 = ss_values[variable+' '+'0'] # Steady-state in the first "normal" condition
                    sa = ((ssn-ss0)/ss0)/((stimuli[index]-stimuli[0])/stimuli[0]) if index != 0 else None # Calculated sensibility change due to changes in a parameter (defined by the percentage change in this parameter, B)
                    sa_values.update({variable+' '+str(index):sa})
            return pd.DataFrame([ss_values, sa_values],index=['SS','SA']).transpose() # Creates a DataFrame from the SS values and the SA values for each variable in each tested condition
        else:
            return ss_values

    def response_rate_max(self,variables):
        # Variables: list of strings containing the variables from the model to be analyzed.
        split = str(self.model_solve).split('\n') #Transform he model results into a list of strings, each element of the list corresponding to a line fo the results (the results are shown in a DataFrame)
        df = []
        for v in variables: #For v in a list containing the variables
            values_plus = []
            for i in range(0,len(self.model_solve[v])-2): 
                x = self.model_solve[v][i+1] - self.model_solve[v][i]
                y = self.model_solve[v][i+2] - self.model_solve[v][i+1]
                if x > 0 and y < 0:
                    values_plus.append(self.model_solve[v][i+1])
            link = {v:values_plus}
            data = pd.DataFrame(link)
            #print(data)
            temp_max=[]
            for line in str(self.model_solve).split('\n'):
                for max in data[v].tolist():
                    if str(max) in line:
                        time = line.split(' ')[0]
                        temp_max.append(time)
            dafr = pd.DataFrame(temp_max,data[v].tolist(),['time ' + str(v)])
            dafr.index.name = 'Max Value'
            df.append(dafr)
        return df

    def response_rate_min(self,variables):
        # Variables: list of strings containing the variables from the model to be analyzed.
        split = str(self.model_solve).split('\n') #Transform he model results into a list of strings, each element of the list corresponding to a line fo the results (the results are shown in a DataFrame)
        df = []
        for v in variables: #For v in a list containing the variables
            values_minus = []
            for i in range(0,len(self.model_solve[v])-2): 
                x = self.model_solve[v][i+1] - self.model_solve[v][i]
                y = self.model_solve[v][i+2] - self.model_solve[v][i+1]
                if x < 0 and y > 0:
                    values_minus.append(self.model_solve[v][i+1])
            link = {v:values_minus}
            data = pd.DataFrame(link)
            temp_min=[]
            for line in str(self.model_solve).split('\n'):
                for max in data[v].tolist():
                    if str(max) in line:
                        time = line.split(' ')[0]
                        temp_min.append(time)
            dafr = pd.DataFrame(temp_min,data[v].tolist(),['time ' + str(v)])
            dafr.index.name = 'Max Value'
            df.append(dafr)
        return df
    
    def script(self,k,value,dependent=False):
        cut=self.model.split('\n')
        #cut[0] = cut[0] + str(k)
        for line in cut:
            if 'title' in line:
                cut[cut.index(line)] = line + ' ' + str(k) + '=' + str(value)
            if str(k) in line and ',' not in line:
                replace = line.split('=')
                replace[-1]=str(value)
                new = ' = '.join(replace)
                cut[cut.index(line)] = new
                final = '\n'.join(cut)
                return SteadyState(final,tf=self.simtime,change=self.parmchange,npoints=self.resolution,outputs=self.outputs,scan=True)
            if dependent == True:
                if 'init' in line:
                    dep_variable = line.split(':')[-1].split(')')[0].split(',')
                    for variable in dep_variable:
                        if str(k) in variable:
                            var_to_change = variable.split('=')
                            var_to_change[-1] = str(value)
                            dep_variable[dep_variable.index(variable)]='='.join(var_to_change)
                    cut[cut.index(line)] = 'init: ' + ','.join(dep_variable) + ')'
                    final = '\n'.join(cut)
                    return SteadyState(final,tf=self.simtime,change=self.parmchange,npoints=self.resolution,outputs=self.outputs,scan=True)
    