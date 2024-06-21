from stimator import read_model
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from collections import Counter
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import math as m

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
# Objects in python: also called instances, objects are general limits that contain some amount of code and data that will be further processed. They make up python's 
#Programming language and the interconversion of objects is precisely what defines coding.
# When creating a method or a constructor, the "self" argument is used so as for the function to act upon the instance. Basically, having "self" as an argument is the same saying
# Class.Method(Instance), which is summarized as Instance.Method() (these are two equivalent forms; the method acts upon the instance.
class SteadyState():
    def __init__(self,model,change={},tf=1000,npoints=1000,scan=False): # A constructor that is used to insert instance variables (i.e., variables that are particular to the created object)
        self.model = model
        self.model_solve = read_model(model).solve(tf=tf,npoints=npoints) #Creates the "model" instance variable that is the solved model 
        if scan == True:
            self.model = model
            self.model_scan = read_model(model).scan(change,tf=tf,npoints=npoints) #Creates the "model" instance variable that is the scanned model
            self.model_solve = read_model(model).solve(tf=tf,npoints=npoints) #Creates the "model" instance variable that is the solved model 
        #print(f'St+eadyState constructed')
    def steady_scan(self,lists,stimuli,sensibility=False):
        sa_values = {}#Empty dictionary for sensibility values 
        ss_values = {}#Empty dictionary for steady-state values
    # The scan object is a result of applying the scan function to the model
    #The list requires is a list containing a string with the variables to be analyzed
    #The stimuli list is a list contaning different percentage values of change for a paramter as to determine the sensibility of the system to that parameter
    #The sensibility term is an option so as to decide if the sensibility analysis will be performed
        for n in lists: #For a variable in a list containing all variable in strings
            for x in range(0,len(self.model_scan)): #For a number assigined to all tested conditions
                ss=self.model_scan[x][n][-1] #Iterates through each test condition and each variable, evaluating it's Steady-State(SS) value
                ss_values.update({n+str(x):ss})  #Adds this value to a dictionary
        if sensibility == True: #Condition to peform a sensibility analysis
            for n in lists:
                for x in range(0,len(self.model_scan)):
                    ssn=ss_values[n+str(x)] #Steady-state in a specific condition
                    ss0=ss_values[n+'0'] #Steady-state in the first "normal" condition
                    #print(ssn,ss0)
                    sa=((ssn-ss0)/ss0)/(stimuli[x]-1) #Calculated sensibility change due to changes in a parameter (defined by the percentage change in this parameter, B)
                    sa_values.update({n+str(x):sa})
            return pd.DataFrame([ss_values, sa_values],index=['SS','SA']).transpose() #Creates a DataFrame from the SS values and the SA values for each variable in each tested condition
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
    
    def script(self,k,value,tf=1000,npoints=10000): #Changes a parameter of the model to any number requested and creates another object with the new parameter
        # k = parameter to be changed
        # value = new value of the parameter
        # tf = final time of the simulation
        # npoints = total amount of points calculated
        cut=self.model.split('\n')
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
            elif 'init' in n:
                init = n.split(',')
                for element in init:
                    if str(k) in element:
                        var_change = element.strip().split('=')
                        print(var_change)
                        var_change[-1] = str(value)
                        new = '='.join(var_change)
                        print(new)
                        replace[replace.index(element)] = new
                        replace = ','.join(replace)
                        cut[cut.index(n)] = replace
                        print(cut)
                        final = '\n'.join(cut)
            
        return SteadyState(final,tf=tf,npoints=npoints)
        
