from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

def union(list1,list2): # performs a union operation between both lists passed as arguments
    unique_list = []
    for element in list1:
        if element not in unique_list:
            unique_list.append(element)
        else:
            continue
    for element in list2:
        if element not in unique_list:
            unique_list.append(element)
        else:
            continue
    return unique_list

def unique(doubled_list):
    unique_list = []
    for element in doubled_list:
        if element in unique_list:
            continue
        unique_list.append(element)
    return unique_list

def exclude(list1,list2): # retrieve elements in list1 that are not present in list2
    # list1 - list2
    return unique([element for element in list1 if element not in list2])

def intersect(list1,list2): # intersection between two sets of values
    intersection_list = list()
    for element in list1:
        if element in list2 and element not in intersection_list:
            intersection_list.append(element)

    return unique([element for element in list1 if element in list2 and element])

def value_count(my_list):
    list_count = {}
    for index in range(len(my_list)):
        count = 0
        for element in list:
            if element == my_list[index]:
                count += 1
        list_count[my_list[index]] = count
    
    return list_count

#______________________________________________________________________________

def add_interactors(graph,subgraph,int_list,attribute:str,show_statistic=False):
    maingraph_intersect = intersect(graph.vs[attribute],int_list) # intersection between interactors and the general graph
    int_not_in_graph = exclude(int_list,maingraph_intersect)
    print(f"{len(int_not_in_graph)} experimental interactors have no interaction described in the current model")
    subgraph_intersect = intersect(subgraph.vs[attribute],int_list) # intersection between interactors and the specific graph
    interactors_to_add = exclude(maingraph_intersect,subgraph_intersect) # remove already existing interactors
    seed_set = which_attribute(union(maingraph_intersect,subgraph_intersect),graph,"name")
    if show_statistic:
        cent_measure = np.array(subgraph.betweenness())[seed_set]
        fig, ax = plt.subplots(1,1,figsize=(7,5))

        ax.boxplot(cent_measure)
        plt.show()
    new_int = {}
    i = 0
    for interactor in interactors_to_add:
        index = graph.vs.find(name=interactor).index # get the interactor index
        name = graph.vs[index][attribute]
  
        neighbors = graph.neighbors(index,mode="ALL") # retrieve the interactor's neighbors in the general graph
        neig_in_subgraph = intersect(subgraph.vs[attribute],graph.vs[neighbors][attribute]) # check which neighbors are already in the subgraph

        if not neig_in_subgraph:
            print(f"{graph.vs[index][attribute]} has no neighbor in the selected subgraph")
            continue
        else:
            subgraph.add_vertices(name)
            edge_list = []
            for neighbor in neig_in_subgraph:
                edge_list.append((name,neighbor))
            subgraph.add_edges(edge_list)
            i += 1
        new_int[name] = neig_in_subgraph
    print(f"{i} nodes were added")
    return new_int

#______________________________________________________________________________
# function to obtain the index of the biggest component of the current graph
def which_max(graph):
    biggest_component = max([len(x) for x in graph.components()])
    for element in graph.components():
        if len(element) == biggest_component:
            index = list(graph.components()).index(element)
            return index
        
def which_index(node_list,graph,attribute_name:str): # retrieve the attribute of the nodes in a list based on a specific index in the graph
    attribute_list = []
    for node in node_list:
        attribute_list.append(graph.vs[node][attribute_name])
    return attribute_list

def which_attribute(attribute_list,graph,attribute): # retrieve the indexes of the nodes in a list based on a specific attribute in the graph
    node_index = [v.index for v in graph.vs if v[attribute] in attribute_list]
    return node_index

def get_seeds(seed_list,graph,percentile=95):
    centralities = pd.Series(graph.degree())[seed_list]
    threshold = np.percentile(centralities,percentile)
    filtered_seeds = centralities[centralities < threshold].index.to_list()
    return filtered_seeds

# annotate the nodes in an igraph object with the pathways it belongs to according to its gene label
def add_pathway_annotation(graph,pathway_df,graph_field:str,df_gene_field,df_path_field:str,inplace=False,voice=False):
    annotated_nodes = {}
    for node in graph.vs:
        if node[graph_field] in pathway_df[df_gene_field].to_list():
            annotated_nodes[node.index] = pathway_df.loc[pathway_df[df_gene_field]==node[graph_field],df_path_field].to_list()[0]
        else:
            continue
    if inplace:
        graph.vs[annotated_nodes.keys()][df_path_field] = [name for name in annotated_nodes.values()] # in-place modification of graph fields
    
    print(f"{len(annotated_nodes.values())} genes were annotated")
    return annotated_nodes
