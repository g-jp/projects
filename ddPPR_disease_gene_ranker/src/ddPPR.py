import igraph as ig
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from tqdm import tqdm
import plotly.express as px
import plotly.graph_objects as go
from statsmodels.stats.multitest import multipletests

def ddPPR(g:ig.Graph,ctrl_dict,test_dict,background,directed,iterations:int,alpha:float,gamma:float,attribute_display:str,pathway=False,path_df=pd.DataFrame(),graph_path_field=""):

    seed_set_c = ctrl_dict["seed"]
    seed_set_t = test_dict["seed"]
    
    seeds = seed_set_c + seed_set_t
    node_ids = g.vs[attribute_display] # get the node IDs from the graph

    # Initiate control and test restart vectors
    restart_vector_c = np.full(
        g.vcount(),
        background
    )
    restart_vector_t = np.full(
        g.vcount(),
        background
    )

    # Load the seed restart weights in both vectors
    if seed_set_c:
        restart_vector_c[seed_set_c] = ctrl_dict["score"] # attribute membrane seeds full restart p
    if seed_set_t:   
        restart_vector_t[seed_set_t] = test_dict["score"] # attribute non-membrane seeds lower restart p 
        
    print(f"The control non-normalized restart vector entries sum to {restart_vector_c.sum()}")
    print(f"The test non-normalized restart vector entries sum to {restart_vector_t.sum()}")

    #norm = max(restart_vector_c.sum(),restart_vector_t.sum())
    restart_vector_c /= restart_vector_c.sum() # Normalize to sum to 1
    restart_vector_t /= restart_vector_t.sum()
    
    # Run the PPR algorithm
    es_pagerank_scores_c = np.array( # convert to np array for flexibility
        g.personalized_pagerank( # personalized page rank
        reset=restart_vector_c,
        damping=alpha, # damping probability
        directed=directed))
    es_pagerank_scores_t = np.array( 
        g.personalized_pagerank( 
        reset=restart_vector_t,
        damping=alpha,
        directed=directed))

    # Run a control PPR with unbiased restart probabilities
    ctrl_restart_vector = np.full(g.vcount(),
                              background) # restart vector for the random walk

    ctrl_pagerank_scores = np.array(g.personalized_pagerank(reset=ctrl_restart_vector, damping=alpha, directed=directed))

    dppr_c = (es_pagerank_scores_c - ctrl_pagerank_scores) # calculate the difference between the two pagerank scores
    dppr_t = (es_pagerank_scores_t - ctrl_pagerank_scores)
    ddppr = (dppr_t - dppr_c) / np.abs(dppr_c)
    
    # Perform pathway analysis 
    if pathway:
        path_df["Which_genes"] = pd.Series([[] for _ in range(len(path_df))]) # initialize a column to store the genes in each pathway
        path_ddppr = np.zeros(len(path_df))
        for node in g.vs:
            pathways = node[graph_path_field] # fetch the pathways in which the node is included
            pathways = [pathways] if not isinstance(pathways,list) else pathways # check whether it is a list of pathways
            annotation = len(pathways)            
            pathways_with_node = path_df[graph_path_field].isin(pathways) # mask to retireve the associated pathways
            path_ddppr[pathways_with_node] += ddppr[node.index] # add the score per pathway
            for idx in path_df.index[pathways_with_node]:
                path_df.at[idx,"Which_genes"].append(node_ids[node.index]) # add the node ID to the pathway genes list

    # Run multiple network rewirings to calculate significance scores

    null_ddppr_scores = np.zeros((g.vcount(), iterations)) # Matrix initialization to collect null difference scores from the rewired graph
    null_path_ddppr_scores = np.zeros((len(path_df),iterations))
    for i in tqdm(range(iterations)):
        # Rewire the graph while preserving node degrees
        rewired_graph = g.copy()
        rewired_graph.rewire(mode="simple", n=10 * rewired_graph.ecount())
        
        # Run PPR on the rewired graph using the same restart vector
        null_ppr_c = np.array(rewired_graph.personalized_pagerank(reset=restart_vector_c, damping=alpha, directed=False))
        null_ppr_t = np.array(rewired_graph.personalized_pagerank(reset=restart_vector_t, damping=alpha, directed=False))
        null_ctrl_ppr = np.array(rewired_graph.personalized_pagerank(reset=ctrl_restart_vector, damping=alpha, directed=False))
        null_dppr_c = (null_ppr_c - null_ctrl_ppr) 
        null_dppr_t = (null_ppr_t - null_ctrl_ppr) 
        null_ddppr = (null_dppr_t - null_dppr_c) / np.abs(null_dppr_c)
        null_ddppr_scores[:, i] = null_ddppr

        # perform pathway enrichment analysis significance evaluation
        if pathway:
            null_path_ddppr = np.zeros(len(path_df))
            for node in g.vs:
                pathways = node[graph_path_field] # fetch the pathways in which the node is included
                pathways = [pathways] if not isinstance(pathways,list) else pathways # check whether it is a list of pathways
                annotation = len(pathways)
                pathways_with_node = path_df[graph_path_field].isin(pathways)
                null_path_ddppr[pathways_with_node] += null_ddppr[node.index]
            null_path_ddppr_scores[:,i] = null_path_ddppr

    # Compute z-scores and p-values
    p_val_two_tailed = np.sum(np.abs(null_ddppr_scores) >= np.abs(ddppr[:,np.newaxis]),axis=1)/iterations
    _, fdr_corrected_pvalue, _, _ = multipletests(p_val_two_tailed, method='fdr_bh')

    if pathway:
        # Compute z-scores and p-values for pathway enrichment
        path_p_val_two_tailed = np.sum(np.abs(null_path_ddppr_scores) >= np.abs(path_ddppr[:,np.newaxis]),axis=1)/iterations
        _, path_fdr_corrected_pvalue, _, _ = multipletests(path_p_val_two_tailed, method='fdr_bh')   
        path_df["Enrichment"] = path_ddppr
        path_df["P_value"] = path_p_val_two_tailed
        path_df["FDR"] = path_fdr_corrected_pvalue

    final_score = pd.DataFrame(data={ # create a df with the scores and p-values
        "UniprotID":node_ids,
        "Control_PPR":es_pagerank_scores_c,
        "Test_PPR":es_pagerank_scores_t,
        "dPPR_C":dppr_c,
        "dPPR_t":dppr_t,
        "ddPPR":ddppr,
        "P_value_2t":p_val_two_tailed,
        "FDR":fdr_corrected_pvalue
    })
    
    # add an Is Seed annotation to each node
    final_score["Is seed"] = 0 
    for node in final_score.index:
        if node in seed_set_c and node not in seed_set_t:
            final_score.loc[node,"Is seed"] = "Control"
        elif node in seed_set_t and node not in seed_set_c:
            final_score.loc[node,"Is seed"] = "Test"
        elif node in seed_set_c and node in seed_set_t:
            final_score.loc[node,"Is seed"] = "Control + Test"
        else:
            final_score.loc[node,"Is seed"] = "Nonseed"
    
    if pathway:
        return {"ddPPR":final_score, "Null Score":null_ddppr_scores,"path_ddPPR":path_df} 
    
    return {"ddPPR":final_score, "Null Score":null_ddppr_scores}
