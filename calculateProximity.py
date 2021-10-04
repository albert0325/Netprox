import shortestPath as shpath #import shortestPath methods
import pandas as pd
import statistics
import networkx as nx

#Calculate the proximity ds of a lncRNA in module N1 (lncN1_1)
def get_proximity_shortest(G_interactome, dic_DrugTarget, ls_mod_PPI):
    dic_Tx_score = {}
    for drug in dic_DrugTarget:
        dsN1_1 = 0
        for i in dic_DrugTarget[drug]:
            ls = [] #To collect protein of dic_DugTarget inside G_interactome.node
            for j in ls_mod_PPI:
                if i in G_interactome.nodes and j in G_interactome.nodes and nx.has_path(G_interactome, i,j) == True:
                    dsN1_1 += shpath.get_shortest_path_length (G_interactome, i, j)
                    ls.append(i)           
        if len(ls) != 0: #To exclude proteins not in the interactome
            dic_Tx_score[drug] = dsN1_1/(len(dic_DrugTarget[drug])*len(ls_mod_PPI))

    #sort the proximity ds according to ds
    sort_orders = sorted(dic_Tx_score.items(), key=lambda x: x[1], reverse=False)

    df_shortest = pd.DataFrame.from_dict(sort_orders)
    return df_shortest
def get_proximity_shortest_1(G_interactome, dic_DrugTarget, ls_mod_PPI): #no sort
    dic_Tx_score = {}
    for drug in dic_DrugTarget:
        dsN1_1 = 0
        ls = [] #To collect protein of dic_DugTarget inside G_interactome.node
        for i in dic_DrugTarget[drug]:
            
            for j in ls_mod_PPI:
                if i in G_interactome.nodes and j in G_interactome.nodes and nx.has_path(G_interactome, i,j) == True:
                    dsN1_1 += shpath.get_shortest_path_length (G_interactome, i, j)
                    ls.append(i)           
        if len(ls) != 0: #To exclude proteins not in the interactome
            dic_Tx_score[drug] = [dsN1_1/(len(dic_DrugTarget[drug])*len(ls_mod_PPI))]

    #sort the proximity ds according to ds
    # sort_orders = sorted(dic_Tx_score.items(), key=lambda x: x[1], reverse=False)

    df_shortest = pd.DataFrame.from_dict(dic_Tx_score, orient='index', columns=["proximity"])
    return df_shortest
# Calculate the closest proximity 
def get_proximity_closest(G_interactome,dic_DrugTarget,ls_mod_PPI):
    dic_Tx_score={}
    for drug in dic_DrugTarget:
        d_min=[]
        for i in dic_DrugTarget[drug]:
            d=[]
            for j in ls_mod_PPI:
                if i in G_interactome.nodes and j in G_interactome.nodes and nx.has_path(G_interactome, i,j) == True:
                    d.append(shpath.get_shortest_path_length (G_interactome, i, j))
            if len(d)!=0:
                d_min.append(min(d))
        if len(d_min) > 1:
            dic_Tx_score[drug]=statistics.mean(d_min)
        
        elif len(d_min) == 1:
            dic_Tx_score[drug]= d_min[0]
    
    sort_orders = sorted(dic_Tx_score.items(), key=lambda x: x[1], reverse=False)
    df_shortest = pd.DataFrame.from_dict(sort_orders)
    return df_shortest

def get_proximity_closest_1(G_interactome,dic_DrugTarget,ls_mod_PPI):#no sort
    dic_Tx_score={}
    for drug in dic_DrugTarget:
        d_min=[]
        for i in dic_DrugTarget[drug]:
            d=[]
            for j in ls_mod_PPI:
                if i in G_interactome.nodes and j in G_interactome.nodes and nx.has_path(G_interactome, i,j) == True:
                    d.append(shpath.get_shortest_path_length (G_interactome, i, j))
            if len(d)!=0:
                d_min.append(min(d))
        if len(d_min) > 1:
            dic_Tx_score[drug]=[statistics.mean(d_min)]
        
        elif len(d_min) == 1:
            dic_Tx_score[drug]= [d_min[0]]
    
    # sort_orders = sorted(dic_Tx_score.items(), key=lambda x: x[1], reverse=False)
    df_shortest = pd.DataFrame.from_dict(dic_Tx_score, orient='index',columns=["proximity"])
    return df_shortest
# Calculate the closest proximity 
def get_proximity_closest_for_z(G_interactome,ls_grand_drug_rand_protein,ls_lnc_rand_protein, ls_drug):
    dic_Tx_score={}
    d_min = []
    for position in range(len(ls_drug)):
        for i in ls_grand_drug_rand_protein[position]:
            d=[]
            for j in ls_lnc_rand_protein:
                d.append(shpath.get_shortest_path_length (G_interactome, i, j))
            if len(d)!=0:
                d_min.append(min(d))
        if len(d_min) > 1:
            dic_Tx_score[ls_drug[position]]=statistics.mean(d_min)
        
        elif len(d_min) == 1:
            dic_Tx_score[ls_drug[position]]= d_min[0]
    
    df_shortest = pd.DataFrame.from_dict(dic_Tx_score.items())
    return list(df_shortest.iloc[:,0]), list(df_shortest.iloc[:,1]) #[0]為drug; [1]為proximity    

# Calculate the weighted closest proximity
def get_proximity_weight_closest(G_weight_interactome,dic_DrugTarget,ls_mod_PPI):
    dic_Tx_score={}
    for drug in dic_DrugTarget:
        d_min=[]
        for i in dic_DrugTarget[drug]:
            d=[]
            for j in ls_mod_PPI:
                if i in G_weight_interactome.nodes and j in G_weight_interactome.nodes and nx.has_path(G_weight_interactome, i,j) == True:
                    d.append(nx.shortest_path_length(G_weight_interactome,i,j,weight='weight'))
            if len(d)!=0:d_min.append(min(d))
        if len(d_min) > 1:
            dic_Tx_score[drug]=statistics.mean(d_min)
        elif len(d_min) == 1:
            dic_Tx_score[drug]= d_min[0]

    sort_orders = sorted(dic_Tx_score.items(), key=lambda x: x[1], reverse=False)
    df_shortest = pd.DataFrame.from_dict(sort_orders)
    return df_shortest    

# Calculate the weighted closest proximity. single source dijkstra_path_length (calculate the shortest path one by one)
def get_proximity_weight_closest_dijkstra_path_length(G_weight_interactome,dic_DrugTarget,ls_mod_PPI):
    dic_Tx_score={}
    for drug in dic_DrugTarget:
        d_min=[]
        for i in dic_DrugTarget[drug]:
            d=[]
            for j in ls_mod_PPI:
                if i in G_weight_interactome.nodes and j in G_weight_interactome.nodes:
                    d.append(nx.dijkstra_path_length(G_weight_interactome,i,j,weight='weight'))
            if len(d)!=0:d_min.append(min(d))
        if len(d_min) > 1:
            dic_Tx_score[drug]=statistics.mean(d_min)
        elif len(d_min) == 1:
            dic_Tx_score[drug]= d_min[0]

    sort_orders = sorted(dic_Tx_score.items(), key=lambda x: x[1], reverse=False)
    df_shortest = pd.DataFrame.from_dict(sort_orders)
    return df_shortest    

#Multi source dijkstra_path_length (calculate the shortest path with multiple source as input)
def get_proximity_weight_closest_dijkstra_path_length_multisource(G_weight_interactome,dic_DrugTarget,ls_mod_PPI):
    #multiple source: lncPPI as a whole; target: drug targets one by one 
        
    dic_Tx_score={}
    for drug in dic_DrugTarget:
        d_min=[]
        set_source_node = set()
        for j in ls_mod_PPI:
            if j in G_weight_interactome.nodes:
                set_source_node.add(j)
        for i in dic_DrugTarget[drug]:
            if i in G_weight_interactome.nodes:
                length, path = nx.multi_source_dijkstra(G_weight_interactome,set_source_node,i,weight='weight')
                d_min.append(length)
        if len(d_min) > 1:
            dic_Tx_score[drug]=statistics.mean(d_min)
        elif len(d_min) == 1:
            dic_Tx_score[drug]= d_min[0]

    sort_orders = sorted(dic_Tx_score.items(), key=lambda x: x[1], reverse=False)
    df_shortest = pd.DataFrame.from_dict(sort_orders)
    return df_shortest    

def get_proximity_weight_closest_dijkstra_path_length_multisource_for_z(G_weight_interactome,ls_grand_drug_rand_protein,ls_lnc_rand_protein, ls_drug):
    #multiple source: lncPPI as a whole; target: drug targets one by one 
        
    dic_Tx_score={}
    set_source_node = set()
    d_min = []
        
    for j in ls_lnc_rand_protein:
        set_source_node.add(j)
    for position in range(len(ls_grand_drug_rand_protein)): #200多個藥物
        for drugTargetProtein in ls_grand_drug_rand_protein[position]:            
            length, path = nx.multi_source_dijkstra(G_weight_interactome,set_source_node,drugTargetProtein,weight='weight')
            d_min.append(length)
        if len(d_min) > 1:
            dic_Tx_score[ls_drug[position]]=statistics.mean(d_min)
        elif len(d_min) == 1:
            dic_Tx_score[ls_drug[position]]= d_min[0]

    #print(dic_Tx_score)
    df_shortest = pd.DataFrame.from_dict(dic_Tx_score.items())
    return list(df_shortest.iloc[:,0]), list(df_shortest.iloc[:,1]) #[0]為drug; [1]為proximity    

