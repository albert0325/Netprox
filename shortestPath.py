import pandas as pd
import numpy as np
import networkx as nx
import csv

def create_graph_from_edgelist(path):
    with open(path) as infile:
        interactome = csv.reader(infile)
        graph = nx.Graph(interactome)
    return graph

#def get_shortest_path (graph, source, target):
#    return nx.shortest_path(graph,source,target)

def get_all_shortest_path (graph, source, target):
    return ([p for p in nx.all_shortest_paths(graph,source,target)])

def get_shortest_path_length (graph, source, target):
    return nx.shortest_path_length(graph,source,target)
