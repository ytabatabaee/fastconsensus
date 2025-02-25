import networkx as nx
import numpy as np
import itertools
import sys
import os
import random
import igraph as ig
from networkx.algorithms import community
import networkx.algorithms.isolate
import community as cm
import math
import random
import argparse
import multiprocessing as mp
import leidenalg
import matplotlib.pyplot as plt


def check_consensus_graph(G, n_p, delta):
    '''
    This function checks if the networkx graph has converged.
    Input:
    G: networkx graph
    n_p: number of partitions while creating G
    delta: if more than delta fraction of the edges have weight != n_p then returns False, else True
    '''



    count = 0

    for wt in nx.get_edge_attributes(G, 'weight').values():
        if wt != 0 and wt != n_p:
            count += 1

    if count > delta*G.number_of_edges():
        return False

    return True



def nx_to_igraph(Gnx):
    '''
    Function takes in a network Graph, Gnx and returns the equivalent
    igraph graph g
    '''
    g = ig.Graph()
    g.add_vertices(sorted(Gnx.nodes()))
    g.add_edges(sorted(Gnx.edges()))
    g.es["weight"] = 1.0
    for edge in Gnx.edges():
        g[edge[0], edge[1]] = Gnx[edge[0]][edge[1]]['weight']
    return g


def group_to_partition(partition):
    '''
    Takes in a partition, dictionary in the format {node: community_membership}
    Returns a nested list of communities [[comm1], [comm2], ...... [comm_n]]
    '''

    part_dict = {}

    for index, value in partition.items():

        if value in part_dict:
            part_dict[value].append(index)
        else:
            part_dict[value] = [index]


    return part_dict.values()

def check_arguments(args):

    if(args.d > 1):
        print('delta is too high. Allowed values are between 0 and 1')
        return False
    if(args.d < 0):
        print('delta is too low. Allowed values are between 0 and 1')
        return False
    if(args.alg not in ('louvain', 'lpm', 'cnm', 'infomap', 'leiden')):
        print('Incorrect algorithm entered. run with -h for help')
        return False
    if (args.t > 1 or args.t < 0):
        print('Incorrect tau. run with -h for help')
        return False

    return True

def communities_to_dict(communities):
    """
    Creates a [node] -> [community] lookup
    """
    result = {}
    community_index = 0
    for c in communities:
        community_mapping = ({str(node):community_index for index, node in enumerate(c)})

        result = {**result, **community_mapping}
        community_index += 1
    return result

def edges_lookup_table_by_node(edges):
    """
    Creates a [node] -> [[u,v]] lookup
    """
    result = {}
    for u, v in edges:
        if u in result:
            result[u].append((u,v))
        else:
            result[u] = [(u,v)]

        if v in result:
            result[v].append((v,u))
        else:
            result[v] = [(v,u)]
    return result


def do_leiden_community_detection(data):
    networkx_graph, seed = data
    return leidenalg.find_partition(nx_to_igraph(networkx_graph), leidenalg.ModularityVertexPartition, weights='weight',  seed=seed, n_iterations=1).as_cover()

def get_graph_and_seed(graph, times):
    for seed in range(times):
        yield graph, seed

def fast_consensus(G,  algorithm = 'louvain', n_p = 20, thresh = 0.2, delta = 0.02):

    graph = G.copy()
    L = G.number_of_edges()
    N = G.number_of_nodes()

    for u,v in graph.edges():
        graph[u][v]['weight'] = 1.0

    iter_count = 0
    while(True):
        iter_count += 1
        print("iter ", iter_count)
        #print(nx.adjacency_matrix(graph, weight='weight'))
        if (algorithm == 'louvain'):

            nextgraph = graph.copy()
            L = G.number_of_edges() # edges in the original network
            for u,v in nextgraph.edges():
                nextgraph[u][v]['weight'] = 0.0

            communities_all = [cm.partition_at_level(cm.generate_dendrogram(graph, randomize = True, weight = 'weight'), 0) for i in range(n_p)]

            print("number of edges in the initial network, ", nextgraph.number_of_edges())
            for node,nbr in graph.edges():
                if (node,nbr) in graph.edges() or (nbr, node) in graph.edges():
                    if graph[node][nbr]['weight'] not in (0,n_p):
                        for i in range(n_p):
                            communities = communities_all[i]
                            if communities[node] == communities[nbr]:
                                nextgraph[node][nbr]['weight'] += 1
                    else:
                        nextgraph[node][nbr]['weight'] = graph[node][nbr]['weight']

            print('step 1 done... Dij built for edges in the network')
            #print(nx.adjacency_matrix(nextgraph, weight='weight'))
            print("number of after Dij, ", nextgraph.number_of_edges())

            remove_edges = []
            for u,v in nextgraph.edges():
                if nextgraph[u][v]['weight'] < thresh*n_p:
                    remove_edges.append((u, v))

            nextgraph.remove_edges_from(remove_edges)

            print('step 2 done... thresholding')
            #print(nx.adjacency_matrix(nextgraph, weight='weight'))
            print("number of edges after thresholding, ", nextgraph.number_of_edges())

            if check_consensus_graph(nextgraph, n_p = n_p, delta = delta):
                break

            for _ in range(L):

                node = np.random.choice(nextgraph.nodes())
                neighbors = [a[1] for a in nextgraph.edges(node)]

                if (len(neighbors) >= 2):
                    a, b = random.sample(set(neighbors), 2)

                    if not nextgraph.has_edge(a, b):
                        nextgraph.add_edge(a, b, weight = 0)

                        for i in range(n_p):
                            communities = communities_all[i]

                            if communities[a] == communities[b]:
                                nextgraph[a][b]['weight'] += 1

            # make sure the graph is connected by fixing singleton nodes
            for node in nx.isolates(nextgraph):
                    nbr, weight = sorted(graph[node].items(), key=lambda edge: edge[1]['weight'])[0]
                    nextgraph.add_edge(node, nbr, weight = weight['weight'])

            print('step 3 done... triadic closure')
            #print(nx.adjacency_matrix(nextgraph, weight='weight'))
            print("number of edges after triadic closure, ", nextgraph.number_of_edges())


            graph = nextgraph.copy()


            if check_consensus_graph(nextgraph, n_p = n_p, delta = delta):
                break

        elif algorithm == 'leiden':

            nextgraph = graph.copy()

            for u,v in nextgraph.edges():
                nextgraph[u][v]['weight'] = 0.0

            #communities = []
            #for i in range(n_p):
            #    communities.append(do_leiden_community_detection((get_graph_and_seed(graph, n_p))))
            with mp.Pool(processes=n_p) as pool:
                communities = pool.map(do_leiden_community_detection, get_graph_and_seed(graph, n_p))

            print("number of edges in the initial network, ", nextgraph.number_of_edges())
            # remember: see if leiden also benefits from the adjustment done in louvain for keeping converged elements in the matrix
            for i in range(n_p):
                node_community_lookup = communities_to_dict(communities[i])
                for community_index, _ in enumerate(communities[i]):
                    for node, nbr in graph.edges():
                        if str(node) in node_community_lookup and str(nbr) in node_community_lookup and node_community_lookup[str(node)] == node_community_lookup[str(nbr)]:
                            if node_community_lookup[str(node)] != community_index:
                                # only count each community once
                                continue
                            nextgraph[node][nbr]['weight'] += 1

            print('step 1 done... Dij built for edges in the network')
            #print(nx.adjacency_matrix(nextgraph, weight='weight'))
            print("number of edges after building Dij, ", nextgraph.number_of_edges())

            remove_edges = []
            for u,v in nextgraph.edges():
                if nextgraph[u][v]['weight'] < thresh*n_p:
                    remove_edges.append((u, v))
            nextgraph.remove_edges_from(remove_edges)

            print('step 2 done... thresholding')
            #print(nx.adjacency_matrix(nextgraph, weight='weight'))
            print("number of edges after thresholding, ", nextgraph.number_of_edges())

            if check_consensus_graph(nextgraph, n_p = n_p, delta = delta):
                break

            for i in range(n_p):
                node_community_lookup = communities_to_dict(communities[i])
                n_graph_nodes = len(nextgraph.nodes())

                edges_lookup_table = edges_lookup_table_by_node(nextgraph.edges)

                for _ in range(L):
                    random_node_index = random.randint(1, n_graph_nodes)
                    neighbors = [a[1] for a in edges_lookup_table.get(str(random_node_index), [])]

                    if (len(neighbors) >= 2):
                        a, b = random.sample(set(neighbors), 2)

                        if not nextgraph.has_edge(a, b):
                            nextgraph.add_edge(a, b, weight = 0)

                            if a in node_community_lookup and b in node_community_lookup and node_community_lookup[a] == node_community_lookup[b]:
                                nextgraph[a][b]['weight'] += 1

            for node in nx.isolates(nextgraph):
                nbr, weight = sorted(graph[node].items(), key=lambda edge: edge[1]['weight'])[0]
                nextgraph.add_edge(node, nbr, weight = weight['weight'])

            print('step 3 done... triadic closure')
            #print(nx.adjacency_matrix(nextgraph, weight='weight'))
            print("number of edges after tradic closure, ", nextgraph.number_of_edges())

            graph = nextgraph.copy()

            if check_consensus_graph(nextgraph, n_p = n_p, delta = delta):
                break

        elif (algorithm in ('infomap', 'lpm')):

            nextgraph = graph.copy()

            for u,v in nextgraph.edges():
                nextgraph[u][v]['weight'] = 0.0

            if algorithm == 'infomap':
                communities = [{frozenset(c) for c in nx_to_igraph(graph).community_infomap().as_cover()} for _ in range(n_p)]
            if algorithm == 'lpm':
                communities = [{frozenset(c) for c in nx_to_igraph(graph).community_label_propagation().as_cover()} for _ in range(n_p)]


            for node, nbr in graph.edges():

                for i in range(n_p):
                    for c in communities[i]:
                        if node in c and nbr in c:
                            if not nextgraph.has_edge(node,nbr):
                                nextgraph.add_edge(node, nbr, weight = 0)
                            nextgraph[node][nbr]['weight'] += 1



            remove_edges = []
            for u,v in nextgraph.edges():
                if nextgraph[u][v]['weight'] < thresh*n_p:
                    remove_edges.append((u, v))
            nextgraph.remove_edges_from(remove_edges)



            for _ in range(L):
                node = np.random.choice(nextgraph.nodes())
                neighbors = [a[1] for a in nextgraph.edges(node)]

                if (len(neighbors) >= 2):
                    a, b = random.sample(set(neighbors), 2)

                    if not nextgraph.has_edge(a, b):
                        nextgraph.add_edge(a, b, weight = 0)

                        for i in range(n_p):
                            if a in communities[i] and b in communities[i]:
                                nextgraph[a][b]['weight'] += 1


            graph = nextgraph.copy()

            if check_consensus_graph(nextgraph, n_p = n_p, delta = delta):
                break

        elif (algorithm == 'cnm'):

            nextgraph = graph.copy()

            for u,v in nextgraph.edges():
                nextgraph[u][v]['weight'] = 0.0

            communities = []
            mapping = []
            inv_map = []


            for _ in range(n_p):

                order = list(range(N))
                random.shuffle(order)
                maps = dict(zip(range(N), order))

                mapping.append(maps)
                inv_map.append({v: k for k, v in maps.items()})
                G_c = nx.relabel_nodes(graph, mapping = maps, copy = True)
                G_igraph = nx_to_igraph(G_c)

                communities.append(G_igraph.community_fastgreedy(weights = 'weight').as_clustering())


            for i in range(n_p):

                edge_list = [(mapping[i][j], mapping[i][k]) for j,k in graph.edges()]

                for node,nbr in edge_list:
                    a, b = inv_map[i][node], inv_map[i][nbr]

                    if graph[a][b] not in (0, n_p):
                        for c in communities[i]:
                            if node in c and nbr in c:
                                nextgraph[a][b]['weight'] += 1

                    else:
                        nextgraph[a][b]['weight'] = graph[a][b]['weight']


            remove_edges = []
            for u,v in nextgraph.edges():
                if nextgraph[u][v]['weight'] < thresh*n_p:
                    remove_edges.append((u, v))

            nextgraph.remove_edges_from(remove_edges)


            for _ in range(L):
                node = np.random.choice(nextgraph.nodes())
                neighbors = [a[1] for a in nextgraph.edges(node)]

                if (len(neighbors) >= 2):
                    a, b = random.sample(set(neighbors), 2)
                    if not nextgraph.has_edge(a, b):
                        nextgraph.add_edge(a, b, weight = 0)

                        for i in range(n_p):
                            for c in communities[i]:
                                if mapping[i][a] in c and mapping[i][b] in c:

                                    nextgraph[a][b]['weight'] += 1

            if check_consensus_graph(nextgraph, n_p, delta):
                break

        else:
            break

    if (algorithm == 'louvain'):
        return [cm.partition_at_level(cm.generate_dendrogram(graph, randomize = True, weight = 'weight'), 0) for _ in range(n_p)]
    if algorithm == 'leiden':
        with mp.Pool(processes=n_p) as pool:
            communities = pool.map(do_leiden_community_detection, get_graph_and_seed(graph, n_p))
        return communities
    if algorithm == 'infomap':
        return [{frozenset(c) for c in nx_to_igraph(graph).community_infomap().as_cover()} for _ in range(n_p)]
    if algorithm == 'lpm':
        return [{frozenset(c) for c in nx_to_igraph(graph).community_label_propagation().as_cover()} for _ in range(n_p)]
    if algorithm == 'cnm':

        communities = []
        mapping = []
        inv_map = []

        for _ in range(n_p):
            order = list(range(N))
            random.shuffle(order)
            maps = dict(zip(range(N), order))

            mapping.append(maps)
            inv_map.append({v: k for k, v in maps.items()})
            G_c = nx.relabel_nodes(graph, mapping = maps, copy = True)
            G_igraph = nx_to_igraph(G_c)

            communities.append(G_igraph.community_fastgreedy(weights = 'weight').as_clustering())

        return communities


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('-f', metavar='filename', type=str, nargs = '?', help='file with edgelist')
    parser.add_argument('-np', metavar='n_p', type=int, nargs = '?', default=20, help='number of input partitions for the algorithm (Default value: 20)')
    parser.add_argument('-t', metavar='tau', type=float, nargs = '?', help='used for filtering weak edges')
    parser.add_argument('-d', metavar='del', type=float,  nargs = '?', default = 0.02, help='convergence parameter (default = 0.02). Converges when less than delta proportion of the edges are with wt = 1')
    parser.add_argument('--alg', metavar='alg', type=str, nargs = '?', default = 'louvain' , help='choose from \'louvain\' , \'cnm\' , \'lpm\' , \'infomap\' ')

    args = parser.parse_args()

    default_tau = {'louvain': 0.2, 'cnm': 0.7 ,'infomap': 0.6, 'lpm': 0.8}
    if (args.t == None):
        args.t = default_tau.get(args.alg, 0.2)

    if check_arguments(args) == False:
        quit()

    G = nx.read_edgelist(args.f, nodetype=int)

    # plot the network
    '''fig = plt.figure(figsize=(30, 30), dpi=80)
    nx.draw(G, with_labels=True, font_size=8, node_color='lightblue')
    fig.savefig("network.pdf")'''

    # relabel nodes
    mapping = dict(zip(G, range(0, G.number_of_nodes())))
    G = nx.relabel_nodes(G, mapping)

    output = fast_consensus(G, algorithm = args.alg, n_p = args.np, thresh = args.t, delta = args.d)

    out_partitions_path = 'out_partitions_t' + str(args.t) + '_d' + str(args.d) + '_np' + str(args.np)
    membership_path = 'memberships_t' + str(args.t) + '_d' + str(args.d) + '_np' + str(args.np)


    if not os.path.exists(out_partitions_path):
        os.makedirs(out_partitions_path)

    if not os.path.exists(membership_path):
        os.makedirs(membership_path)

    if(args.alg == 'louvain'):
        for i in range(len(output)):
            with open(membership_path + '/' + str(i), 'w') as f:
                for k, v in sorted(output[i].items()):
                    f.write(str(k+1) + "\t" + str(v+1) + '\n')
            output[i] = group_to_partition(output[i])

    i = 0
    for partition in output:
        i += 1
        with open(out_partitions_path + '/' + str(i) , 'w') as f:
            for community in partition:
                print(*community, file = f)
        if args.alg == 'leiden':
            with open(out_partitions_path + '/' + str(i), 'w') as f:
                for j in range(len(partition.membership)):
                    f.write(str(j+1) + "\t" + str(partition.membership[j][0]+1) + '\n')
