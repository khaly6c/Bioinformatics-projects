
import pickle
import networkx as nx
 
#on load le fichier
with open('1S72.pickle', 'rb') as f:
    graph = pickle.load(f)

#Comment aller a travers tout le graph ET regarder le dictionaire data associe avec
for node, data in graph.nodes(data=True):
    #Un noeud est TOUJOURS un tuple (position, chain)
    position, chain = node

#On regarde si la position est modifiee chimiquement
if data['nucleotide'] not in ('A', 'C', 'G', 'U'):
    print(f"{node=}, {data['nucleotide']=}, {chain=}, {position=}")

#Maintenant on va chercher les aretes CHS
for source, target, data in graph.edges(data=True):
    if data['label'] == 'CHH':
        print(f"{source=}, {target=}, {data['label']=}")

 