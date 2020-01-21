import csv
import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pylab as plt
from threading import Thread
import math
from random import sample
import statistics as st
import json
import random

def reader(filename):
    master = []
    with open(filename, mode= 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter = ",")
        for row in csv_reader:
            if(row[0] != ''):
                master.append(row)
    return master    

cancer_data  = reader('cancer_data.csv')
edge_list    = reader('edge_list.csv')
healthy_data = reader('healthy_data.csv')
genome_id    = reader('id_uniprot_ensembl_map.csv')

def dict_maker(data):
    master = {}
    firstline = data[0]
    first = True
    for row in data[1:]:
        if int(row[0] not in master.keys()):
            master[int(row[0])] = {}
            for column in firstline:
                try:
                    master[int(row[0])][column] = float(row[firstline.index(column)])
                except:
                    master[int(row[0])][column] = row[firstline.index(column)]
                    
    return master

cancer_dict = dict_maker(cancer_data)
healthy_dict = dict_maker(healthy_data)
genome_id_dict = dict_maker(genome_id)

cancer_att = {}
healthy_att = {}
cancer_fail = set()
healthy_fail = set()

for i in cancer_dict.keys():
    for attribute in cancer_dict[i].keys():
        if(attribute != 'int_id'):
            try:
                if attribute not in cancer_att.keys():
                    cancer_att[attribute] = float(cancer_dict[i][attribute])
                else:
                    cancer_att[attribute]+= float(cancer_dict[i][attribute])
            except:
                cancer_fail.add(attribute)
                pass     
            
for i in healthy_dict.keys():
    for attribute in healthy_dict[i].keys():
        if(attribute != 'int_id'):
            try:
                if(math.isnan(float(healthy_dict[i][attribute])) == False):
                    if attribute not in healthy_att.keys():
                        healthy_att[attribute] = float(healthy_dict[i][attribute])
                    else:
                        healthy_att[attribute] += float(healthy_dict[i][attribute])
                else:
                    healthy_fail.add(attribute)
                    
            except:
                healthy_fail.add(attribute)
                pass

# plt.bar(list(cancer_att.keys()),list(cancer_att.values()))
# plt.xticks(rotation = 'vertical')
# plt.show()
# plt.bar(list(healthy_att.keys()), list(healthy_att.values()))
# plt.xticks(rotation = 'vertical')
# plt.show()


def create_undirected_graph_from_edge_list(edge_list):
    graph = nx.Graph()
    edges = []
    for row in edge_list[1:]:
        edges.append([int(row[0]),int(row[1])])
    graph.add_edges_from(edges)
    return graph



undirected_cancer = create_undirected_graph_from_edge_list(edge_list)
undirected_healthy = create_undirected_graph_from_edge_list(edge_list)

nx.set_node_attributes(undirected_cancer,cancer_dict)
nx.set_node_attributes(undirected_healthy, healthy_dict)


# print("cancer graph analytics")
# print("Number of nodes: ",undirected_cancer.number_of_nodes())
# print("Number of edges: ",undirected_cancer.number_of_edges())
# print("Mean degree: ",2*undirected_cancer.number_of_edges()/ undirected_cancer.number_of_nodes())
# print("Clustering coefficient: ", nx.average_clustering(undirected_cancer))
# print("------------------------")
# print("healthy graph analytics")
# print("Number of nodes: ",undirected_healthy.number_of_nodes())
# print("Number of edges: ",undirected_healthy.number_of_edges())
# print("Mean degree: ",2*undirected_healthy.number_of_edges()/ undirected_healthy.number_of_nodes())
# print("Clustering coefficient: ", nx.average_clustering(undirected_healthy))

def create_subset(input_graph, portion):
    with_attributes = set()
    for node in input_graph.nodes():
        for attribute in input_graph.node[node].keys():
            if(attribute != 'int_id' and attribute != '' and input_graph.node[node][attribute] != None and math.isnan(float(input_graph.node[node][attribute])) == False and input_graph.node[node][attribute] != '' ):
                with_attributes.add(node)
    size = len(list(with_attributes))
    sampling_size = round(size*portion)
    sampled_list = np.random.choice(list(with_attributes),sampling_size)
    return sampled_list
    #return list(with_attributes)

def give_level_1(input_graph,node):
    neighbors = input_graph.neighbors(node)
    return neighbors

def give_level_2(input_graph,node,main_node):
    neighbors = input_graph.neighbors(node)
    return neighbors

def get_uniform_attributes_cancer(cancer_dict,sampled_list):
    uniform_att = {}
    for key in cancer_dict.keys():
        for attribute in cancer_dict[i].keys():
            if(attribute == 'int_id' and cancer_dict[key]['int_id'] in sampled_list):
                continue
            else:
                if(math.isnan(float(cancer_dict[i][attribute])) == False):
                    try:
                        if(attribute not in uniform_att.keys()):
                            uniform_att[attribute] = [float(cancer_dict[i][attribute])]
#                         else:
#                             uniform_att[attribute].append(float(cancer_dict[i][attribute]))
                    except:
                        if(attribute not in uniform_att.keys()):
                            uniform_att[attribute].append(float(cancer_dict[i][attribute]))
                        pass
    for key in uniform_att.keys():
        uniform_att[key] = np.mean(uniform_att[key])
    return uniform_att

def get_uniform_attributes_healthy(healthy_dict,sampled_list):
    uniform_att = {}
    for key in healthy_dict.keys():
        for attribute in healthy_dict[i].keys():
            if(attribute == 'int_id' and healthy_dict[key]['int_id'] in sampled_list):
                continue
            else:
                if(math.isnan(float(healthy_dict[i][attribute])) == False):
                    try:
                        if(attribute not in uniform_att.keys()):
                            uniform_att[attribute] = [float(healthy_dict[i][attribute])]
#                         else:
#                             uniform_att[attribute].append(float(healthy_dict[i][attribute]))
                    except:
                        if(attribute not in uniform_att.keys()):
                            uniform_att[attribute].append(float(healthy_dict[i][attribute]))
#                         pass
    for key in uniform_att.keys():
        uniform_att[key] = np.mean(uniform_att[key])
    return uniform_att

def prediction(input_graph,node_to_predict,level,cancer,sample):
    neighbors = []
    if(level == 1):
        neighbors = give_level_1(input_graph,node_to_predict)
    elif(level == 2):
        direct_neighbors = give_level_1(input_graph,node_to_predict)
        for node in direct_neighbors:
            x = give_level_2(input_graph,node,node_to_predict)
            for i in x:
                neighbors.append(i)
    
    if(neighbors == []):
        if(cancer == True):
            new = get_uniform_attributes_cancer(cancer_dict,sample)
            new['int_id'] = node
            return new
        else:
            new = get_uniform_attributes_healthy(healthy_dict,sample)
            new['int_id'] = node
            return new
    else:
        attribute_D = {}
        for node in neighbors:
            if(node != node_to_predict):
                for attribute in input_graph.node[node].keys():
                    if(attribute != 'int_id' and input_graph.node[node][attribute] and 
                       math.isnan(float(input_graph.node[node][attribute])) == False):
                        if(attribute not in attribute_D.keys()):
                            attribute_D[attribute] = [float(input_graph.node[node][attribute])]
                        else:
                            attribute_D[attribute].append(float(input_graph.node[node][attribute]))
        for key in attribute_D.keys():
            attribute_D[key] = np.mean(attribute_D[key])
        return attribute_D
    return 0

def test_run(graph, portion,cancer):
    returning = []
    sample = create_subset(graph,portion)
    for node in sample:
        #attr = {node:prediction(graph, node, 2,True,sample)}
        #nx.set_node_attributes(graph,attr)
#         print(graph.node[node])
#         print('----')
#         print(attr)
        returning.append((graph.node[node],prediction(graph, node, 1,cancer,sample)))
    return returning

def scalar_values(input_a):
    real_value = []
    calculated = []
    for tup in input_a:
        if 'int_id' in tup[0].keys():
            del tup[0]['int_id']
        original_values = tup[0].values()
        predicted_values = tup[1].values()
        x = 0
        y = 0
        for value in original_values:
            x += float(value)**2
        for value in predicted_values:
            y += float(value)**2
        real_value.append(math.sqrt(x))
        calculated.append(math.sqrt(y))
    return [real_value,calculated]

def as_vector(input_a):
    real_value = {}
    calculated = {}
    for tup in input_a:
        if 'int_id' in tup[0].keys():
            del tup[0]['int_id']
        for key in tup[0].keys():
            if key not in real_value.keys():
                real_value[key] = [tup[0][key]]
            else:
                real_value[key].append(tup[0][key])
        for key in tup[1].keys():
            if key not in calculated.keys():
                calculated[key] = [tup[1][key]]
            else:
                calculated[key].append(tup[1][key])
    return [real_value,calculated]

def numerical_experiment(input_array):
    master = {}
    for item in input_array:
        master[item] = (test_run(undirected_healthy,item,False),test_run(undirected_cancer,item,True))
    return master


def prediction_layer(cancer_graph, healthy_graph, node_to_predict, level, sample, chance):
#     print("in prediction layer")
    neighbors = []
    healthy_neighbors = []
    
    if(level == 1):
        neighbors = give_level_1(cancer_graph, node_to_predict)
        healthy_neighbors = give_level_1(healthy_graph, node_to_predict)
        
    elif(level == 2):
        direct_neighbors = give_level_1(cancer_graph,node_to_predict)
        healthy_neighbors_d = give_level_1(healthy_graph, node_to_predict)
        
        for node in direct_neighbors:
            x = give_level_2(cancer_graph, node, node_to_predict)
            for i in x:
                neighbors.append(i)
        for node in healthy_neighbors_d:
            x = give_level_2(healthy_graph, node, node_to_predict)
            for i in x:
                healthy_neighbors.append(i)
                
    elif(level == 12):
        direct_neighbors = give_level_1(cancer_graph,node_to_predict)
        healthy_neighbors_d = give_level_1(healthy_graph, node_to_predict)
        for node in direct_neighbors:
            x = give_level_2(cancer_graph, node, node_to_predict)
            for i in x:
                neighbors.append(i)
        for node in direct_neighbors:
            neighbors.append(node)
            
        for node in healthy_neighbors_d:
            x = give_level_2(healthy_graph, node, node_to_predict)
            for i in x:
                healthy_neighbors.append(i)
        for node in healthy_neighbors_d:
            healthy_neighbors.append(node)
            
    if(neighbors == []):
#         print("uniformly getting the empirical distribution: ")
        new = get_uniform_attributes_cancer(cancer_dict, sample)
        new['int_id'] = node
        return new
    
    else:
#         print("predicting in general")
        attribute_d = {}
        for node in neighbors:
            if(node != node_to_predict):
                attribute = 'breast_cancer'
                if(attribute in cancer_graph.node[node].keys() and attribute != 'int_id' and cancer_graph.node[node][attribute] and math.isnan(float(cancer_graph.node[node][attribute])) == False):
#                     print('in inner if in cancer predict')
                    if(attribute not in attribute_d.keys()):
                        attribute_d[attribute] = [float(cancer_graph.node[node][attribute])]
                    else:
                        attribute_d[attribute].append(float(cancer_graph.node[node][attribute]))
        if random.randint(0,100) <= chance:
#             print("in the randint if statement")
            for node in healthy_neighbors:
                if( node != node_to_predict):
                    attribute = 'breast'
                    if( attribute in healthy_graph.node[node].keys() and attribute != 'int_id' and healthy_graph.node[node][attribute] and math.isnan(float(healthy_graph.node[node][attribute])) == False):
                        if('breast_cancer' not in attribute_d.keys()):
                            attribute_d['breast_cancer'] = [float(healthy_graph.node[node][attribute])]
                        else:
                            attribute_d['breast_cancer'].append(float(healthy_graph.node[node][attribute]))
        if(attribute_d == {}):
            new = get_uniform_attributes_cancer(cancer_dict, sample)
            new['int_id'] = node
            return new

        attribute_d['breast_cancer'] = np.mean(attribute_d['breast_cancer'])
        return attribute_d
    return 0

        
test_array = [0.05, 0.1, 0.15, 0.2, 0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7]

# final_result = numerical_experiment(test_array)

# with open('final_result.json','w') as outfile:
#     json.dump(final_result, outfile)

def test_run_gen(cancer, healthy, portion, level, chance):
    returning = []
    sample = create_subset(cancer, portion)
    for node in sample:
        returning.append((cancer.node[node],prediction_layer(cancer,healthy,node,level,sample,chance)))
    return returning

def extract_brest_cancer(input_array):
    real = []
    predicted = []
    for tup in input_array:
        real.append(tup[0]['breast_cancer'])
        predicted.append(tup[1]['breast_cancer'])
    return (real,predicted)

with open('final_result.json','r') as json_file:
    data = json.load(json_file)

def find_difference(inputt):
    master = []
    for small_data in inputt:
        for item in small_data:
            if('breast' in item[0].keys() and 'breast' in item[1].keys()):
                master.append(item[0]['breast']-item[1]['breast'])
            elif('breast_cancer' in item[0].keys() and 'breast_cancer' in item[1].keys()):
                master.append(item[0]['breast_cancer']-item[1]['breast_cancer'])
            
#     for item in inputt:
#         print(item[0])
#         if('breast' in item[0].keys()):
#             master.append(item[0]['breast']-item[1]['breast'])
#         else:
#             master.append(item[0]['breast_cancer']-item[1]['breast_cancer'])
    return master

def extract_real_vs_predicted(datat):
    master = {}
    for key in datat.keys():
        master[key] = find_difference(datat[key])
    return master

breast_final = extract_real_vs_predicted(data)

def create_cool_array(number):
    size = int(11482/number)
    x = [0]
    for i in range(0,number-1):
        x.append(i+size)
    return x

def general_plot(input_data):
    colors = ["#ffe6e6","#ffcccc","#ffb3b3","#ff9999","#ff8080",
             "#ff6666","#ff4d4d","#ff3333","#ff1a1a","#ff0000",
             "#e60000","#cc0000","#b30000","#990000"]
    colors_1 = ["#e6194B","#f58231","#ffe119","#bfef45","#3cb44b",
               "#42d4f4","#4363d8","#911eb4","#f032e6","#a9a9a9",
               "#9A6324","#808000","#000000","#469990"]
    color_i = 0
    for i in input_data.keys():
        if(i[-1] == '5'):
            x = create_cool_array(len(input_data[i]))
            print(len(x),len(input_data[i]))
            plt.scatter(x,input_data[i],color=colors_1[color_i],s=[0.25],label = i)
            color_i+=1
    plt.axhline(y=0, color='black', linestyle='-')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.figure(figsize=(20,25))
    plt.show()

# general_plot(breast_final)

with open("final_result_healthy_chance.json",'r') as json_file:
    data = json.load(json_file)

def difference_between_healthy_cancer_chance(healthyd,cancerd):
    difference = []
    if(len(healthyd) != len(cancerd)):
        print("the array size is different")
        return 0
    for i in range(len(healthyd)):
        difference.append(healthyd[i]-cancerd[i])
    return difference

def final_chance_difference(data_input):
    final = {}
    for key in data_input.keys():
        final[key] = difference_between_healthy_cancer_chance(data_input[key][0],data[key][1])
    return final

a = final_chance_difference(data)

x = []
for i in range(0,1865):
    x.append(i)
    
def general_plot(input_data):
    colors = ["#ffe6e6","#ffcccc","#ffb3b3","#ff9999","#ff8080",
             "#ff6666","#ff4d4d","#ff3333","#ff1a1a","#ff0000",
             "#e60000","#cc0000","#b30000","#990000"]
    colors_1 = ["#e6194B","#f58231","#ffe119","#bfef45","#3cb44b",
               "#42d4f4","#4363d8","#911eb4","#f032e6","#a9a9a9"]
              # "#9A6324","#808000","#000000","#469990"]
    color_i = 0
    for i in input_data.keys():
            plt.scatter(x,input_data[i],color=colors_1[color_i],s=[0.1],label = str(i))
            color_i+=1
    plt.axhline(y=0, color='black', linestyle='-')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.figure(figsize=(20,25))
    plt.show()

general_plot(a)

# final_chance = {}
# chance_array = [10,20,30,40,50,60,70,80,90,100]
# for i in chance_array:
#     a = test_run_gen(undirected_cancer,undirected_healthy, 0.2, 1, i )
#     b = extract_brest_cancer(a)
#     final_chance[i] = b 

# with open('final_result_healthy_chance.json','w') as outfile:
#     json.dump(final_chance, outfile)
# for i in test_array:
#     a = test_run_gen(undirected_cancer,undirected_healthy, i, 1, )
# a = test_run_gen(undirected_cancer,undirected_healthy, 0.2, 1, 30)
# b = extract_brest_cancer(a)


