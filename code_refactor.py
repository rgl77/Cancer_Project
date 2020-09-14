import pandas as pd 
import networkx as nx
import numpy as np

class DATA_OBJECT(object):
	def __init__(self,cancer,healthy,genome,edge_list):
		self.cancerous_data  = self.read_csv_file(cancer)
		self.healthy_data    = self.read_csv_file(healthy)
		self.genome_id_data  = self.read_csv_file(genome)
		self.edge_list       = self.read_csv_file(edge_list)

		self.cancer_dict     = self.dictionary_maker(self.cancerous_data)
		self.healthy_dict    = self.dictionary_maker(self.healthy_data)
		self.genome_id_dict  = self.dictionary_maker(self.genome_id_data)

		self.cancer_fail = set()
		self.healthy_fail= set()
		self.cancer_att  = self.populate_attributes(self.cancer_dict,  self.cancer_fail)
		self.healthy_att = self.populate_attributes(self.healthy_dict, self.healthy_fail)

		self.undirected_cancer = self.create_undirected_graph_from_edge_list()
		self.undirected_healthy= self.create_undirected_graph_from_edge_list()

		nx.set_node_attributes(self.undirected_cancer , self.cancer_dict)
		nx.set_node_attributes(self.undirected_healthy, self.healthy_dict)

	def read_csv_file(self,entry):
		df = pd.read_csv(entry)
		all_rows = [df.columns.values.tolist()]+df.values.tolist()
		return all_rows

	def dictionary_maker(self,data):
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

	def populate_attributes(self,dict,fail):
		master = {}
		for i,attribute_dict in dict.items():
			for attribute_key,attribute in attribute_dict.items():
				if attribute_key != 'int_id':
					try:
						if attribute_key not in master.keys():
							master[attribute_key] = float(attribute)
						else:
							master[attribute_key] += float(attribute)
					except:
						fail.add(attribute_key)
						pass
		return master

	def create_undirected_graph_from_edge_list(self):
		graph = nx.Graph()
		edges = []
		for row in self.edge_list[1:]:
			edges.append([int(row[0]),int(row[1])])
		graph.add_edges_from(edges)
		return graph

class graph_analysis(object):
	def __init__(self,data_obj):
		self.data_obj = data_obj

		self.undirected_cancer = create_undirected_graph_from_edge_list()
		self.undirected_healthy= create_undirected_graph_from_edge_list()

		nx.set_node_attributes(self.undirected_cancer , self.data_obj.cancer_dict)
		nx.set_node_attributes(self.undirected_healthy, self.data_obj.healthy_dict)

	def give_level_1(self,input_graph,node):
		neighbors = input_graph.neighbors(node)
		return neighbors

	def give_level_2(self,input_graph,node):
		neighbors = input_graph.neighbors(node)
		return neighbors


	def create_undirected_graph_from_edge_list(self):
		graph = nx.Graph()
		edges = []
		for row in self.edge_list[1:]:
			edges.append([int(row[0]),int(row[1])])
		graph.add_edges_from(edges)
		return graph

	def get_uniform_attribute(self,input_dict, sample_list):
		uniform_att = {}
		for key,item in input_dict.items():
			for attribute in item:
				if attribute == 'int_id' and item[attribute] in sample_list:
					continue
				else:
					if math.isnan(float(item[attribute])):
						if attribute not in uniform_att.keys():
							try:
								uniform_att[attribute] = [float(item[attribute])]
							except:
								uniform_att[attribute].append(float(item[attribute]))

		for key, item in uniform_att.items():
			uniform_att[key] = np.mean(item)

		return uniform_att

	def prediction(self,input_graph,node_to_predict,level,cancer,sample):
		neighbors = []
		if(level == 1):
			neighbors = self.give_level_1(input_graph,node_to_predict)
		elif(level == 2):
			direct_neighbors = self.give_level_1(input_graph,node_to_predict)
			for node in direct_neighbors:
				x = self.give_level_2(input_graph,node,node_to_predict)
				for i in x:
					neighbors.append(i)
		if(neighbors == []):
			if(cancer == True):
				new = self.get_uniform_attributes_cancer(cancer_dict,sample)
				new['int_id'] = node
				return new
			else:
				new = self.get_uniform_attributes_healthy(healthy_dict,sample)
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


main = DATA_OBJECT('cancer_data.csv', 'healthy_data.csv', 'id_uniprot_ensembl_map.csv', 'edge_list.csv')