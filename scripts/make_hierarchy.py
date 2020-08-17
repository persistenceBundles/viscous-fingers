from helper import *
from functools import cmp_to_key
import json

tree_type = TREE_TYPE_SPLIT

file_name = ''
file_path = ''

split_scalars = {}
visited = {}
adjacency = {}
pairs = {}

x_coordinates = {}
y_coordinates = {}
z_coordinates = {}

index_map = {}
postorder_map = {}
preorder_map = {}

birth = {}
death = {}

ratio_threshold = 0.0

class Tree(object):
	def __init__(self):
		self.index = None
		self.children = []
		self.parent = None
		self.label = None
		self.pair = None
		self.birth = None
		self.death = None
		self.preorder = None
		self.postorder = None
		self.merged_nodes = []
		self.merged = False
		self.x = None
		self.y = None
		self.z = None

	def __str__(self):
		values = [self.preorder, self.pair.preorder, self.birth.preorder, self.death.preorder]
		return ' '.join(map(str, values))

	def is_root(self):
		return (self.parent is None)

	def is_leaf(self):
		return (len(self.children) == 0)

	def is_saddle(self):
		return ((not(self.is_root())) and (not(self.is_leaf())))

	def merge_node(self, parent):
		# sign papers for adoption
		for child in self.children:
			child.parent = parent
			child.pair = parent
			child.birth = child
			child.death = parent
		# transfer children of node to parent
		parent.children.extend(self.children)
		# transfer merged nodes of node to parent
		parent.merged_nodes.extend(self.merged_nodes)
		# remove children of current node after transfer
		self.children = []
		# make node to be merged with parent
		parent.merged_nodes.append(self)
		# remove node from children of parent
		parent.children.remove(self)
		# mark node as merged
		# the way to access other outer nodes would be using node.parent
		self.merged = True



class Hierarchy(object):
	def __init__(self, node):
		self.index = str(node.index)
		self.children = []
		self.birth = node.label
		self.death = node.pair.label
		self.persistence = round(float(self.death - self.birth), 5)


class HierarchyEncoder(json.JSONEncoder):
    def default(self, obj):
        if not isinstance(obj, Hierarchy):
            return super(HierarchyEncoder, self).default(obj)

        return obj.__dict__

def initialize_tree(index):
	root = Tree()
	root.index = index
	root.label = split_scalars[index]
	root.pair = pairs[index]

	# add coordinates as well
	root.x = x_coordinates[index]
	root.y = y_coordinates[index]
	root.z = z_coordinates[index]

	# add mapping to dictionary
	index_map[index] = root

	return root

def add_node(index, parent):
	node = Tree()
	node.index = index
	parent.children.append(node)
	node.parent = parent
	node.label = split_scalars[index]
	node.pair = pairs[index]

	# add coordinates as well
	node.x = x_coordinates[index]
	node.y = y_coordinates[index]
	node.z = z_coordinates[index]

	# add mapping to dictionary
	index_map[index] = node

	return node


def compare_vertices(a, b):
	# try to sort using the split_scalars
	# if they are equal, sort using index value
	if split_scalars[a] > split_scalars[b]:
		return 1
	elif split_scalars[a] == split_scalars[b]:
		if a > b:
			return 1
		else:
			return -1
	else:
		return -1

# same function as above, but this time compare the tree nodes
def compare_nodes(a, b):
	# try to sort using the split_scalars
	# if they are equal, sort using index value
	if a.label > b.label:
		return 1
	elif a.label == b.label:
		if a.index > b.index:
			return 1
		else:
			return -1
	else:
		return -1

def traverse(index, parent):
	#print (index, split_scalars[index])
	visited[index] = True

	# cmp_to_key is for Python 3
	sorted(adjacency[index], key=cmp_to_key(compare_vertices))
	for node in adjacency[index]:
		if not visited[node]:
			current = add_node(node, parent)
			traverse(node, current)

def add_pairs(node):
	if node is None:
		return
	else:
		node.pair = index_map[pairs[node.index]]
		node.birth = index_map[birth[node.index]]
		node.death = index_map[death[node.index]]
		for child in node.children:
			add_pairs(child)

def postorder(node):
	# python needs a mutable object for updation
	order = {'index': 1}

	def set_order(node):
		if node is None:
			return
		else:
			for child in node.children:
				set_order(child)

			node.postorder = order['index']
			postorder_map[order['index']] = node
			order['index'] += 1

	set_order(node)

def preorder(node):
	# python needs a mutable object for updation
	order = {'index': 1}

	def set_order(node):
		if(node == None):
			return
		else:
			node.preorder = order['index']
			preorder_map[order['index']] = node
			order['index'] += 1

			for child in node.children:
				set_order(child)

	set_order(node)


# return the height for a given node
def height(node):
	if node is None:
		return 0

	max_depth = 0
	for child in node.children:
		max_depth = max(max_depth, height(child))

	return max_depth + 1

# Function to traverse level order traversal of tree
def merge_unstable_saddles(node):
	h = height(node)
	# traverse from below
	for i in range(h+1, 0, -1):
		traverse_level(node, i)
		#print ('')

# find the ratio from bottom-up
def traverse_level(node, level):
	if node is None:
		return
	if level == 1:
		# do this merging only for saddles
		if node.is_saddle():
			# lets not tinker with the root
			if not node.parent.is_root():
				parent = node.parent
				ratio =  (node.label - parent.label)

				# in this case I'm merging every saddle out there
				if ratio < ratio_threshold:
					node.merge_node(parent)
					#print (round(ratio,4), node.label, parent.label, 'merge')
				else:
					pass
					#print (node.label, ratio, 'chill :)')
	elif level > 1:
		# traverse backwards
		for child in reversed(node.children):
			traverse_level(child, level - 1)

# return the right most node for a given node
def right_leaf(node):
	if node is None:
		return
	if not node.children:
		return node
	else:
		# we consider the last child as right most
		return node.children[-1]

# sort the children at each level
def sort_tree(node):
	if node is None:
		return
	else:
		if not node.is_leaf():
			sorted(node.children, key=cmp_to_key(compare_nodes))
			for child in node.children:
				sort_tree(child)


def get_merge_tree():
	# Get merge tree path
	tree_file_arguments = [tree_type, TREE_INFIX, file_name, CSV_EXTENSION]
	tree_file_path = get_output_path(file_path, tree_file_arguments, folder_name = TREES_FOLDER)

	# Read merge tree file
	with open(tree_file_path, 'r') as csvfile:
		csvfile.readline()
		spamreader = csv.reader(csvfile, delimiter=' ')
		for r in spamreader:
			row = r[0].split(',')
			node1 = int(row[0])
			node2 = int(row[1])

			split_scalars[node1] = float(row[2])
			split_scalars[node2] = float(row[3])

			visited[node1] = False
			visited[node2] = False

			if node1 not in adjacency.keys():
				adjacency[node1] = []

			if node2 not in adjacency.keys():
				adjacency[node2] = []

			adjacency[node1].append(node2)
			adjacency[node2].append(node1)

	for i in adjacency.keys():
		if len(adjacency[i]) == 1:
			if (split_scalars[i] < split_scalars[adjacency[i][0]]):
				root = i

	return root

def get_persistent_pairs():
	# Get persistence pairs
	pairs_file_arguments = [tree_type, PAIRS_INFIX, file_name, CSV_EXTENSION]
	pairs_file_path = get_output_path(file_path, pairs_file_arguments, folder_name = PAIRS_FOLDER)

	with open(pairs_file_path, 'r') as persistence_pairs:
		persistence_pairs.readline()
		spamreader = csv.reader(persistence_pairs, delimiter=' ')
		for r in spamreader:
			row = r[0].split(',')
			node1 = int(row[0])
			node2 = int(row[1])

			#if (node1 in split_scalars.keys()) and (node2 in split_scalars.keys()):
			# there will be pairs that do not exist in the merge tree
			# they will be removed/ignored subsequently

			pairs[node1] = node2
			pairs[node2] = node1

			# add birth and death values of nodes to dictionaries
			birth[node1] = node1
			death[node1] = node2

			birth[node2] = node1
			death[node2] = node2


def get_coordinates():
	coordinates_file_arguments = [tree_type, COORDINATES_INFIX, file_name, CSV_EXTENSION]
	coordinates_file_path = get_output_path(file_path, coordinates_file_arguments, folder_name = COORDINATES_FOLDER)

	with open(coordinates_file_path, 'r') as coordinates_file:
		coordinates_file.readline()
		spamreader = csv.reader(coordinates_file, delimiter=' ')
		for r in spamreader:
			row = r[0].split(',')
			identifier = int(row[0])
			x = round(float(row[1]),3)
			y = round(float(row[2]),3)
			z = round(float(row[3]),3)

			x_coordinates[identifier] = x
			y_coordinates[identifier] = y
			z_coordinates[identifier] = z

def print_tree(node):
	if node is None:
		return
	else:
		print (node)
		for child in node.children:
			print_tree(child)

def write_postorder(node):
	postorder_file_arguments = [tree_type, POSTORDER_INFIX, file_name, CSV_EXTENSION]
	postorder_file_path = get_output_path(file_path, postorder_file_arguments, folder_name = POSTORDER_FOLDER)

	postorder_file = open(postorder_file_path, 'w')
	fieldnames = ['node', 'order']

	writer = csv.writer(postorder_file, delimiter=',')
	writer.writerow(fieldnames)

	for order_index in postorder_map.keys():
		content = [postorder_map[order_index].index, order_index]
		writer.writerow(content)

	postorder_file.close()


def print_tree_label(node):
	print (str(node.preorder) + " [label=\""+ str(node.preorder) + " \\n["+ str(node.pair.preorder) + "] \\n" +str(node.label) +"\"]")

def print_tree_edge(node):
	print (str(node.parent.preorder) + "->" + str(node.preorder))

def print_tree_dot(node):
	if node is None:
		return
	else:
		print_tree_label(node)
		for child in node.children:
			print_tree_edge(child)
			print_tree_dot(child)

def save_dictionaries(node):
	right_dictionary = {}
	parent_dictionary = {}
	label_dictionary = {}
	difference_dictionary = {}
	pairs_dictionary = {}
	inverse_index_map = {}

	def populate_dictionaries(node):
		if node is None:
			return
		else:
			right_dictionary[node.preorder] = (right_leaf(node)).preorder
			label_dictionary[node.preorder] = node.label
			if node.parent is None:
				parent_dictionary[node.preorder] = 0
			else:
				parent_dictionary[node.preorder] = node.parent.preorder

			difference_dictionary[node.preorder] = abs(node.label - node.pair.label)
			pairs_dictionary[node.preorder] = node.pair.preorder
			inverse_index_map[node.index] = node.preorder

			#print node.index, node.preorder, tree.label, tree.parent.preorder

			for child in node.children:
				populate_dictionaries(child)

	populate_dictionaries(node)

	# save dictionaries to respective files
	save_dictionary(right_dictionary, file_path, file_name, RIGHT_NODE_SUFFIX)
	save_dictionary(parent_dictionary, file_path, file_name, PARENT_NODE_SUFFIX)
	save_dictionary(label_dictionary, file_path, file_name, LABEL_NODE_SUFFIX)
	save_dictionary(difference_dictionary, file_path, file_name, DIFFERENCE_NODE_SUFFIX)
	save_dictionary(pairs_dictionary, file_path, file_name, PAIRS_NODE_SUFFIX)
	save_dictionary(inverse_index_map, file_path, file_name, MAPPING_NODE_SUFFIX)

def write_recon_file(node):
	recon_file_arguments = [file_name, TXT_EXTENSION]
	recon_file_path = get_output_path(file_path, recon_file_arguments, folder_name = RECON_FOLDER)

	recon_file = open(recon_file_path, 'w')

	writer = csv.writer(recon_file, delimiter=',')
	writer.writerow([len(postorder_map.keys())])

	def recon_print_tree(node):
		if(node == None):
			return
		else:
			if node.parent == None:
				parent_order = -1
			else:
				parent_order = node.parent.postorder

			values = [node.postorder, parent_order, node.death.postorder, -1, node.label]
			num_children = len(node.children)
			values.append(num_children)
			if num_children != 0:
				for child in node.children:
					values.append(child.postorder)

			writer.writerow(values)

			for child in node.children:
				recon_print_tree(child)

	recon_print_tree(node)
	recon_file.close()

def save_tree_string(tree):
	def stringify_tree(node):
		tree_string = ""

		def make_tree_string(node):
			nonlocal tree_string
			if(node == None):
				return
			else:
				tree_string += '{'
				tree_string += str(node.postorder) + '|'
				tree_string += str(node.index) + '|'
				tree_string += str(node.label) + '|'
				tree_string += str(node.birth.label) + '|'
				tree_string += str(node.death.label)

				for child in node.children:
					make_tree_string(child)

				tree_string += '}'

		make_tree_string(node)
		return tree_string

	string_file_arguments = [file_name, TXT_EXTENSION]
	string_file_path = get_output_path(file_path, string_file_arguments, folder_name = STRINGS_FOLDER)
	save_string(stringify_tree(tree), string_file_path)

def get_persistence_hierarchy(node):
	processed_nodes = {}
	hierarchy_map = {}

	hierarchy_root = Hierarchy(node)
	hierarchy_map[node.preorder] = hierarchy_root

	def find_hierarchy_parent(node):
		if node.preorder not in processed_nodes:
			parent_node = node.parent
			parent_found = False

			# keep going up the merge tree till you find a parent that itself and its pair within the range
			while((parent_node != None) and (parent_found == False)):
				if parent_node.preorder < node.preorder < parent_node.pair.preorder:
					parent_found = True
				else:
					parent_node = parent_node.parent

			if not parent_found:
				pass
			else:
				hierarchy_node = Hierarchy(node)
				hierarchy_map[node.preorder] = hierarchy_node

				hierarchy_parent = hierarchy_map[parent_node.preorder]
				#hierarchy_node.parent = hierarchy_parent
				hierarchy_parent.children.append(hierarchy_node)

			processed_nodes[node.preorder] = True
			processed_nodes[node.pair.preorder] = True

	def get_hierarchy_parent(node):
		if node is None:
			return
		else:
			find_hierarchy_parent(node)
			for child in node.children:
				get_hierarchy_parent(child)

	get_hierarchy_parent(node)
	return hierarchy_root

def print_treemap(node):
	processed_nodes = {}
	treemap_string = {}
	treemap_value = {}
	treemap_parent = {}
	treemap_container = {}

	def find_treemap_parent(node):
		if node.preorder not in processed_nodes:
			parent_node = node.parent
			paired_node = node.pair
			parent_found = False

			# keep going up the merge tree till you find a parent that itself and its pair within the range
			while((parent_node != None) and (parent_found == False)):
				if parent_node.preorder < node.preorder < parent_node.pair.preorder:
					parent_found = True
				else:
					parent_node = parent_node.parent

			if not parent_found:
				treemap_container[node.preorder] = str(node.preorder)
				treemap_parent[node] = None
				treemap_parent[node.pair] = node
			else:
				treemap_container[node.preorder] = treemap_container[parent_node.preorder] + "." + str(node.preorder)
				treemap_parent[node.pair] = node
				treemap_parent[node] = parent_node

			treemap_string[node.preorder] = treemap_container[node.preorder] + "." + str(node.preorder)
			treemap_string[node.pair.preorder] = treemap_container[node.preorder] + "." + str(node.pair.preorder)

			treemap_value[node.pair.preorder] = node.pair.label
			treemap_value[node.preorder] = node.label

			processed_nodes[node.preorder] = True
			processed_nodes[node.pair.preorder] = True

	def get_tree_structure(node):
		if(node == None):
			return
		else:
			find_treemap_parent(node)
			for child in node.children:
				get_tree_structure(child)

	get_tree_structure(node)
	for key in treemap_container.keys():
		print str(treemap_container[key]) + ","

	for key in treemap_string.keys():
		print str(treemap_string[key]) + ","+ str(int((treemap_value[key]+0.05)*1000))

def save_hierarchy(hierarchy_root):
	hierarchy_file_arguments = [file_name, JSON_EXTENSION]
	hierarchy_file_path = get_output_path(file_path, hierarchy_file_arguments, folder_name = HIERARCHY_FOLDER)

	hierarchy_file = open(hierarchy_file_path, 'w')
	json.dump(hierarchy_root, hierarchy_file, cls=HierarchyEncoder, indent=1, sort_keys=True)
	hierarchy_file.close()

def make_hierarchy(name, path):
	global file_name, file_path
	file_name = name
	file_path = path
	root = get_merge_tree()
	get_persistent_pairs()
	get_coordinates()

	tree = initialize_tree(root)
	traverse(root, tree)
	add_pairs(tree)

	# merge saddles -- use the ratio_threshold at the top
	merge_unstable_saddles(tree)
	sort_tree(tree)

	# assign postorder indices after merging of saddles
	# postorder is for making the thresholding easier after segmentation in Paraview
	# preorder is used for making the hierarchy
	postorder(tree)
	preorder(tree)

	# just write all the postorder indices to a file
	write_postorder(tree)

	# write the dictionaries for Xu edit distance comparison
	save_dictionaries(tree)

	# write the file for ordered APTED edit distance comparison
	save_tree_string(tree)

	# write the file for unordered constrained edit distance comparison
	write_recon_file(tree)

	# print to console for quick checking on webgraphviz
	print ('digraph {')
	print_tree_dot(tree)
	print ('}')

	print("")

	hierarchy_root = get_persistence_hierarchy(tree)
	save_hierarchy(hierarchy_root)
	

# if you make changes to this file, make sure to stash it to Paraview's directory
