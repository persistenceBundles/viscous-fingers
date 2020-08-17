import os, time
import itertools, sys, pickle, inspect
from helper import *
from scipy import spatial

# Remove the try-catch once testing is complete
try:
	__file__
except:
	sys.argv = [sys.argv[0], 'tv_102', 'tv_103']

colors = ['white', 'khaki', 'seagreen', 'IndianRed', 'mediumpurple', 'darkorange', 'yellowgreen', 'gold', 'orchid2', 'PeachPuff', 'skyblue', 'coral', 'plum', 'darkolivegreen', 'crimson', 'DarkGoldenrod3', 'mediumvioletred', 'sienna3', 'cyan', 'darkseagreen', 'rosybrown', 'honeydew']

def traverse(dictionary, indices):
	indices = printIndices(dictionary, indices)
	if indices is None:
		return
	if len(indices) == 2:
		print ('The matrix has been opened :|')
	if (indices[0] == S_MATRIX_IDENTIFIER):
		traverse(S,indices)
	elif(indices[0] == S1_MATRIX_IDENTIFIER):
		traverse(S1, indices)
	elif(indices[0] == S2_MATRIX_IDENTIFIER):
		traverse(S2, indices)
	else:
		print ('Error!')

def get_style(index):
	if index == GAP_NODE:
		style_attribute = "shape=diamond style=filled fillcolor=lightslategray"
	else:
		 style_attribute = "shape=circle style=filled fillcolor=" + colors[index]
	return style_attribute

def get_node(index, pairs, mappings, labels, costs, color_index):
	label_attribute = get_label(index, pairs, mappings, labels)
	style_attribute = get_style(color_index)
	return str(mappings[index]) + " [" + label_attribute  + " " + style_attribute + "]\n"

def print_tree():
	file_path = os.path.abspath(inspect.getfile(inspect.currentframe()))
	filenames = join_strings([filename1, filename2])
	compare_path = get_output_path(file_path, [COMPARE_PREFIX, filenames, DOT_EXTENSION], folder_name=COMPARE_GRAPHS_FOLDER)
	compare_file = open(compare_path, 'w')

	# print header for first tree
	compare_file.write("digraph T1 {\n")

	# print nodes dictionary for first tree
	for i in range(1, right1[1]+1):
		if map1[i] != GAP_NODE:
			compare_file.write(get_node(i, pairs1, index_mapping1, label1, cost1, i))
		else:
			compare_file.write(get_node(i, pairs1, index_mapping1, label1, cost1, GAP_NODE))

	# print nodes connectivity for first tree
	for i in range(1, right1[1]+1):
		compare_file.write(get_connectivity(i, parent1, index_mapping1))

	compare_file.write("}\n\n")

	# print header for second tree
	compare_file.write("digraph T2 {\n")

	# print nodes dictionary for second tree
	for j in range(1, right2[1]+1):
		compare_file.write(get_node(j, pairs2, index_mapping2, label2, cost2, map2[j]))

	# print nodes connectivity for second tree
	for j in range(1, right2[1]+1):
		compare_file.write(get_connectivity(j, parent2, index_mapping2))

	compare_file.write("}\n")

	compare_file.close()

	# get comparison image
	run_shell_script(MAKE_IMAGE_SCRIPT, [compare_path])

def printIndices(dictionary, k):
	try:
		try:
			# parse values from array
			i = k[2]
			j = k[4]
			operation = k[5]
			message = k[6]
			cost = round(k[7],3)

			#print (i, j, operation, message, cost)
			#print (message)

			# according to each operation store associated mapping and cost
			if operation in [RELABEL_IDENTIFIER]:
				# i+1 is relabelled to j+1
				map1[i+1] = j+1
				map2[j+1] = i+1
				cost1[i+1] = cost
				cost2[j+1] = 0
			elif operation in [T1_STARTING_GAP_IDENTIFIER, T1_CONTINUING_GAP_IDENTIFIER]:
				# i+1 is a gap in T1
				map1[i+1] = GAP_NODE
				# Adhitya may be wrong here. Ask Vijay
				cost1[i+1] = cost * -1
			elif operation in [T1_GENERIC_GAP_IDENTIFIER]:
				# i is a gap in T1
				map1[i] = GAP_NODE
				cost1[i] = cost
			elif operation in [T2_STARTING_GAP_IDENTIFIER, T2_CONTINUING_GAP_IDENTIFIER]:
				# j+1 is a gap in T2
				map2[j+1] = GAP_NODE
				cost2[j+1] = cost
			elif operation == T2_GENERIC_GAP_IDENTIFIER:
				# j is a gap in T2
				map2[j] = GAP_NODE
				cost2[j] = cost
			else:
				print ('Govinda')

		except IndexError:
			# This should only happen only during the first call of this function
			pass

		next_comparison = dictionary[k[1]][k[2]][k[3]][k[4]]

		return next_comparison
	except:
		#print ('Done! :)')
		return None


# Get the paths
file_path = os.path.abspath(inspect.getfile(inspect.currentframe()))

# Take the names of files as arguments
filename1 =  sys.argv[1]
filename2 = sys.argv[2]

try:
	# Get right-most node for each node
	right1 = get_dictionary(file_path, [filename1, RIGHT_NODE_SUFFIX])
	right2 = get_dictionary(file_path, [filename2, RIGHT_NODE_SUFFIX])

	# Get parent of each node
	parent1 = get_dictionary(file_path, [filename1, PARENT_NODE_SUFFIX])
	parent2 = get_dictionary(file_path, [filename2, PARENT_NODE_SUFFIX])

	# Get function value of each node
	label1 = get_dictionary(file_path, [filename1, LABEL_NODE_SUFFIX])
	label2 = get_dictionary(file_path, [filename2, LABEL_NODE_SUFFIX])

	# Get persistence of each node
	difference1 = get_dictionary(file_path, [filename1, DIFFERENCE_NODE_SUFFIX])
	difference2 = get_dictionary(file_path, [filename2, DIFFERENCE_NODE_SUFFIX])

	pairs1 = get_dictionary(file_path, [filename1, PAIRS_NODE_SUFFIX])
	pairs2 = get_dictionary(file_path, [filename2, PAIRS_NODE_SUFFIX])

	index_mapping1 = get_dictionary(file_path, [filename1, MAPPING_NODE_SUFFIX])
	index_mapping2 = get_dictionary(file_path, [filename2, MAPPING_NODE_SUFFIX])

	# Get size of both the trees
	size1 = len(parent1.keys())
	size2 = len(parent2.keys())

	# How much are we planning to compare?
	extents = [0, 1, right1[1], 1, right2[1]]

	# Get the intermediate matrices
	Q = get_matrix(file_path, [filename1, filename2, Q_IDENTIFIER])
	Q1 = get_matrix(file_path, [filename1, filename2, Q1_IDENTIFIER])
	Q2 = get_matrix(file_path, [filename1, filename2, Q2_IDENTIFIER])

	S = get_matrix(file_path, [filename1, filename2, S_IDENTIFIER])
	S1 = get_matrix(file_path, [filename1, filename2, S1_IDENTIFIER])
	S2 = get_matrix(file_path, [filename1, filename2, S2_IDENTIFIER])

	# store mappings from one tree to another
	map1 = {}
	map2 = {}

	# store costs from one tree to another
	cost1 = {}
	cost2 = {}

	# mapping from vertex indice to preorder/map index
	inverse1 = {v: k for k, v in index_mapping1.iteritems()}
	inverse2 = {v: k for k, v in index_mapping2.iteritems()}

	# store list of vertex indices in each tree
	vertex_indices1 = inverse1.keys()
	vertex_indices2 = inverse2.keys()

except:
	print ("Something bad happened :(", filename1, filename2)

print ("*******")

print (filename1, filename2)

# Start tracing back from the back
traverse(S, extents)
print_tree()

# get coordinates for each node using vertex index
coords1 = [find_coords(vertex_index) for vertex_index in [index_mapping1[i] for i in map1.keys()]]
coords2 = [find_coords(vertex_index) for vertex_index in [index_mapping2[j] for j in map2.keys()]]

# data structure to hold the coordinates
spatial_tree1 = spatial.KDTree(coords1)
spatial_tree2 = spatial.KDTree(coords2)

# iterate through mappings of first tree
for i in map1.keys():
	# get the coordinate which closest to this node in the second tree
	distance, ideal_index = spatial_tree2.query([coords1[i-1]])
	ideal_index = ideal_index[0]+1
	# get vertex index of ideal mapping
	ideal_mapping = index_mapping2[ideal_index]
	# get preorder index of mapping in second tree
	j = map1[i]
	# get vertex index of node
	index1 = index_mapping1[i]
	# if node is not a gap
	if j != GAP_NODE:
		actual_mapping = index_mapping2[j]
		if ideal_mapping == actual_mapping:
			#continue
			print (i, index1, j, actual_mapping, 'cool :)')
		else:
			print (i, index1, j, actual_mapping, ideal_index, ideal_mapping,'relabel?')
	else:
		print (i, index1, ideal_index, ideal_mapping, 'gap?')

# same as above but only iterate for gaps in second tree
for j in map2.keys():
	i = map2[j]
	if i == GAP_NODE:
		# do fancy checking only if node is a gap

		# get the coordinate which closest to this node in the first tree
		distance, ideal_index = spatial_tree1.query([coords2[j-1]])
		ideal_index = ideal_index[0]+1

		# get vertex index of ideal mapping and current node
		ideal_mapping = index_mapping1[ideal_index]
		index2 = index_mapping2[j]

		print (j, index2, ideal_index, ideal_mapping, 'gap?')

print ("*******")

# Write the Merge Tree to file
debug_file_arguments = [join_strings([filename1,filename2]), CSV_EXTENSION]
debug_file_path = get_output_path(file_path, debug_file_arguments, folder_name = DEBUG_FOLDER)

debug_file = open(debug_file_path, 'w')
fieldnames = ['tree', 'vertex', 'index', 'parent', 'scalar', 'pair', 'persistence', 'match']
writer = csv.writer(debug_file, delimiter=',')
writer.writerow(fieldnames)

for i in map1.keys():
	content = ['T1', i, index_mapping1[i], parent1[i], round(label1[i],3), pairs1[i], round(difference1[i],3), map1[i]]
	writer.writerow(content)
for j in map2.keys():
	content = ['T2', j, index_mapping2[j], parent2[j], round(label2[j],3), pairs2[j], round(difference2[j],3), map2[j]]
	writer.writerow(content)
