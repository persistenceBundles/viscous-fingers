import os, re, shutil, pickle, inspect, csv, sys, math, io
import numpy as np
from collections import defaultdict
from collections import Counter
from scipy.spatial import distance

# List of constants
CSV_EXTENSION = '.csv'
TXT_EXTENSION = '.txt'
BIN_EXTENSION = '.bin'
VTP_EXTENSION = '.vtp'
PNG_EXTENSION = '.png'
DOT_EXTENSION = '.dot'
JT_EXTENSION = '.jt'
VTK_EXTENSION = '.vtk'
JSON_EXTENSION = '.json'

PAIRS_INFIX = '-pairs-'
TREE_INFIX = '-tree-'
NODES_INFIX = '-nodes-'
ARCS_INFIX = '-arcs-'
MATRIX_INFIX = '-matrix-'
SEGMENTATION_INFIX = '-segmentation-'
SCREENSHOT_INFIX = '-screenshot-'
PARENTS_INFIX = '-parents-'
POSTORDER_INFIX = '-postorder-'
TRIANGLES_INFIX = '-triangles-'
COMPARE_PREFIX = 'compare-'
EDIT_DISTANCE_RESULT = 'results'

RIGHT_NODE_SUFFIX = '-right'
PARENT_NODE_SUFFIX = '-parent'
LABEL_NODE_SUFFIX = '-labels'
DIFFERENCE_NODE_SUFFIX = '-difference'
PAIRS_NODE_SUFFIX = '-pairs'
MAPPING_NODE_SUFFIX = '-mapping'
BIRTH_NODE_SUFFIX = '-birth'
DEATH_NODE_SUFFIX = '-death'
CONTOUR_TREE_SUFFIX = '-contour'
SPLIT_TREE_SUFFIX = '-split'
STRING_SUFFIX = '-string'
TRIANGLES_UNMAPPED_INFIX = '-triangles-unmapped-'
SEGMENTATION_TRIANGLES_PREFIX = 'segmentation-triangles-'
COMPARISON_APTED_INFIX = '-comparison-apted-'
COMPARISON_UNORDERED_INFIX = '-comparison-unordered-'
COORDINATES_INFIX = '-coordinates-'
CONFIDENCE_INFIX = '-confindence-'

TREE_TYPE_SPLIT = 'split'
TREE_TYPE_CONTOUR = 'contour'

Q_IDENTIFIER = '-Q'
Q1_IDENTIFIER = '-Q1'
Q2_IDENTIFIER = '-Q2'

S_IDENTIFIER = '-S'
S1_IDENTIFIER = '-S1'
S2_IDENTIFIER = '-S2'

INPUT_FOLDER = 'input'
OUTPUT_FOLDER = 'output'
PAIRS_FOLDER = 'pairs'
TREES_FOLDER = 'trees'
POSTORDER_FOLDER = 'postorder'
TRIANGLES_FOLDER = 'triangles'
TRIANGLES_UNMAPPED_FOLDER = 'triangles-unmapped'
SCRIPTS_FOLDER = 'scripts'
SCREENSHOT_FOLDER = 'screenshots'
SEGMENTATION_FOLDER = 'segmentation'
INTERMEDIATE_FOLDER = 'intermediate'
PERSISTENCE_FOLDER = 'persistence'
DICTIONARY_FOLDER = 'dictionary'
RESULTS_FOLDER = 'results'
MATRICES_FOLDER = 'matrices'
IMAGES_FOLDER = 'images'
GRAPHS_FOLDER = 'graph'
COMPARE_GRAPHS_FOLDER = 'compare-graphs'
COMPARE_IMAGES_FOLDER = 'compare-images'
MERGED_GRAPHS_FOLDER = 'merged-graphs'
MERGED_IMAGES_FOLDER = 'merged-images'
DEBUG_FOLDER = 'debug'
JT_FOLDER = 'jt'
CLIQUE_GRAPHS_FOLDER = 'clique-graphs'
STRINGS_FOLDER = 'strings'
STABILITY_FOLDER = 'stability'
PARENTS_FOLDER = 'parents'
MATRIX_FOLDER = 'matrix'
SEGMENTATION_TRIANGLES_FOLDER = 'segmentation-triangles'
COMPARISON_APTED_FOLDER = 'comparison-apted'
COMPARISON_UNORDERED_FOLDER = 'comparison-unordered'
COORDINATES_FOLDER = 'coordinates'
TUPLES_FOLDER = 'tuples'
RECON_FOLDER = 'recon'
CONFIDENCE_FOLDER = 'confidence' # what a name :P
HIERARCHY_FOLDER = 'hierarchy'

PYTHON_COMMAND = 'python'
PARAVIEW_COMMAND = 'paraview'
JAR_COMMAND = 'java -jar'

COMPUTE_SCRIPT = 'compute.py'
SPLIT_MAKE_GRAPH_SCRIPT = 'split-make-graph.py'
SPLIT_MAKE_GRAPH_LEFT_SCRIPT = 'split-make-graph-left.py'
MAKE_IMAGE_SCRIPT = 'make-image.sh'
MAKE_STABLE_SCRIPT = 'make-stable.py'
GENERATE_JT_FILES_SCRIPT = 'generate-jt-files.py'

CONTOUR_MAKE_TREE_SCRIPT = 'contour-make-tree.py'
SPLIT_MAKE_TREE_SCRIPT = 'split-make-tree.py'

TIMESTEP_PREFIX = 'tv_'

INFINITY = float('inf')

UNKNOWN_COST = float('-inf')
RELABEL_IDENTIFIER = 0

T1_STARTING_GAP_IDENTIFIER = 1
T1_CONTINUING_GAP_IDENTIFIER = 2
T1_GENERIC_GAP_IDENTIFIER = 3
T1_RIGHT_GAP_IDENTIFIER = 4

T2_STARTING_GAP_IDENTIFIER = 5
T2_CONTINUING_GAP_IDENTIFIER = 6
T2_GENERIC_GAP_IDENTIFIER = 7
T2_RIGHT_GAP_IDENTIFIER = 8

S_MATRIX_IDENTIFIER = 0
S1_MATRIX_IDENTIFIER = 1
S2_MATRIX_IDENTIFIER = 2

GAP_NODE = -1

# in addition to the normal isoband we have few special ones in case equality
# the following cases will not have an isovalue in-between
# indices are as [<lower, =lower, <lower, upper>, =upper, >upper]
# which are translated as [0, 1, 2, 3, 4] respectively
#
# 0 - 1, 1 - 0, 1 - 1, 3 - 3, 3 - 4 are not possible now
# 0 - 0, 4 - 4 were not possible even without equality check

lower_intolerable_indices = ['000', '001', '010', '011', '100', '101', '110', '111']
upper_intolerable_indices = ['333', '334', '343', '344', '433', '434', '443', '444']

cheap_intolerable_indices = ['000', '222']

# get working directory and add '/' at the end
def cwd():
	return os.path.join(os.getcwd(), '')

# get path of current file
def current_path():
	return os.path.abspath(inspect.getfile(inspect.currentframe()))

# Replace a pattern with another in a file
def replace_wildcard(fname, pat, s_after):
    # first, see if the pattern is even in the file.
    with open(fname) as f:
        if not any(re.search(pat, line) for line in f):
            return # pattern does not occur in file so we are done.

    # pattern is in the file, so perform replace operation.
    with open(fname) as f:
        out_fname = fname + ".tmp"
        out = open(out_fname, "w")
        for line in f:
            out.write(re.sub(pat, s_after, line))
        out.close()
        os.rename(out_fname, fname)

# Return the just the name when path is False
# Return the name appended with the parent path when path is True
def get_file_name(file_path, path = False):
	file_name = os.path.splitext(os.path.basename(file_path))[0]
	parent_path = get_parent_path(file_path)
	if path:
		return join_file_path(parent_path, file_name)
	else:
		return file_name

def get_file_extension(file_path):
	file_basename = os.path.basename(file_path)
	file_text = os.path.splitext(file_basename)
	return file_text[1]

def get_parent_path(file_path):
	return os.path.abspath(os.path.join(file_path, os.pardir))

def get_folder(name):
	return get_output_path(cwd(), [], folder_name=name)

def get_input_path(file_path):
	return os.path.join(get_parent_path(file_path), INPUT_FOLDER)

def get_triangles_path(file_path):
	return os.path.join(get_parent_path(file_path), TRIANGLES_FOLDER)

# Takes contourForest.TreeType parameter and returns a string
def get_tree_type(tree_type):
	return tree_type.split(' ')[0].lower()

# join two strings
def join_strings(strings):
	return '-'.join(strings)

def get_timestep(string):
	# remove extension if it exists
	file_name = string.split('.')[0]
	timestep = file_name.split(TIMESTEP_PREFIX)[1]
	return timestep

# Get a new filename in the output directory
# This takes in a file from input directory and gives out a string with ../output/file_name
# arguments can be sent in a list to the aforementioned string
# output_folder can be set to False for ../file_name
# pass the folder name if you don't want the default output_folder

def get_output_path(file_path, arguments, output_folder = True, folder_name = None):
	if folder_name is None:
		output_path = os.path.join(get_parent_path(get_parent_path(file_path)), OUTPUT_FOLDER)
	else:
		output_path = os.path.join(get_parent_path(get_parent_path(file_path)), folder_name)

	# Add '/' at the end
	output_path = os.path.join(output_path, '')

	if output_folder == False:
		output_path = get_parent_path(file_path)

	# Magic: Prefix with file name! :facepalm:
	#output_path = os.path.join(output_path, get_file_name(file_path))

	for argument in arguments:
		output_path += argument

	return output_path

def join_file_path(file_path, file_name):
	return os.path.join(file_path, file_name)

def run_python_script(script_name, arguments):
	# python <script_name>
	command = PYTHON_COMMAND + ' ' + script_name

	# python <script_name> <arguments>
	# Add all the arguments to the command
	for argument in arguments:
		command += ' ' + argument

	print (command)
	os.system(command)

# You can't send in arguments here!
def run_paraview_script(script_name):
	# paraview --script= <script_name>
	command = PARAVIEW_COMMAND + ' --script=' + script_name
	os.system(command)

def run_shell_script(script_name, arguments):
	command = './' + script_name

	for argument in arguments:
		command += ' ' + argument

	os.system(command)

def get_dictionary(file_path, arguments):
	# extension is always BIN
	arguments.append(BIN_EXTENSION)
	dictionary_path = get_output_path(file_path, arguments, folder_name = DICTIONARY_FOLDER)

	with open(dictionary_path, 'rb') as handle:
		current_dictionary = pickle.loads(handle.read())

	return current_dictionary

def get_matrix(file_path, arguments):
	# extension is always BIN
	arguments.append(BIN_EXTENSION)

	# combine the filenames together
	arguments[0:2] = [join_strings(arguments[0:2])]

	dictionary_path = get_output_path(file_path, arguments, folder_name = MATRICES_FOLDER)

	with open(dictionary_path, 'rb') as handle:
		current_dictionary = pickle.loads(handle.read())

	return current_dictionary


def get_folder(file_path, folder_name = None):
	# This will give you a path with the folder name appended
	output_path = get_output_path(file_path, [], folder_name = folder_name)
	return output_path

# Create an output folder if it does not exist
def create_output_folder(file_path):
	output_path = get_output_folder(file_path)
	# Check if it exists
	if not os.path.exists(output_path):
		os.makedirs(output_path)

# Delete the output folder
def delete_output_folder(file_path):
	shutil.rmtree(get_output_folder(file_path))

# remove ds_store
def remove_ds_store(file_list):
	try:
		file_list.remove('.DS_Store')
	except:
		pass

def get_comparison_indices(s):
	# seperate on '-'
	current_timestep, sep, comparison_timestep = s.partition('-')
	# seperate on 'comparison_timestep.extension'
	comparison_timestep, sep, extension = comparison_timestep.partition('.')

	# idea is that both if you have both values in the filename,
	# and you join just the indice values together,
	# and sort them, it is equivalent to sorting the files

	# seperate current_timestep based on it's prefix
	prefix, sep, current_timestep = current_timestep.partition('_')
	prefix, sep, comparison_timestep = comparison_timestep.partition('_')

	return [current_timestep, comparison_timestep]

# used in run: sort files according to characters and numeric
def sort_comparison_files(s):
	[current_timestep, comparison_timestep] = get_comparison_indices(s)
	index = current_timestep + comparison_timestep

	# isnumeric only works for unicode strings
	if unicode(index, 'utf-8').isnumeric():
		return int(index)
	return float('inf')

# used in run: sort files according to characters and numeric
def sort_files(s):
	# seperate on '_'
	timestep, sep, index = s.partition('_')
	# seperate on 'index.extension'
	index, sep, extension = index.partition('.')
	# isnumeric only works for unicode strings
	if str(index).isnumeric():
		return int(index)
	return float('inf')

# used in run: sort files according to characters and numeric
def sort_files_numeric(s):
	# seperate on 'index.extension'
	index, sep, extension = s.partition('.')
	# isnumeric only works for unicode strings
	if str(index).isnumeric():
		return int(index)
	return float('inf')

# Pretty print the time taken
def pretty_print_time(seconds):
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	time_taken =  "%d:%02d:%02d" % (h, m, s)
	return time_taken

# save dictionary to file
def save_dictionary(dictionary, file_path, file_name, identifier):
	dictionary_file_arguments = [file_name, identifier, BIN_EXTENSION]
	dictionary_file_path = get_output_path(file_path, dictionary_file_arguments, folder_name = DICTIONARY_FOLDER)
	with open(dictionary_file_path, 'wb') as handle:
		pickle.dump(dictionary, handle)

# save dictionary to file
def save_object(dictionary, file_path, file_name):
	dictionary_file_arguments = [file_name, BIN_EXTENSION]
	dictionary_file_path = get_output_path(file_path, dictionary_file_arguments, folder_name = DICTIONARY_FOLDER)
	with open(dictionary_file_path, 'wb') as handle:
		pickle.dump(dictionary, handle)

# used in edit-distance to save intermediate matrices
def save_matrix(dictionary, file_path, filenames, identifier):
	matrix_file_arguments = [join_strings(filenames), identifier, BIN_EXTENSION]
	matrix_file_path = get_output_path(file_path, matrix_file_arguments, folder_name = MATRICES_FOLDER)
	with open(matrix_file_path, 'wb') as handle:
		pickle.dump(dictionary, handle)

def find_coords(vertex_index):
	x_dim = y_dim = z_dim = 300
	z = vertex_index/(x_dim*y_dim)
	xy = vertex_index - z * x_dim * y_dim
	y = xy/x_dim
	x = xy - y * x_dim
	return (x, y, z)

def find_distance(vertex_index):
	[x, y, z] = find_coords(vertex_index)
	return math.sqrt(x*x + y*y + z*z)

# used for printing in split-make-graph-*.py
def get_label(index, pairs, mappings, labels):
	scalar = round(labels[index], 3)
	label = str(index) + " [" + str(pairs[index]) + "]" + "\\n"
	label += str(mappings[index]) + "\\n"
	label += str(scalar)
	label_attribute = "label=""\"" + label +"\""
	return label_attribute

def get_node(index, pairs, mappings, labels):
	label_attribute = get_label(index, pairs, mappings, labels)
	return str(mappings[index]) + " [" + label_attribute + "]\n"

def get_connectivity(index, parent, mapping):
	if parent[index] != 0:
		node1 = str(mapping[index])
		connector = ' -> '
		node2 = str(mapping[parent[index]])
		line = node2 + connector + node1 + '\n'
		return line
	else:
		return ''

def save_string(string, file_path):
	with open(file_path, "w") as string_file:
		string_file.write(string)

class cfile:
	def __init__(self, name, mode = 'w'):
		self.f = open(name, mode)
	def __enter__(self):
		return self.f
	def __exit__(self, exc_type, exc_value, traceback):
		self.f.close()
	def wl(self, string):
		self.f.writelines(string + '\n')
		return None
	def close(self):
		self.f.close()

# i remember writing countless System.out.println()
# to output to a HTML file almost three years back
def save_stability(scalars, birth_pairs, stability_file_path):
	stability_file = cfile(stability_file_path, 'w')
	coordinates = {}
	cell_indices = {}
	cell_index = 0
	num_edges = len(birth_pairs.keys())
	num_nodes = 2 * len(birth_pairs.keys())

	stability_file.wl('# vtk DataFile Version 4.2')
	stability_file.wl('persistence pairs')
	stability_file.wl('ASCII')
	stability_file.wl('DATASET UNSTRUCTURED_GRID')
	stability_file.wl('POINTS ' + str(num_nodes) + ' float')

	for vertex in birth_pairs.keys():
		for vertex_index in [vertex, birth_pairs[vertex]]:
			z = 0
			y = vertex_index/400
			x = vertex_index - (400 * y)
			coordinates[vertex_index] = [x, y, z]
			cell_indices[vertex_index] = cell_index
			cell_index += 1
			stability_file.wl(" ".join(map(str, coordinates[vertex_index])))

	stability_file.wl('')

	stability_file.wl('CELLS ' + str(num_edges) + ' ' + str(3 * num_edges))
	for birth_vertex in birth_pairs.keys():
		death_vertex = birth_pairs[birth_vertex]
		cell_data = [2, cell_indices[birth_vertex], cell_indices[death_vertex]]
		stability_file.wl(" ".join(map(str, cell_data)))

	stability_file.wl('')

	stability_file.wl('CELL_TYPES ' + str(num_edges))
	for edge_type in range(0, num_edges):
		stability_file.wl('3')

	stability_file.wl('')

	stability_file.wl('CELL_DATA ' + str(num_edges))
	stability_file.wl('SCALARS Persistence float')
	stability_file.wl('LOOKUP_TABLE default')
	for birth_vertex in birth_pairs.keys():
		death_vertex = birth_pairs[birth_vertex]
		stability_file.wl(str(scalars[death_vertex] - scalars[birth_vertex]))
	stability_file.wl('')

	stability_file.wl('FIELD FieldData 2')
	stability_file.wl('birth_vertex 1 ' + str(num_edges) + ' int')
	for birth_vertex in birth_pairs.keys():
		stability_file.wl(str(birth_vertex))
	stability_file.wl('')

	stability_file.wl('death_vertex 1 '+ str(num_edges) + ' int')
	for birth_vertex in birth_pairs.keys():
		stability_file.wl(str(birth_pairs[birth_vertex]))
	stability_file.wl('')

	stability_file.wl('POINT_DATA ' +  str(num_nodes))
	stability_file.wl('FIELD FieldData 2')
	stability_file.wl('Scalar 1 ' + str(num_nodes) + ' float')
	for vertex in birth_pairs.keys():
		for vertex_index in [vertex, birth_pairs[vertex]]:
				stability_file.wl(str(scalars[vertex_index]))

	stability_file.wl('')

	stability_file.wl('VertexIndex 1 ' + str(num_nodes) + ' int')
	for vertex in birth_pairs.keys():
		for vertex_index in [vertex, birth_pairs[vertex]]:
				stability_file.wl(str(vertex_index))
	stability_file.wl('')
	stability_file.close()


# i remember writing countless System.out.println()
# to output to a HTML file almost three years back
def save_segmentation(point_to_coordinates, triangle_index_to_points, processed_triangles, segmentation_file_path):
	stability_file = cfile(segmentation_file_path, 'w')

	num_points = len(point_to_coordinates.keys())
	num_triangles = len(processed_triangles.keys())

	stability_file.wl('# vtk DataFile Version 4.2')
	stability_file.wl('split tree segmentation')
	stability_file.wl('ASCII')
	stability_file.wl('DATASET UNSTRUCTURED_GRID')
	stability_file.wl('POINTS ' + str(num_points) + ' float')

	for point in sorted(point_to_coordinates.keys()):
		stability_file.wl(" ".join(map(str, point_to_coordinates[point])))
	stability_file.wl('')

	stability_file.wl('CELLS ' + str(num_triangles) + ' ' + str(4 * num_triangles))
	for triangle in sorted(processed_triangles.keys()):
		triangle_points = map(int, list(triangle_index_to_points[triangle]))
		cell_data = [3]
		cell_data.extend(triangle_points)
		stability_file.wl(" ".join(map(str, cell_data)))

	stability_file.wl('')

	stability_file.wl('CELL_TYPES ' + str(num_triangles))
	for edge_type in range(0, num_triangles):
		stability_file.wl('5')

	stability_file.wl('')

	stability_file.wl('CELL_DATA ' + str(num_triangles))
	stability_file.wl('SCALARS segmentation int')
	stability_file.wl('LOOKUP_TABLE default')

	for triangle in sorted(processed_triangles.keys()):
		segmentation = processed_triangles[triangle]
		stability_file.wl(str(segmentation))
	stability_file.wl('')

	stability_file.close()


# i remember writing countless System.out.println()
# to output to a HTML file almost three years back
# the only difference from the above method and this is that
# cell data type is int above and it is float now
# what a facepalm
def save_confidence_region(point_to_coordinates, triangle_index_to_points, processed_triangles, segmentation_file_path):
	stability_file = cfile(segmentation_file_path, 'w')

	num_points = len(point_to_coordinates.keys())
	num_triangles = len(processed_triangles.keys())

	stability_file.wl('# vtk DataFile Version 4.2')
	stability_file.wl('split tree segmentation')
	stability_file.wl('ASCII')
	stability_file.wl('DATASET UNSTRUCTURED_GRID')
	stability_file.wl('POINTS ' + str(num_points) + ' float')

	for point in sorted(point_to_coordinates.keys()):
		stability_file.wl(" ".join(map(str, point_to_coordinates[point])))
	stability_file.wl('')

	stability_file.wl('CELLS ' + str(num_triangles) + ' ' + str(4 * num_triangles))
	for triangle in sorted(processed_triangles.keys()):
		triangle_points = map(int, list(triangle_index_to_points[triangle]))
		cell_data = [3]
		cell_data.extend(triangle_points)
		stability_file.wl(" ".join(map(str, cell_data)))

	stability_file.wl('')

	stability_file.wl('CELL_TYPES ' + str(num_triangles))
	for edge_type in range(0, num_triangles):
		stability_file.wl('5')

	stability_file.wl('')

	stability_file.wl('CELL_DATA ' + str(num_triangles))
	stability_file.wl('SCALARS segmentation float')
	stability_file.wl('LOOKUP_TABLE default')

	for triangle in sorted(processed_triangles.keys()):
		confidence = processed_triangles[triangle]
		stability_file.wl(str(confidence))
	stability_file.wl('')

	stability_file.close()

# indices are as [<lower, =lower, <lower, upper>, = upper, >upper]
def calculate_isoband_index(scalar, lower, upper):
	if (scalar < lower):
		return '0'
	# check for equality with lower
	elif(abs(scalar-lower) < 1e-8):
		return '1'
	elif (scalar > upper):
		return '4'
	# check for equality with upper
	elif(abs(upper-scalar) < 1e-8):
		return '3'
	else:
		return '2'

# indices are as [<lower, =lower, <lower, upper>, = upper, >upper]
def calculate_cheap_isoband_index(scalar, lower, upper):
	if (scalar < lower):
		return '0'
	elif (scalar > upper):
		return '2'
	else:
		return '1'

