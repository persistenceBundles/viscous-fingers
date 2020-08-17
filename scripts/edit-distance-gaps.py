import os, time
os.system('date')
start_time = time.time()

import itertools
from helper import *

# Only me and god knows how this code works
# In a few days only god will know

# Take the names of files as arguments
filename1 =  sys.argv[1]
filename2 = sys.argv[2]

# Get the paths
file_path = os.path.abspath(inspect.getfile(inspect.currentframe()))
dictionary1 = get_output_path(file_path, [filename1], folder_name = DICTIONARY_FOLDER)
dictionary2 = get_output_path(file_path, [filename2], folder_name = DICTIONARY_FOLDER)
results_path = get_output_path(file_path, [EDIT_DISTANCE_RESULT, CSV_EXTENSION], folder_name = RESULTS_FOLDER)

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

# Get size of both the trees
size1 = len(parent1.keys())
size2 = len(parent2.keys())

# Dictionaries to initialize the distance
Q = {}
Q1 = {}
Q2 = {}

# Dictionaries to initialize the sequence
S = {}
S1 = {}
S2 = {}

# How much are we planning to compare?
extents = [right1[1], right2[1]]
print (extents)

# Gap cost associated with the first tree
def gap1(i):
	return difference1[i]/2.0

# Gap cost associated with the second tree
def gap2(j):
	return difference2[j]/2.0

# Relabel cost
def relabel(i, j):
	return abs(label1[i] - label2[j])

# Check for bad extents in the first tree
def checkBadExtents1(i1, i):
	return (i < 1 or i > size1 or i1 < 1 or i1 > size1 or i < i1)

# Check for bad extents in the second tree
def checkBadExtents2(j1, j):
	return (j < 1 or j > size2 or j1 < 1 or j1 > size2 or j < j1)

# Check for bad extents in both the trees
def checkGoodExtents(i1, i, j1, j):
	t1 = not checkBadExtents1(i1, i)
	t2 = not checkBadExtents2(j1, j)
	return (t1, t2)

# Create an entry in corresponding dictionary
def createEntry(dictionary, i1, i, j1, j):
	# Dirty check to see if entry already present
	try:
		dictionary[i1][i][j1][j] += 0
	except:
		if i1 not in dictionary.keys():
			dictionary[i1] = {}
		if i not in dictionary[i1].keys():
			dictionary[i1][i] = {}
		if j1 not in dictionary[i1][i].keys():
			dictionary[i1][i][j1] = {}
		if j not in dictionary[i1][i][j1].keys():
			dictionary[i1][i][j1][j] = {}

def isComputedEarlier(dictionary, i1, i, j1, j):
	# simple check to see if value of key is not zero
	return bool(dictionary[i1][i][j1][j])

# Return index of minimum if there are multiple minima present
def locate_minimum(elements):
	smallest = min(elements)
	return min([index for index, element in enumerate(elements) if smallest == element])

def Qf(i1, i, j1, j):
	# create an empty key if already not present
	createEntry(Q, i1, i, j1, j)

	# check if extents are fine
	[extent1, extent2] = checkGoodExtents(i1, i, j1, j)
	
	# Boundary Condition: if both are outside bounds, assign zero
	if ((not extent1) and (not extent2)):
		Q[i1][i][j1][j] = 0
		return Q[i1][i][j1][j]

	# if either one is outside bounds, assign infinity
	elif ((not extent1) or (not extent2)):
		Q[i1][i][j1][j] = INFINITY
		return Q[i1][i][j1][j]

	# if already created earlier, return that value instead of computing again
	if (isComputedEarlier(Q, i1, i, j1, j)):
		return Q[i1][i][j1][j]
	else:
		# None of i or j is a gap point, then i must be matched to j
		zeroth = Qf(i1, i - 1, j1, j - 1) + relabel(i, j)

		# i is a gap point
		first = Q1f(i1, i, j1, j)

		# j is a gap point
		second = Q2f(i1, i, j1, j)

		# store minimum value of all three cases
		Q[i1][i][j1][j] = min(zeroth, first, second)

		#print ('Q', i1, i, j1, j, Q[i1][i][j1][j])

		# Adhitya: thought about optimising and decided not to. 
		# Adhitya from future: Do *not* try to optimise!
		cases = [zeroth, first, second]
		# Store the case index [whatever that means!]
		case_index = locate_minimum(cases)
		# Create a entry for this new sequence element
		createEntry(S, i1, i, j1, j)
		# Store previous element
		if(case_index == 0):
			message = str(i) + ' relabeled to ' + str(j)
			cost = relabel(i, j)
			S[i1][i][j1][j] = [S_MATRIX_IDENTIFIER, i1, i - 1, j1, j - 1, RELABEL_IDENTIFIER, message, cost]
		elif (case_index == 1):
			message = str(i) + ' is a gap point in T1'
			cost = UNKNOWN_COST
			S[i1][i][j1][j] = [S1_MATRIX_IDENTIFIER, i1, i, j1, j, T1_GENERIC_GAP_IDENTIFIER, message, cost]
		else:
			message = str(j) + ' is a gap point in T2'
			cost = UNKNOWN_COST
			S[i1][i][j1][j] = [S2_MATRIX_IDENTIFIER, i1, i, j1, j, T2_GENERIC_GAP_IDENTIFIER, message, cost]

		return Q[i1][i][j1][j]

# i is a gap point
def Q1f(i1, i, j1, j):
	# create an empty key if already not present
	createEntry(Q1, i1, i, j1, j)

	# check if extents are fine
	[extent1, extent2] = checkGoodExtents(i1, i, j1, j)

	# there is a unique matching between T1[1][i] with an empty tree: we have i gap points
	if (not extent2):
		Q1[i1][i][j1][j] = gap1(i)
		return Q1[i1][i][j1][j]

	# it is impossible to match an empty tree with T2[1][j] such that the former ends with a gap node
	elif (not extent1):
		Q1[i1][i][j1][j] = INFINITY
		return Q1[i1][i][j1][j]

	# if already created earlier, return that value instead of computing again
	if (isComputedEarlier(Q1, i1, i, j1, j)):
		return Q1[i1][i][j1][j]
	else:
		# initialize a minimum
		minimum = INFINITY
		
		# initialize a minimum index	
		minimum_k = None

		# if parent(i) is a gap node and i is its right child
		# then i is continuing a preexisting gap, hence gets penalized
		for k in range(j1, j + 1):

			# Note: left child of parent(i) has index parent(i) + 1
			# For some k in [j1, j]; with parent(i) being a gap node, 
			# the subforest T1[i1][parent(i)] is matched to a subforest T2[j1][k]
			# therefore, subtree rooted at [parent(i) + 1] will be matched with
			# the remaining part of T2[j1][j] = T2[k+1][j]

			zeroth = Qf(parent1[i] + 1, i - 1, k + 1, j)
			first = Q1f(i1, parent1[i], j1, k)

			# store the values for each case
			cases = [minimum, first + zeroth + gap1(i)]
			# store the case index first [whatever that means :P]
			case_index = locate_minimum(cases)
			# only then does the index of k change
			if case_index == 1: minimum_k = k
			minimum = min(cases)

		# if [parent(i) does not exist] or [parent(i) exists and is not a gap node]
		# then i is starting a new gap, hence gets penalized
		zeroth = Qf(i1, i - 1, j1, j) + gap1(i)

		# if parent(i) is a gap node and i is its left child
		# then i is continuing a preexisting gap, hence gets penalized
		first = Q1f(i1, i - 1, j1, j) + gap1(i)

		Q1[i1][i][j1][j] = min(first, zeroth, minimum)

		# store the sequence [This preferential order is important]
		cases = [first, zeroth, minimum]
		case_index = locate_minimum(cases)
		# Create a entry for this new sequence element
		createEntry(S1, i1, i, j1, j)
		if (case_index == 0):
			message = str(i) + ' is continuing a preexisting gap in T1'
			cost = gap1(i)
			S1[i1][i][j1][j] = [S1_MATRIX_IDENTIFIER, i1, i - 1, j1, j, T1_CONTINUING_GAP_IDENTIFIER, message, cost]
		elif(case_index == 1):
			message = str(i) + ' is starting a new gap in T1'
			cost = gap1(i)
			S1[i1][i][j1][j] = [S_MATRIX_IDENTIFIER, i1, i - 1, j1, j, T1_STARTING_GAP_IDENTIFIER, message, cost]
		else:
			k = minimum_k
			S1[i1][i][j1][j] = [[S_MATRIX_IDENTIFIER, parent1[i] + 1, i - 1, k + 1, j], [S1_MATRIX_IDENTIFIER, i1, parent1[i], j1, k]]

		#print ('Q1', i1, i, j1, j, Q1[i1][i][j1][j])
		return Q1[i1][i][j1][j]

# j is a gap point
def Q2f(i1, i, j1, j):
	# create an empty key if already not present
	createEntry(Q2, i1, i, j1, j)

	# check if extents are fine
	[extent1, extent2] = checkGoodExtents(i1, i, j1, j)

	# it is impossible to match an empty tree with T1[1][j] such that the former ends with a gap node
	if (not extent2):
		Q2[i1][i][j1][j] = INFINITY
		return Q2[i1][i][j1][j]

	# there is a unique matching between T2[1][i] with an empty tree: we have j gap points
	elif (not extent1):
		Q2[i1][i][j1][j] = gap2(j)
		return Q2[i1][i][j1][j]

	# if already created earlier, return that value instead of computing again
	if (isComputedEarlier(Q2, i1, i, j1, j)):
		return Q2[i1][i][j1][j]

	# if parent(j) is a gap node and j is its right child
	# then j is continuing a preexisting gap, hence gets penalized
	else:
		# initialize a minimum
		minimum = INFINITY

		# initialize a minimum index
		minimum_k = None

		for k in range(i1, i + 1):

			# Note: left child of parent(j) has index parent(j) + 1
			# For some k in [i1, i]; with parent(j) being a gap node, 
			# the subforest T1[i1][k] is matched to a subforest T2[j1][parent(j)]
			# therefore, subtree rooted at [parent(j) + 1] will be matched with
			# the remaining part of T1[i1][i] = T1[k+1][i]

			zeroth = Qf(k + 1, i, parent2[j] + 1, j - 1)
			second = Q2f(i1, k, j1, parent2[j])

			# store the values for each case
			cases = [minimum, second + zeroth + gap2(j)]
			# store the case index first [whatever that means :P]
			case_index = locate_minimum(cases)
			# only then does the index of k change
			if case_index == 1: minimum_k = k

			minimum = min(cases)

		# if [parent(j) does not exist] or [parent(j) exists and is not a gap node]
		# then j is starting a new gap, hence gets penalized
		zeroth = Qf(i1, i, j1, j - 1) + gap2(j)

		# if parent(j) is a gap node and j is its left child
		# then j is continuing a preexisting gap, hence gets penalized
		second = Q2f(i1, i, j1, j - 1) + gap2(j)

		Q2[i1][i][j1][j] = min(second, zeroth, minimum)

		# store the sequence [This preferential order is important]
		cases = [second, zeroth, minimum]
		case_index = locate_minimum(cases)
		# Create a entry for this new sequence element
		createEntry(S2, i1, i, j1, j)
		if (case_index == 0):
			message = str(j) + ' is continuing a preexisting gap in T2'
			cost = gap2(j)
			S2[i1][i][j1][j] = [S2_MATRIX_IDENTIFIER, i1, i, j1, j - 1, T2_CONTINUING_GAP_IDENTIFIER, message, cost]
		elif(case_index == 1):
			message = str(j) + ' is starting a new gap in T2'
			cost = gap2(j)
			S2[i1][i][j1][j] = [S_MATRIX_IDENTIFIER, i1, i, j1, j - 1, T2_STARTING_GAP_IDENTIFIER, message, cost]
		else:
			k = minimum_k
			S2[i1][i][j1][j] = [[S_MATRIX_IDENTIFIER, k + 1, i, parent2[j] + 1, j - 1], [S2_MATRIX_IDENTIFIER, i1, k, j1, parent2[j]]]

		#print ('Q2', i1, i, j1, j, Q2[i1][i][j1][j])
		return Q2[i1][i][j1][j]

# Iterate across both the trees!
for i1, j1 in itertools.product(range(1, size1 + 1), range(1, size2 + 1)):
	for i in range(i1, right1[i1] + 1):
		for j in range(j1, right2[j1] + 1):
			#print (i1, i, j1, j)
			Q1[i1][i][j1][j] = Q1f(i1, i, j1, j)
			Q2[i1][i][j1][j] = Q2f(i1, i, j1, j)
			Q[i1][i][j1][j] = Qf(i1, i, j1, j)
			#print ('')

# Complete distance is given at the end the matrix
difference = Q[1][right1[1]][1][right2[1]]
difference = "%.9f" % difference

# Stop timer
seconds = time.time() - start_time
time_taken = pretty_print_time(seconds)

# Write distance to file
csvfile = open(results_path, 'a')

# Change this when the filename format changes
timestep_value = filename1.split('tv_')[1]
csvfile.write(timestep_value + ',' + str(difference) +'\n')
csvfile.close()

# Save intermediate matrices
save_matrix(Q, [filename1, filename2], Q_IDENTIFIER)
save_matrix(Q1, [filename1, filename2], Q1_IDENTIFIER)
save_matrix(Q2, [filename1, filename2], Q2_IDENTIFIER)

save_matrix(S, [filename1, filename2], S_IDENTIFIER)
save_matrix(S1, [filename1, filename2], S1_IDENTIFIER)
save_matrix(S2, [filename1, filename2], S2_IDENTIFIER)

print ('Difference: ', difference, 'Time Taken: ', time_taken)
os.system('date')
