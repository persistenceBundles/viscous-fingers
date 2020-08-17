from paraview.simple import *
import csv, os
from datetime import datetime
from helper import *
from collections import defaultdict
from make_hierarchy import *
from collections import Counter

paraview.simple._DisableFirstRenderCameraReset()

# start timer
startTime = datetime.now()

# initialize Path variables
# create a new 'Legacy VTK Reader'
full_file_name = 'tv_2.vtk'
num_subdivisions = 2
simplification_percentage = 1.2

parent_path = cwd()
data_path = get_input_path(parent_path)
file_path = join_file_path(data_path, full_file_name)
file_name = get_file_name(full_file_name)

# create a new 'Legacy VTK Reader'
inputFile = LegacyVTKReader(FileNames=[file_path])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
inputFileDisplay = Show(inputFile, renderView1)
inputFileDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera()
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [199.5, 24.5, 10000.0]
renderView1.CameraFocalPoint = [199.5, 24.5, 0.0]
inputFileDisplay.SetScalarBarVisibility(renderView1, False)
renderView1.Update()

# get color transfer function/color map for 'magnitude'
magnitudeLUT = GetColorTransferFunction('magnitude')

# create a new 'Extract Surface'
extractSurface1 = ExtractSurface(Input=inputFile)
extractSurface1Display = Show(extractSurface1, renderView1)
extractSurface1Display.Representation = 'Surface'
Hide(inputFile, renderView1)
extractSurface1Display.SetScalarBarVisibility(renderView1, False)
renderView1.Update()

# create a new 'Triangulate'
triangulate1 = Triangulate(Input=extractSurface1)
triangulate1Display = Show(triangulate1, renderView1)
triangulate1Display.Representation = 'Surface'
Hide(extractSurface1, renderView1)
triangulate1Display.SetScalarBarVisibility(renderView1, False)
renderView1.Update()

# create a new 'Loop Subdivision'
loopSubdivision1 = LoopSubdivision(Input=triangulate1)
loopSubdivision1.NumberofSubdivisions = num_subdivisions
loopSubdivision1Display = Show(loopSubdivision1, renderView1)
loopSubdivision1Display.Representation = 'Surface'
Hide(triangulate1, renderView1)
loopSubdivision1Display.SetScalarBarVisibility(renderView1, False)
renderView1.Update()

# create a new 'Clean to Grid'
vtkFile = CleantoGrid(Input=loopSubdivision1)
vtkFileDisplay = Show(vtkFile, renderView1)
vtkFileDisplay.Representation = 'Surface'
Hide(loopSubdivision1, renderView1)
vtkFileDisplay.SetScalarBarVisibility(renderView1, False)
renderView1.Update()

# create a new 'Extract Subset'
#vtkFile = ExtractSubset(Input=inputFile)
#vtkFile.VOI = [0, 399, 0, 49, 0, 0]
#vtkFileDisplay = Show(vtkFile, renderView1)
#vtkFileDisplay.Representation = 'Slice'
#vtkFileDisplay.SetScalarBarVisibility(renderView1, True)
#renderView1.Update()
#Hide(inputFile, renderView1)

# get layout
layout1 = GetLayout()
layout1.SplitHorizontal(0, 0.5)
SetActiveView(None)

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [1038, 1179]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.StereoType = 0
renderView2.Background = [0.32, 0.34, 0.43]

# place view in the layout
layout1.AssignView(2, renderView2)

# create a new 'TTK PersistenceDiagram'
persistenceDiagram = TTKPersistenceDiagram(Input=vtkFile)
#persistenceDiagram.UseInputOffsetField = 1
#persistenceDiagram.InputOffsetField = ''
persistenceDiagramDisplay = Show(persistenceDiagram, renderView2)
persistenceDiagramDisplay.Representation = 'Surface'

# reset view to fit data
renderView2.ResetCamera()
renderView2.InteractionMode = '2D'
renderView2.CameraPosition = [0.6912000179290771, 0.710545003414154, 10000.0]
renderView2.CameraFocalPoint = [0.6912000179290771, 0.710545003414154, 0.0]
renderView2.Update()

# create a new 'Threshold'
persistencePairsThreshold = Threshold(Input=persistenceDiagram)
persistencePairsThreshold.Scalars = ['CELLS', 'PairType']
persistencePairsThreshold.ThresholdRange = [-1, 1]

# show data in view
persistencePairsThresholdDisplay = Show(persistencePairsThreshold, renderView2)
persistencePairsThresholdDisplay.Representation = 'Surface'
Hide(persistenceDiagram, renderView2)
renderView2.Update()

# find max persistence by iterating across the diagram!
persistence_data = servermanager.Fetch(persistencePairsThreshold)

# Get the number of persistent points and arcs
num_persistent_points = persistencePairsThreshold.GetDataInformation().GetNumberOfPoints()
num_persistent_cells = persistencePairsThreshold.GetDataInformation().GetNumberOfCells()

max_persistence = 0
for index in range(num_persistent_cells):
	current_persistence = persistence_data.GetCellData().GetArray('Persistence').GetValue(index)
	max_persistence = max(current_persistence, max_persistence)

# filter all persistent pairs above minimum persistence threshold
min_persistence = (simplification_percentage *  max_persistence) / 100.0

# create a new 'Threshold'
persistenceThreshold = Threshold(Input=persistencePairsThreshold)
persistenceThreshold.Scalars = ['CELLS', 'Persistence']
persistenceThreshold.ThresholdRange = [min_persistence, max_persistence]
persistenceThresholdDisplay = Show(persistenceThreshold, renderView2)
persistenceThresholdDisplay.Representation = 'Surface'
Hide(persistencePairsThreshold, renderView2)
renderView2.Update()

# set active view
SetActiveView(renderView1)

# create a new 'TTK TopologicalSimplification'
topologicalSimplification = TTKTopologicalSimplification(Domain=vtkFile,
    Constraints=persistenceThreshold)
#topologicalSimplification.UseInputOffsetField = 1
#topologicalSimplification.InputOffsetField = ''
topologicalSimplificationDisplay = Show(topologicalSimplification, renderView1)
topologicalSimplificationDisplay.Representation = 'Surface'
Hide(vtkFile, renderView1)
topologicalSimplificationDisplay.SetScalarBarVisibility(renderView1, False)
renderView1.Update()

# create a new 'TTK Merge and Contour Tree (FTM)'
contourTree = TTKMergeandContourTreeFTM(Input=topologicalSimplification)
#contourTree.UseInputOffsetScalarField = 1
contourTree.TreeType = 'Split Tree'
contourTreeDisplay = Show(contourTree, renderView1)
contourTreeDisplay.Representation = 'Surface'
Hide(topologicalSimplification, renderView1)
contourTreeDisplay.SetScalarBarVisibility(renderView1, False)
contourTreeDisplay_1 = Show(OutputPort(contourTree, 1), renderView1)
contourTreeDisplay_1.Representation = 'Surface'
Hide(topologicalSimplification, renderView1)
contourTreeDisplay_1.SetScalarBarVisibility(renderView1, False)
contourTreeDisplay_2 = Show(OutputPort(contourTree, 2), renderView1)
contourTreeDisplay_2.Representation = 'Surface'
Hide(topologicalSimplification, renderView1)
contourTreeDisplay_2.SetScalarBarVisibility(renderView1, False)
renderView1.Update()

# create a new 'Threshold'
segmentationThreshold = Threshold(Input=OutputPort(contourTree, 2))
segmentationThreshold.Scalars = ['POINTS', 'RegionType']
segmentationThreshold.ThresholdRange = [3, 4]
segmentationThresholdDisplay = Show(segmentationThreshold, renderView1)
segmentationThresholdDisplay.Representation = 'Surface'
#Hide(vtkFile, renderView1)
renderView1.Update()

# create a new 'Threshold'
segmentationThreshold2 = Threshold(Input=OutputPort(contourTree, 2))
segmentationThreshold2.Scalars = ['POINTS', 'RegionType']
segmentationThreshold2.ThresholdRange = [1, 1]
segmentationThresholdDisplay2 = Show(segmentationThreshold2, renderView1)
segmentationThresholdDisplay2.Representation = 'Surface'
#Hide(vtkFile, renderView1)
Hide(segmentationThreshold2, renderView1)
renderView1.Update()

# get type for usage below
tree_type = get_tree_type(contourTree.TreeType)

# get color transfer function/color map for 'NodeType'
nodeTypeLUT = GetColorTransferFunction('NodeType')

# hide data in view
Hide(OutputPort(contourTree, 2), renderView1)
SetActiveSource(vtkFile)
#vtkFileDisplay = Show(vtkFile, renderView1)
#vtkFileDisplay.SetScalarBarVisibility(renderView1, False)

# set active source
SetActiveSource(contourTree)
contourTreeDisplay_1.ScaleTransferFunction.RescaleTransferFunction(0.0, 1.17578133675e-38)
contourTreeDisplay_1.OpacityTransferFunction.RescaleTransferFunction(0.0, 1.17578133675e-38)
segmentationIdLUT = GetColorTransferFunction('RegionSize')

# create a new 'Extract Surface'
extractSurface = ExtractSurface(Input=OutputPort(contourTree, 1))
extractSurfaceDisplay = Show(extractSurface, renderView1)
extractSurfaceDisplay.Representation = 'Surface'
Hide(OutputPort(contourTree, 1), renderView1)
extractSurfaceDisplay.SetScalarBarVisibility(renderView1, False)
renderView1.Update()

# create a new 'Tube'
tube = Tube(Input=extractSurface)
tube.Vectors = [None, '']
tube.Radius = 1.7
tubeDisplay = Show(tube, renderView1)
tubeDisplay.Representation = 'Surface'
Hide(extractSurface, renderView1)
tubeDisplay.SetScalarBarVisibility(renderView1, False)
HideScalarBarIfNotNeeded(segmentationIdLUT, renderView1)
tubeDisplay.DiffuseColor = [1.0, 1.0, 0.0]
renderView1.Update()

SetActiveSource(contourTree)

# create a new 'TTK SphereFromPoint'
tTKSphereFromPoint = TTKSphereFromPoint(Input=contourTree)
tTKSphereFromPoint.Radius = 2.0
tTKSphereFromPointDisplay = Show(tTKSphereFromPoint, renderView1)
tTKSphereFromPointDisplay.Representation = 'Surface'
Hide(contourTree, renderView1)
tTKSphereFromPointDisplay.SetScalarBarVisibility(renderView1, False)
renderView1.Update()

# create a new 'Threshold'
tTKSphereFromPointThreshold = Threshold(Input=tTKSphereFromPoint)
tTKSphereFromPointThreshold.Scalars = ['POINTS', 'NodeType']
tTKSphereFromPointThreshold.ThresholdRange = [0, 4]
tTKSphereFromPointThresholdDisplay = Show(tTKSphereFromPointThreshold, renderView1)
tTKSphereFromPointThresholdDisplay.Representation = 'Surface'
Hide(tTKSphereFromPoint, renderView1)
renderView1.Update()

# split cell
layout1.SplitVertical(1, 0.5)
SetActiveView(None)

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
#spreadSheetView1.BlockSize = 1024L

# place view in the layout
layout1.AssignView(4, spreadSheetView1)

# show data in view
contourForestDisplayArcs = Show(OutputPort(contourTree, 1), spreadSheetView1)
spreadSheetView1.FieldAssociation = 'Cell Data'
arcs_file_arguments = [tree_type, ARCS_INFIX, file_name, CSV_EXTENSION]
arcs_file_path = get_output_path(file_path, arcs_file_arguments, folder_name = INTERMEDIATE_FOLDER)
ExportView(arcs_file_path, view=spreadSheetView1)

# Write node of contour tree to file
contourForestDisplayNodes = Show(contourTree, spreadSheetView1)
spreadSheetView1.FieldAssociation = 'Point Data'
nodes_file_arguments = [tree_type, NODES_INFIX, file_name, CSV_EXTENSION]
nodes_file_path = get_output_path(file_path, nodes_file_arguments, folder_name = INTERMEDIATE_FOLDER)
ExportView(nodes_file_path, view=spreadSheetView1)

SetActiveView(renderView2)
layout1.SplitVertical(2, 0.5)
SetActiveView(None)

# Create a new 'SpreadSheet View'
spreadSheetView2 = CreateView('SpreadSheetView')
spreadSheetView2.ColumnToSort = ''
#spreadSheetView2.BlockSize = 1024L
layout1.AssignView(6, spreadSheetView2)

# show data in view
contourTreeDisplay_4 = Show(contourTree, spreadSheetView2)
SetActiveView(renderView1)
SetActiveSource(contourTree)
SetActiveSource(topologicalSimplification)
SetActiveView(renderView2)
Hide(persistenceThreshold, renderView2)

# create a new 'TTK PersistenceDiagram'
simplififedPersistenceDiagram = TTKPersistenceDiagram(Input=topologicalSimplification)
#simplififedPersistenceDiagram.UseInputOffsetField = 1
simplififedPersistenceDiagramDisplay = Show(simplififedPersistenceDiagram, renderView2)
simplififedPersistenceDiagramDisplay.Representation = 'Surface'

SetActiveView(renderView1)
Hide(tTKSphereFromPointThreshold, renderView1)
Hide(tube, renderView1)

SetActiveSource(topologicalSimplification)
renderView1.Update()

# initialize dictionaries
scalars = {}
nodes = {}
birth_pairs = {}
arcs = defaultdict(list)

# Read the merge tree nodes
with open(nodes_file_path, 'r') as csvfile:
	csvfile.readline()
	spamreader = csv.reader(csvfile, delimiter=' ')
	for r in spamreader:
		row = r[0].split(',')
		node_id = int(row[0])
		scalar_value = float(row[1])
		vertex_index = int(row[2])

		# store above values to dictionaries
		nodes[node_id] = vertex_index
		scalars[vertex_index] = scalar_value

# Write the merge tree to file
tree_file_arguments = [tree_type, TREE_INFIX, file_name, CSV_EXTENSION]
tree_file_path = get_output_path(file_path, tree_file_arguments, folder_name = TREES_FOLDER)

tree_file = open(tree_file_path, 'w')
fieldnames = ['Node:0', 'Node:1', 'Scalar:0', 'Scalar:1']
writer = csv.writer(tree_file, delimiter=',')
writer.writerow(fieldnames)

# Read the intermediate arcs file
with open(arcs_file_path, 'r') as csvfile:
	csvfile.readline()
	spamreader = csv.reader(csvfile, delimiter=' ')
	for index, r in enumerate(spamreader):
		row = r[0].split(',')
		upNodeId = int(row[1])
		downNodeId = int(row[2])

		up_vertex_index = nodes[upNodeId]
		down_vertex_index = nodes[downNodeId]

		# store arcs in a dictionary
		arcs[up_vertex_index].append(down_vertex_index)
		arcs[down_vertex_index].append(up_vertex_index)

		content = [up_vertex_index, down_vertex_index, scalars[up_vertex_index], scalars[down_vertex_index]]
		writer.writerow(content)

tree_file.close()

# Write persistent pairs after thresholding to file
pairs_file_arguments = [tree_type, PAIRS_INFIX, file_name, CSV_EXTENSION]
pairs_file_path = get_output_path(file_path, pairs_file_arguments, folder_name = PAIRS_FOLDER)

pairs_file = open(pairs_file_path, 'w')
fieldnames = ['Birth', 'Death']
writer = csv.writer(pairs_file, delimiter=',')
writer.writerow(fieldnames)
birth_vertex = None

# The persistent pairs are one after the other
# First comes birth; immediately followed by death [Adhitya getting philosophical :P]

# Iterate over all the points in the persistent diagram
persistence_threshold_data = servermanager.Fetch(simplififedPersistenceDiagram)
# Get the number of persistent points and arcs
num_persistent_threshold_points = simplififedPersistenceDiagram.GetDataInformation().GetNumberOfPoints()

for index in range(num_persistent_threshold_points):
	vertex_id = persistence_threshold_data.GetPointData().GetArray('ttkVertexScalarField').GetValue(index)

	# If index is even, we are processing death; else just store attributes of birth
	# When death occurs, find persistence and moksha.
	if index & 1:
		death_vertex = vertex_id
		content = [birth_vertex, death_vertex]
		writer.writerow(content)

		# add to pairs dictionary
		# I store only one instance and not the death
		if (birth_vertex in scalars.keys()) and (death_vertex in scalars.keys()):
			birth_pairs[birth_vertex] = death_vertex

	else:
		birth_vertex = vertex_id

pairs_file.close()

# Write the persistence diagram after thresholding of only the super/sub level-set [for usage by TDA]
pairs_file_arguments = [tree_type, PAIRS_INFIX, file_name, CSV_EXTENSION]
pairs_file_path = get_output_path(file_path, pairs_file_arguments, folder_name = PERSISTENCE_FOLDER)

pairs_file = open(pairs_file_path, 'w')
fieldnames = ['dimension', 'Death', 'Birth']
writer = csv.writer(pairs_file, delimiter=',')
writer.writerow(fieldnames)
birth_vertex = None
birth_scalar = None
write_row = True

# The persistent pairs are one after the other
# First comes birth; immediately followed by death [Adhitya getting philosophical :P]

# Iterate over all the points in the persistent diagram
persistence_threshold_data = servermanager.Fetch(simplififedPersistenceDiagram)
# Get the number of persistent points and arcs
num_persistent_threshold_points = simplififedPersistenceDiagram.GetDataInformation().GetNumberOfPoints()

# Iterate across all points in diagram and write persistent pairs
for index in range(num_persistent_threshold_points):
	vertex_id = persistence_threshold_data.GetPointData().GetArray('ttkVertexScalarField').GetValue(index)
	try:
		vertex_scalar = scalars[vertex_id]
		if index & 1:
			death_vertex = vertex_id
			death_scalar = vertex_scalar
			# There exist values which are not present in the merge-tree
			if write_row:
				content = [0, round(death_scalar,4), round(birth_scalar,4)]
				writer.writerow(content)
			write_row = True
		else:
			birth_vertex = vertex_id
			birth_scalar = vertex_scalar
			write_row = True
	except:
		# This row contains a value not present in the merge-tree
		write_row = False
		pass

pairs_file.close()

# take screenshot of scalar field
screen_file_arguments = [tree_type, SCREENSHOT_INFIX, file_name, PNG_EXTENSION]
screen_file_path = get_output_path(file_path, screen_file_arguments, folder_name = SCREENSHOT_FOLDER)
SaveScreenshot(screen_file_path, magnification=1, quality=100, view=renderView1)


#######################################################################


# Start marching! :)

# simulation of simplices for nodes in the tree
sorted_nodes = list(reversed([v[0] for v in sorted(scalars.items(), key=lambda kv: (-kv[1], kv[0]))]))

# sort the nodes associated with each arc node
# the sorting is based on simulation of simplices
for node in sorted_nodes:
	arcs[node] = list(reversed(sorted(arcs[node], key= lambda x: -sorted_nodes.index(x))))

# retain only sub-arcs of the merge tree
# only retain the nodes that are indiced lower in the above ordering
# since we are working with the merge tree
# we should end up with each critical point in the tree except the global minima
# having a single sub-arc corresponding to a unique segment in the domain
for node in sorted_nodes:
	arcs[node][:] = [arc_node for arc_node in arcs[node] if sorted_nodes.index(arc_node) < sorted_nodes.index(node)]

# just for satisfaction
#for node in sorted_nodes:
#	if arcs[node]:
#		print (node, arcs[node])

# fetch the list of vertices and cells from the Domain
# use them for hashing
simplificationData = servermanager.Fetch(topologicalSimplification)
numTriangles = simplificationData.GetNumberOfCells()
numPoints = simplificationData.GetNumberOfPoints()

point_to_coordinates = {}
point_scalars = {}

# this is probably the most weirdest naming scheme I have ever used :P
triangle_index_to_points = defaultdict(list)
point_index_to_link_triangles = defaultdict(list)
triangle_points_to_triangle_index = {}

for triangle_index in xrange(numTriangles):
	# process every cell and get the identifiers of each point
	cell = simplificationData.GetCell(triangle_index)
	point1 = cell.GetPointId(0)
	point2 = cell.GetPointId(1)
	point3 = cell.GetPointId(2)

	point_index_to_link_triangles[point1].append(triangle_index)
	point_index_to_link_triangles[point2].append(triangle_index)
	point_index_to_link_triangles[point3].append(triangle_index)

	# get the coordinates of each point of the cell
	coordinate1 = simplificationData.GetPoint(point1)
	coordinate2 = simplificationData.GetPoint(point2)
	coordinate3 = simplificationData.GetPoint(point3)

	# hash the cell points together for comparison later
	cell_points = frozenset([point1, point2, point3])
	triangle_points_to_triangle_index[cell_points] = triangle_index

	triangle_index_to_points[triangle_index] = cell_points

	# hash the coordinates to PointID
	point_to_coordinates[point1] = coordinate1
	point_to_coordinates[point2] = coordinate2
	point_to_coordinates[point3] = coordinate3

# similarly iterate over all points in the surface and store their scalar values
for point_index in xrange(numPoints):
	point_scalar = simplificationData.GetPointData().GetArray('magnitude').GetValue(point_index)
	point_scalars[point_index] = point_scalar

# store the triangles of the segmented regions
segmented_triangle_regions = defaultdict(list)

# point: segment
processed_points = {}
# triangle: segment
processed_triangles = {}

# need to have a loop over the nodes here
for current_critical in sorted_nodes:
    # the arcs are already sorted
    # next_critical refers to the next lower critical point in the merge tree
    for next_critical in arcs[current_critical]:
        # find the lower and upper bounds of the current isoband
        lower_scalar_bound = point_scalars[next_critical]
        upper_scalar_bound = point_scalars[current_critical]

    	#print (current_critical, next_critical, upper_scalar_bound, lower_scalar_bound)

        # add the current critical point to the processing list
        current_processing_points = [current_critical]

        # process all points between current_critical and next_critical
        while (len(current_processing_points) > 0):
            # process the first element of the processing list
            simplex_point_identifier = current_processing_points.pop(0)

            # mark point as being processed
            processed_points[simplex_point_identifier] = current_critical

            # get all triangles in the the link of the simplex
            simplex_link_triangles = point_index_to_link_triangles[simplex_point_identifier]

            # iterate across each triangle
            for simplex_link_triangle in simplex_link_triangles:
                # make sure triangle has not been processed earlier
                if simplex_link_triangle not in processed_triangles:
                    isoband_value = ''
                    # get all points in the triangle
                    point_identifiers = triangle_index_to_points[simplex_link_triangle]
                    # iterate across each point in the triangle
                    # using the scalar of each point find the isoband value
                    for point_index in point_identifiers:
                        simplex_scalar = point_scalars[point_index]
                        isoband_value += calculate_cheap_isoband_index(simplex_scalar, lower_scalar_bound, upper_scalar_bound)

                    #if ((isoband_value in lower_intolerable_indices) or (isoband_value in upper_intolerable_indices)):
                    # check if the triangle has an isoband passing through it
                    if ((isoband_value in cheap_intolerable_indices)):
                        isoband = False
                    else:
                        processed_triangles[simplex_link_triangle] = current_critical
                        segmented_triangle_regions[current_critical].append(simplex_link_triangle)

                        # since the triangle has an isoband passing through it
                        for point_index in point_identifiers:
                            if point_index not in processed_points:
                                current_processing_points.append(point_index)
		# unnecessary, but just to be on the safer side
        del current_processing_points

#for current_critical in sorted_nodes:
#    print (current_critical, arcs[current_critical], len(segmented_triangle_regions[current_critical]), point_to_coordinates[current_critical])


segmentation_mapping = {}

# write the coordinates to a file
coordinates_file_arguments = [tree_type, COORDINATES_INFIX, file_name, CSV_EXTENSION]
coordinates_file_path = get_output_path(file_path, coordinates_file_arguments, folder_name = COORDINATES_FOLDER)
coordinates_file = open(coordinates_file_path, 'w')
fieldnames = ['identifier', 'x', 'y', 'z']
writer = csv.writer(coordinates_file, delimiter=',')
writer.writerow(fieldnames)

for critical_point in scalars.keys():
	[x, y, z] = point_to_coordinates[critical_point]
	writer.writerow([critical_point, x, y, z])

coordinates_file.close()

#*******************************************

make_hierarchy(file_name, file_path)

#*********************************************

# get the postorder mapping file
postorder_file_arguments = [tree_type, POSTORDER_INFIX, file_name, CSV_EXTENSION]
postorder_file_path = get_output_path(file_path, postorder_file_arguments, folder_name = POSTORDER_FOLDER)

# read the postorder mapping file
with open(postorder_file_path, 'r') as csvfile:
    csvfile.readline()
    spamreader = csv.reader(csvfile, delimiter=' ')
    for r in spamreader:
        row = r[0].split(',')
        node_id = int(row[0])
        order = int(row[1])
        segmentation_mapping[node_id] = order


# Here is what happens next
# triangle -> node -> postorder_index
# I iterate through every triangle and change it's segmentation identifier
# so that it corresponds to regions to that of the scalar field
# with which we are comparing
# Each node is given postorder index
# therefore each triangle should have such a mapping

for triangle in processed_triangles:
	segmentation = processed_triangles[triangle]
	#print ('segmentation', segmentation)
	# when you start merging saddles, it happens that the merged saddles
	# and their corresponding triangles don't have a mapping
	try:
		postorder_segmentation = segmentation_mapping[segmentation]
		#print ('postorder', postorder_segmentation)
		processed_triangles[triangle] = postorder_segmentation
	except:
		# this happens when the triangles associated with the mapped saddle gets merged
		processed_triangles[triangle] = None

# remove the triangles that were marked to be None
# these are triangles that were part of the segmentation of the initial merge tree
# but now not present because of the merging
processed_triangles = {k: v for k, v in processed_triangles.items() if v is not None}

# save the segmented scalar field
segmentation_file_arguments = [tree_type, SEGMENTATION_INFIX, file_name, VTK_EXTENSION]
segmentation_file_path = get_output_path(file_path, segmentation_file_arguments, folder_name = SEGMENTATION_FOLDER)
save_segmentation(point_to_coordinates, triangle_index_to_points, processed_triangles, segmentation_file_path)


# save the segmentation associated with each face

triangle_file_arguments = [tree_type, TRIANGLES_INFIX, file_name, CSV_EXTENSION]
triangle_file_path = get_output_path(file_path, triangle_file_arguments, folder_name = TRIANGLES_FOLDER)

triangle_file = open(triangle_file_path, 'w')
fieldnames = ['triangle', 'segmentation']
writer = csv.writer(triangle_file, delimiter=',')
writer.writerow(fieldnames)

# iterate through all triangles
# for triangles that do not have a segmentation, they are saved as 0
for triangle_index in xrange(numTriangles):
	try:
		content = [triangle_index, processed_triangles[triangle_index]]
		writer.writerow(content)
	except:
		content = [triangle_index, 0]
		writer.writerow(content)

triangle_file.close()


print (datetime.now() - startTime, 'Done! :)')

#os._exit(0)
