
import bpy
from math import sqrt, pi, sin, cos
import numpy as np
from bpy.props import BoolProperty, FloatProperty, IntProperty, PointerProperty
from bpy.types import Operator, Panel, PropertyGroup
import bmesh

test_vertices = []
test_edges = []

class Circle:
    def __init__(self, center, radius, vertices, edges):
        self.center = center # the location of the center of the vertex
        self.radius = radius # the radius of the circle
        self.vertices = vertices # vertices forming the circle
        self.edges = edges # edges connecting vertices on the circle

class Cylinder:
    def __init__(self, v1, v2, radius):
        self.endpoint1 = v1
        self.endpoint2 = v2
        self.radius = radius # the radius of the cylinder
        self.vertices = [] # vertices of the cylinder
        self.edges = [] # edges of the cylinder
        self.vertex_indices = {}
        self.edge_indices = {}
        
        self.create_cylinder(self.endpoint1, self.endpoint2, self.radius) # test create cylinder
        
    
    def create_cylinder(self, v1, v2, radius):
        num_cuts = 20 # the two circles must have the same number of vertices
        num_vertices = len(self.vertices) # current number of vertices in cylinder object
        num_edges = len(self.edges) # current number of edges in cylinder object
        # create circle endcaps
        circle_1 = create_circle(np.asarray(v1) - np.asarray(v2), v1, radius, num_cuts, 1) # create circle1
        circle_2 = create_circle(np.asarray(v1) - np.asarray(v2), v2, radius, num_cuts, 2)  # create circle2
        self.vertices = self.vertices + circle_1.vertices + circle_2.vertices
        circle_1.edges = [(edge[0] + num_vertices, edge[1] + num_vertices) for edge in circle_1.edges] # update edges with vertex indices in the cylinder vertices array
        circle_2.edges = [(edge[0] + num_vertices + len(circle_1.vertices), edge[1] + num_vertices + len(circle_1.vertices)) for edge in circle_1.edges] # update edges with vertex indices in the cylinder vertices array
        self.edges = self.edges + circle_1.edges + circle_2.edges
        
        # remember vertex indices of the circle endcaps
        self.vertex_indices["circle_1"] = [i + num_vertices for i in range(len(circle_1.vertices))] # remember the vertex indices of the vertices making up circle1
        self.vertex_indices["circle_2"] = [i + num_vertices + len(circle_1.vertices) for i in range(len(circle_2.vertices))] # remember the vertex indices of the vertices making up circle2
        
        # remember edge indices of the circle endcaps
        self.edge_indices["circle_1"] = [i + num_edges for i in range(len(circle_1.edges))] # remember the edge indices of the edges making up circle1
        self.edge_indices["circle_2"] = [i + num_edges + len(circle_1.edges) for i in range(len(circle_2.edges))] # remember the edge indices of the edges making up circle2
        
        # connect circles
        num_edges = len(self.edges) # update current number of edges
        for i in range(num_cuts):
            self.edges.append((self.vertex_indices["circle_1"][i], self.vertex_indices["circle_2"][i]))
        self.edge_indices["connect_circles"] = [i + num_edges for i in range(num_cuts)] # remember the edge indices of the edges connecting the two circles
        

# takes a list of cylinders and creates an object to add to "collection"
def cylinders_to_obj(cylinders, collection):
    vertices = [] # list of all vertices in the object
    edges = [] # list of all edges in the object
    for i in range(len(cylinders)):
        cylinder = cylinders[i]
        num_vertices = len(vertices) # current number of vertices
        transformed_edges = [(edge[0] + num_vertices, edge[1] + num_vertices) for edge in cylinder.edges] # edges with new vertex indices corresponding to index in final object vertices array
        for vertex in cylinder.vertices:
            vertices.append(vertex)
        for edge in transformed_edges:
            edges.append(edge)
    
    # create one object out of all of the cylinders
    cylinder_mesh = bpy.data.meshes.new('cylinders') # create Mesh object
    cylinder_mesh.from_pydata(vertices, edges, [])
    cylinder_mesh.update()
    bm = bmesh.new() # create BMesh object
    bm.from_mesh(cylinder_mesh) # take Mesh object and turn to BMesh
    bmesh.ops.contextual_create(bm, geom=bm.edges) # create faces from edges
    bm.to_mesh(cylinder_mesh) # take BMesh object and turn to Mesh
    cylinder_mesh.update()
    # make cylinder object from cylinder mesh
    cylinder_object = bpy.data.objects.new('cylinder_object', cylinder_mesh)
    # add it to the cylinder collection
    collection.objects.link(cylinder_object)
    
    
def get_all_values(dictionary, int_values): # given a dictionary (with potentially multiple layers), get all integer values using recursion; int_values is initially an empty list
    for key, value in dictionary.items():
        if type(value) is dict:
            get_all_values(value, int_values)
        elif type(value) is list:
            for integer in value:
                int_values.append(integer)
        elif isinstance(value, int):
            int_values.append(value)
    return int_values
  

# returns the coordinates of the vertex along the edge connecting vertex_1 and vertex_2, according to a weight value (0.5 weight value gets halfway point)
def interpolate_edge_vertex(vertex_1, vertex_2, weight):
    vertex_x = vertex_1[0] + (vertex_2[0] - vertex_1[0]) * weight 
    vertex_y = vertex_1[1] + (vertex_2[1] - vertex_1[1]) * weight 
    vertex_z = vertex_1[2] + (vertex_2[2] - vertex_1[2]) * weight 
    return (vertex_x, vertex_y, vertex_z)

def create_circle(normal, center, radius, num_vertices, circle_num): # determine plane orthogonal to a given normal vector containing "center" point and create a circle of "radius" with vertices on that plane around the center point; returns a list of vertices defining the circle
        theta = (2 * pi) / num_vertices # calculate angle of each slice of the circle
        circle_center = np.asarray(center) # the center of the circle, given
        # compute the plane with normal "normal" and that includes the point "center"
        plane_A = normal[0]
        plane_B = normal[1]
        plane_C = normal[2]
        plane_D = (normal[0] * center[0] + normal[1] * center[1] + normal[2] * center[2])
        plane = (plane_A, plane_B, plane_C, plane_D) # stores (a, b, c, d) for plane ax + by + cz = d orthogonal to a given vector "normal"
        # compute two other points that lay on this plane 
        point1 = (0.0, 0.0, 0.0) # holder for point1 on the plane
        point2 = (0.0, 0.0, 0.0) # holder for point2 on the plane
        # find two points on the plane
        if plane_C != 0:
            # compute point1
            point1_a = (center[0] + 1)
            point1_b = center[1]
            point1_c = center[2] + (-1.0 * plane_A / plane_C)
            point1 = (point1_a, point1_b, point1_c)
            # compute point2
            point2_a = center[0]
            point2_b = center[1] + 1
            point2_c = center[2] + (-1.0 * plane_B / plane_C)
            point2 = (point2_a, point2_b, point2_c)
        elif plane_A != 0:
            # compute point1
            point1_a = center[0] + (-1.0 * plane_B / plane_A)
            point1_b = center[1] + 1
            point1_c = center[2]
            point1 = (point1_a, point1_b, point1_c)
            # compute point2
            point2_a = center[0] + (-1.0 * plane_C / plane_A)
            point2_b = center[1]
            point2_c = center[2] + 1
            point2 = (point2_a, point2_b, point2_c)
        else: # plane_B is not 0
            # compute point1
            point1_a = center[0] + 1
            point1_b = center[1] + (-1.0 * plane_A / plane_B)
            point1_c = center[2]
            point1 = (point1_a, point1_b, point1_c)
            # compute point2
            point2_a = center[0]
            point2_b = center[1] + (-1.0 * plane_C / plane_B)
            point2_c = center[2] + 1
            point2 = (point2_a, point2_b, point2_c)
        
        vector_1 = np.asarray(point1) - np.asarray(center) # vector 1 defining plane containing circle
        print("vector_1")
        print(vector_1)
        vector_2 = np.asarray(point2) - np.asarray(center) # vector 2 defining plane containing circle
        print("vector_2")
        print(vector_2)
        v1v2_cross = np.cross(vector_1, vector_2) # get cross product of two vectors to get a vector normal to the plane with the circle
        print("v1v2_cross")
        print(v1v2_cross)
        v1v2_cross_normalized = v1v2_cross / np.linalg.norm(v1v2_cross) # normalize the vector that's normal to the plane containing the circle
        print("v1v2_cross_normalized")
        print(v1v2_cross_normalized)
        u = vector_1 / np.linalg.norm(vector_1) # get unit vector to serve as "x-axis" of the plane
        v = np.cross(u, v1v2_cross_normalized) # get vector to serve as "y-axis" of the plane
        v = v / np.linalg.norm(v) # normalize v to get the unit "y-axis" of the plane
        print("u")
        print(u)
        print("v")
        print(v)
        print("dot product: " + str(np.dot(u, v))) # make sure dot product of u and v = 0 so the axes are perpendicular
        
        # create vertices forming the circle
        circle_vertices = [] # save the vertices lying on the circle
        for i in range(0, num_vertices):
            vertex = tuple(circle_center + radius * cos(theta * i) * u + radius * sin(theta * i) * v) 
            circle_vertices.append(vertex)
            print("circle_vertex: " + str(vertex))
            
        # create edges between adjacent circle vertices
        circle_edges = [] # save the edges creating the circle
        for i in range(0, num_vertices):
            if i < (num_vertices - 1): # all vertices besides the last one
                circle_edges.append((i, i + 1)) # create an edge between this vertex and the next one
            else:
                circle_edges.append((i, 0)) # create an edge between the last vertex and the first one to close the circle
        
        return Circle(center, radius, circle_vertices, circle_edges)
            


class Booth: 
    def __init__(self, base_shape, top_style, top_border, purpose, width, length):
        self.base_shape = base_shape # square, rectangle, hexagon
        self.top_style = top_style # flat, pointy
        self.top_border = top_border # flat, mini triangle flags, mini half hexagon flags
        self.purpose = purpose # food, games
        self.width = width # helps define the size of the booth (where the structural poles supporting it go)
        self.length = length # helps define the size of the booth (where the structural poles supporting it go)
        self.height = 3 # total height of the booth
        self.center = (20.0, 0.0, 0.0) # test
        
        self.vertices = []
        self.edges = []
        self.faces = []
        
        self.top = {"vertex_indices": {}, "edge_indices": {}}
        self.bottom = {"vertex_indices": {}, "edge_indices": {}}
        self.poles = {"vertex_indices": {}, "edge_indices": {}}
        
        self.create_top(self.base_shape, self.top_style, self.top_border, self.width, self.length, self.height, self.center) # test
        self.create_poles(self.base_shape, self.width, self.length, self.height) # test
        
    
    def create_top(self, shape, style, border, width, length, height, center):
        top_vertices = []
        top_edges = []
        num_vertices = len(self.vertices) # current number of vertices in the object
        num_edges = len(self.edges) # current number of edges in the object
        if style == "pointy":
            # pointy booth top
            top_v1 = (center[0] + 0.5 * width, center[1] + 0.5 * length, center[2] + 0.5 * height)
            # top_v1_index = num_vertices
            top_v2 = (center[0] - 0.5 * width, center[1] + 0.5 * length, center[2] + 0.5 * height)
            # top_v2_index = num_vertices + 1
            top_v3 = (center[0] - 0.5 * width, center[1] - 0.5 * length, center[2] + 0.5 * height)
            # top_v3_index = num_vertices + 2
            top_v4 = (center[0] + 0.5 * width, center[1] - 0.5 * length, center[2] + 0.5 * height)
            # top_v4_index = num_vertices + 3
            top_pointy_vertex = (center[0], center[1], center[2] + 0.75 * height)
            # top_pointy_vertex_index = num_vertices + 4
            top_vertices = [top_pointy_vertex, top_v1, top_v2, top_v3, top_v4]
            top_edges = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4), (4, 1)]
            
        # flat booth top
        elif shape == "rectangle":
            # the very top rectangular prism of the booth
            top_v1 = (center[0] + 0.5 * width, center[1] + 0.5 * length, center[2] + 0.5 * height)
            top_v2 = (center[0] - 0.5 * width, center[1] + 0.5 * length, center[2] + 0.5 * height)
            top_v3 = (center[0] - 0.5 * width, center[1] - 0.5 * length, center[2] + 0.5 * height)
            top_v4 = (center[0] + 0.5 * width, center[1] - 0.5 * length, center[2] + 0.5 * height)
            top_rect_vertices = [top_v1, top_v2, top_v3, top_v4]
            
            # the bottom of the top rectangular prism of the booth
            top_v5 = (center[0] + 0.5 * width, center[1] + 0.5 * length, center[2] + 0.45 * height)
            top_v6 = (center[0] - 0.5 * width, center[1] + 0.5 * length, center[2] + 0.45 * height)
            top_v7 = (center[0] - 0.5 * width, center[1] - 0.5 * length, center[2] + 0.45 * height)
            top_v8 = (center[0] + 0.5 * width, center[1] - 0.5 * length, center[2] + 0.45 * height)
            bottom_rect_vertices = [top_v5, top_v6, top_v7, top_v8]
            
            top_vertices = top_rect_vertices + bottom_rect_vertices # all vertices forming the top of the booth
            self.top["vertex_indices"]["top_rect"] = [i + num_vertices for i in range(len(top_rect_vertices))] # remember the vertex indices of the vertices making up the top rectangle
            self.top["vertex_indices"]["bottom_rect"] = [i + num_vertices + len(top_rect_vertices) for i in range(len(bottom_rect_vertices))] # remember the vertex indices of the vertices making up the bottom rectangle
            
            top_rect_edges = [(0, 1), (1, 2), (2, 3), (3, 0)] # form the top of the top rectangular prism of the booth
            bottom_rect_edges = [(4, 5), (5, 6), (6, 7), (7, 0)] # form the bottom of the top rectangular prism of the booth
            connect_rect_edges = [(0, 4), (1, 5), (2, 6), (3, 7)] # connect the top and bottom of the top rectangular prism
            
            top_edges = top_rect_edges + bottom_rect_edges + connect_rect_edges
            
            
        elif shape == "hexagon":
             # the very top hexagonal prism of the booth
            top_v1 = (center[0] + 0.5 * width, center[1], center[2] + 0.5 * height)
            top_v2 = (center[0] + 0.35 * width, center[1] + 0.5 * length, center[2] + 0.5 * height)
            top_v3 = (center[0] - 0.35 * width, center[1] + 0.5 * length, center[2] + 0.5 * height)
            top_v4 = (center[0] - 0.5 * width, center[1], center[2] + 0.5 * height)
            top_v5 = (center[0] - 0.35 * width, center[1] - 0.5 * length, center[2] + 0.5 * height)
            top_v6 = (center[0] + 0.35 * width, center[1] - 0.5 * length, center[2] + 0.5 * height)
            top_hexagon_vertices = [top_v1, top_v2, top_v3, top_v4, top_v5, top_v6]
            # the bottom of the top hexagonal prism of the booth
            top_v7 = (center[0] + 0.5 * width, center[1], center[2] + 0.30 * height)
            top_v8 = (center[0] + 0.35 * width, center[1] + 0.5 * length, center[2] + 0.30 * height)
            top_v9 = (center[0] - 0.35 * width, center[1] + 0.5 * length, center[2] + 0.30 * height)
            top_v10 = (center[0] - 0.5 * width, center[1], center[2] + 0.30 * height)
            top_v11 = (center[0] - 0.35 * width, center[1] - 0.5 * length, center[2] + 0.30 * height)
            top_v12 = (center[0] + 0.35 * width, center[1] - 0.5 * length, center[2] + 0.30 * height)
            bottom_hexagon_vertices = [top_v7, top_v8, top_v9, top_v10, top_v11, top_v12]
            top_vertices = top_hexagon_vertices + bottom_hexagon_vertices
            self.top["vertex_indices"]["top_hexagon"] = [i + num_vertices for i in range(len(top_hexagon_vertices))] # remember the vertex indices of the vertices making up the top hexagon
            self.top["vertex_indices"]["bottom_hexagon"] = [i + num_vertices + len(top_hexagon_vertices) for i in range(len(bottom_hexagon_vertices))] # remember the vertex indices of the vertices making up the bottom hexagon
            
            top_hexagon_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)] # form the top of the hexagonal prism of the booth
            bottom_hexagon_edges = [(6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 6)] # form the bottom of the hexagonal prism of the booth
            connect_hexagon_edges = [(0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)] # connect the top and bottom of the hexagonal prism
            top_edges = top_hexagon_edges + bottom_hexagon_edges + connect_hexagon_edges
            self.top["edge_indices"]["top_hexagon"] = [i + num_edges for i in range(len(top_hexagon_edges))] # remember the edge indices of the edges making up the top hexagon
            self.top["edge_indices"]["bottom_hexagon"] = [i + num_edges + len(top_hexagon_edges) for i in range(len(bottom_hexagon_edges))] # remember the edge indices of the edges making up the bottom hexagon
            self.top["edge_indices"]["connect_hexagons"] = [i + num_edges + len(top_hexagon_edges) + len(bottom_hexagon_edges) for i in range(len(connect_hexagon_edges))] # remember the edge indices of the edges connecting the top and bottom hexagons

        top_edges = [(edge[0] + num_vertices, edge[1] + num_vertices) for edge in top_edges]
        
        self.vertices = self.vertices + top_vertices
        self.edges = self.edges + top_edges
    
    
    def create_poles(self, shape, width, length, height):
        poles = []
        if shape == "rectangle":
            # get the vertices and edges forming the surface of the top of the booth that the pole connects to
            top_vertices = self.top["vertex_indices"]["bottom_rect"]
            
            # radius of the pole is 1/12 the length of the shortest side of the rectangle
            radius = 0
            if self.length > self.width:
                radius = 1/12 * self.width
            else:
                radius = 1/12 * self.length
            
            # create cylinder poles near each of the vertices in top_vertices
            for i in range(len(top_vertices)):
                top_vertex = self.vertices[top_vertices[i]]
                if top_vertex[0] < self.center[0]:
                    # vertex has a smaller x than the center of the booth
                    x_displacement = radius
                else:
                    x_displacement = -1.0 * radius
                
                if top_vertex[1] < self.center[1]:
                    # vertex has a smaller y than the center of the booth
                    y_displacement = radius
                else:
                    y_displacement = -1.0 * radius
                
                
                top_vertex = (top_vertex[0] + x_displacement, top_vertex[1] + y_displacement, top_vertex[2])
                bottom_vertex = (top_vertex[0], top_vertex[1], top_vertex[2] - height)
                new_pole = Cylinder(top_vertex, bottom_vertex, radius) # create the cylinder pole
                poles.append(new_pole)
            
        elif shape == "hexagon":
            # get the vertices and edges forming the surface of the top of the booth that the pole connects to
            top_vertices = self.top["vertex_indices"]["bottom_hexagon"]
            print("top_vertices")
            print(top_vertices)
            
            # radius of the pole is 1/12 the length of the shortest side of the rectangle that fits inside the hexagon
            radius = 0
            if self.length > self.width:
                radius = 1/12 * self.width
            else:
                radius = 1/12 * self.length
            
            # create cylinder poles near each of the rectangle vertices in top_vertices
            rect_vertex_indices = [1, 2, 4, 5] # in the array of bottom_hexagon vertices, these are the indices of the vertices forming a rectangle within the hexagonal shape
            for i in rect_vertex_indices:
                top_vertex = self.vertices[top_vertices[i]]
                if top_vertex[0] < self.center[0]:
                    # vertex has a smaller x than the center of the booth
                    x_displacement = radius
                else:
                    x_displacement = -1.0 * radius
                
                if top_vertex[1] < self.center[1]:
                    # vertex has a smaller y than the center of the booth
                    y_displacement = radius
                else:
                    y_displacement = -1.0 * radius
                
                top_vertex = (top_vertex[0] + x_displacement, top_vertex[1] + y_displacement, top_vertex[2])
                bottom_vertex = (top_vertex[0], top_vertex[1], top_vertex[2] - height)
                new_pole = Cylinder(top_vertex, bottom_vertex, radius) # create the cylinder pole
                poles.append(new_pole)

        # make cylinder collection
        cylinder_collection = bpy.data.collections.new('cylinder_collection')
        bpy.context.scene.collection.children.link(cylinder_collection)
        cylinders_to_obj(poles, cylinder_collection) # create object out of list of cylinder poles

        # TEST CREATE BOOTH    
        booth_top_mesh = bpy.data.meshes.new('booth') # create Mesh object
        booth_top_mesh.from_pydata(self.vertices, self.edges, [])
        booth_top_mesh.update()
        bm = bmesh.new() # create BMesh object
        bm.from_mesh(booth_top_mesh) # take Mesh object and turn to BMesh
        bmesh.ops.contextual_create(bm, geom=bm.edges)
        bm.to_mesh(booth_top_mesh) # take BMesh object and turn to Mesh
        booth_top_mesh.update()
        
        # make booth object from booth mesh
        booth_object = bpy.data.objects.new('booth_object', booth_top_mesh)
        # make booth collection
        booth_collection = bpy.data.collections.new('booth_collection')
        bpy.context.scene.collection.children.link(booth_collection)
        # add booth object to scene collection
        booth_collection.objects.link(booth_object)


class Cart:
    def __init__(self, center, width, height):
        self.center = center
        self.width = width # half the size of every ring on the ferris wheel (the radius of the circle extending along that ring)
        self.height = 0.9 * height # leave some room between carts
        self.obj = None # hold the Blender object
        
        self.vertices = [] # vertices of the ferris wheel cart
        self.edges = [] # edges of the ferris wheel cart
        self.faces = [] # faces of the ferris wheel cart
        self.top = {"vertex_indices": {}, "edge_indices": {}}
        self.base = {"vertex_indices": {}, "edge_indices": {}}
        self.pole = {"vertex_indices": {}, "edge_indices": {}}
        
        self.create_cart(self.center, self.width, self.height)
        
        
    def create_horizontal_circle(self, center, radius, num_vertices, circle_dict, circle_num): # add vertices of circle with center "center" and radius "radius"
            theta = (2 * pi) / num_vertices # calculate angle of each slice of the circle
            circle_center = np.asarray(center)
            
            self.vertices.append(center) # add vertex for circle center
            circle_dict["vertex_indices"]["center_" + str(circle_num)] = len(self.vertices) - 1
            
            u = np.asarray((1.0, 0.0, 0.0)) # get unit vector of x-direction to serve as "x-axis" of the plane
            v = np.asarray((0.0, 1.0, 0.0)) # get up unit vector of y-direction to serve as "y-axis" of the plane
            print("dot product: " + str(np.dot(u, v))) # make sure dot product of u and v = 0 so the axes are perpendicular
            num_vertices_total = len(self.vertices) # get current # of vertices
            circle_vertices = [] # save indices of the vertices lying on the circle
            for i in range(0, num_vertices):
                vertex = tuple(circle_center + radius * cos(theta * i) * u + radius * sin(theta * i) * v) 
                self.vertices.append(vertex)
                print("circle_vertex: " + str(vertex))
                circle_vertices.append(num_vertices_total + i) # save index in vertices list of the vertex being added
            circle_dict["vertex_indices"]["outer_circle_" + str(circle_num)] = circle_vertices
            
            
            # create edges between adjacent outer circle vertices
            circle_edges = [] # save indices of the edges creating the circle in the edges list
            total_edges = len(self.edges) # the total number of edges
            for i in range(0, num_vertices):
                if i < (num_vertices - 1): # all vertices besides the last one
                    self.edges.append([circle_vertices[i], circle_vertices[i + 1]]) # create an edge between this vertex and the next one
                    circle_edges.append(total_edges + i) # store what index that newly added edge is within the entire edges list
                else:
                    self.edges.append([circle_vertices[i], circle_vertices[0]]) # create an edge between the last vertex and the first one to close the circle
                    circle_edges.append(total_edges + i) # store what index that newly added edge is within the entire edges list
            circle_dict["edge_indices"]["outer_circle_" + str(circle_num)] = circle_edges
            
    def create_cart_obj(self):
        cart_mesh = bpy.data.meshes.new('cart') # create Mesh object
        cart_mesh.from_pydata(self.vertices, self.edges, self.faces)
        cart_mesh.update()
        bm = bmesh.new() # create BMesh object
        bm.from_mesh(cart_mesh) # take Mesh object and turn to BMesh
        bmesh.ops.contextual_create(bm, geom=bm.edges) # create faces from edges
        print("bm.faces") # TEST
        print([face for face in bm.faces]) # TEST: figure out what face to delete
        faces_to_delete = [bm.faces[7], bm.faces[14]] # bm.faces[7] is the bottom face of the top hood, bm.faces[14] is the top of the basket
        bmesh.ops.delete(bm, geom=faces_to_delete, context='FACES') # remove the unnecessary faces
        bm.to_mesh(cart_mesh) # take BMesh object and turn to Mesh
        cart_mesh.update()
        # make cart object from cart mesh
        cart_object = bpy.data.objects.new('cart_object', cart_mesh)
        # make cart collection
        self.obj = cart_object
            
              
    def create_cart(self, center, width, height):
        # create top cylinder
        self.create_horizontal_circle((center[0], center[1], center[2]), width, 6, self.top, 1) # create top hexagon of the top cylinder
        self.create_horizontal_circle((center[0], center[1], center[2] -(0.05 * height)), width, 6, self.top, 2) # create top hexagon of the top cylinder
        top_circle1_vertices = self.top["vertex_indices"]["outer_circle_1"] # get list of numbers corresponding to the indices of the vertices of the top hexagon of the top cylinder
        top_circle2_vertices = self.top["vertex_indices"]["outer_circle_2"] # get list of numbers corresponding to the indices of the vertices of the bottom hexagon of the top cylinder
        top_connect_edges = [] # hold edges connecting the hexagons
        total_edges = len(self.edges) # the total number of edges
        for i in range(6): # go through every vertex of the hexagon circle
            self.edges.append([top_circle1_vertices[i], top_circle2_vertices[i]]) # create an edge between the corresponding vertices on the two hexagons
            top_connect_edges.append(total_edges + i) # store what index that newly added edge is within the entire edges list 
        self.top["edge_indices"]["top"] = top_connect_edges 
        
        # create base basket
        self.create_horizontal_circle((center[0], center[1], center[2] - (0.4 * height)), width, 6, self.base, 1) # create top hexagon of the basket
        self.create_horizontal_circle((center[0], center[1], center[2] - height), width / 1.5, 6, self.base, 2) # create bottom hexagon of the basket
        base_circle1_vertices = self.base["vertex_indices"]["outer_circle_1"] # get list of numbers corresponding to the indices of the vertices of the top hexagon of the basket
        base_circle2_vertices = self.base["vertex_indices"]["outer_circle_2"] # get list of numbers corresponding to the indices of the vertices of the bottom hexagon of the basket
        base_connect_edges = [] # hold edges connecting the hexagons
        total_edges = len(self.edges) # the total number of edges
        for i in range(6): # go through every vertex of the hexagon circle
            self.edges.append([base_circle1_vertices[i], base_circle2_vertices[i]]) # create an edge between the corresponding vertices on the two hexagons
            base_connect_edges.append(total_edges + i) # store what index that newly added edge is within the entire edges list
        self.base["edge_indices"]["basket"] = base_connect_edges 
            
        # create the cylinder attached to the ring holding the cart up
        top_center = self.vertices[self.top["vertex_indices"]["center_1"]] # get the vertex coordinates of the center of the cylinder top
        bottom_center = self.vertices[self.base["vertex_indices"]["center_2"]] # get the vertex coordinates of the center of the basket bottom
        self.create_horizontal_circle(top_center, width / 16, 12, self.pole, 1) # create top circle of the pole
        self.create_horizontal_circle(bottom_center, width / 16, 12, self.pole, 2) # create bottom of the pole
        pole_circle1_vertices = self.pole["vertex_indices"]["outer_circle_1"] # get list of numbers corresponding to the indices of the vertices of the top circle of the pole
        pole_circle2_vertices = self.pole["vertex_indices"]["outer_circle_2"] # get list of numbers corresponding to the indices of the vertices of the bottom circle of the pole
        pole_connect_edges = [] # hold edges connecting the top and bottom circles of the pole
        total_edges = len(self.edges) # the total number of edges
        for i in range(12): # go through every vertex of the pole circle
            self.edges.append([pole_circle1_vertices[i], pole_circle2_vertices[i]]) # create an edge between the corresponding vertices on the two circles
            pole_connect_edges.append(total_edges + i) # store what index that newly added edge is within the entire edges list
        self.pole["edge_indices"]["pole"] = pole_connect_edges 
        
        self.create_cart_obj() # create Blender object


class Wheel:
    def __init__(self, center, num_carts, size):
        self.center = center # center of wheel, user-specified
        self.num_carts = num_carts # number of carts, user-specified
        self.size = size # size of wheel, user-specified
        
        self.vertices = [] # vertices of the ferris wheel
        self.edges = [] # edges of the ferris wheel
        self.faces = [] # faces of the ferris wheel
        self.wheel1 = {"vertex_indices": {}, "edge_indices": {}}
        self.wheel2 = {"vertex_indices": {}, "edge_indices": {}}
        self.posts = {"vertex_indices": {}, "edge_indices": {}}
        self.cart_width = 0.20 * self.size * 0.5 # size of each ring (which is 0.30 * self.size) and then divide by 2
        self.cart_height = (self.size / (num_carts / 4)) # approximate the height each cart can take up
        
        self.create_wheel((self.center[0], self.center[1] - (0.10 * self.size), self.center[2]), self.size, self.num_carts, self.wheel1) # 2.0 is the radius of a ring holding a cart
        self.create_wheel((self.center[0], self.center[1] + (0.10 * self.size), self.center[2]), self.size, self.num_carts, self.wheel2)
        self.connect_wheels()
        self.create_posts(self.wheel1, 1)
        self.create_posts(self.wheel2, 2)
        self.create_internal_support()
        self.create_carts()
        self.create_wheel_obj()
    
    def create_wheel(self, center, radius, num_vertices, wheel_dict): # add vertices of circle with center "center" and radius "radius"
            theta = (2 * pi) / num_vertices # calculate angle of each slice of the circle
            circle_center = np.asarray(center)
            
            self.vertices.append(center) # add vertex for circle center
            wheel_dict["vertex_indices"]["center"] = len(self.vertices) - 1
            
            u = np.asarray((1.0, 0.0, 0.0)) # get unit vector of x-direction to serve as "x-axis" of the plane
            v = np.asarray((0.0, 0.0, 1.0)) # get up unit vector of z-direction to serve as "y-axis" of the plane
            print("dot product: " + str(np.dot(u, v))) # make sure dot product of u and v = 0 so the axes are perpendicular
            num_vertices_total = len(self.vertices) # get current # of vertices
            circle_vertices = [] # save indices of the vertices lying on the circle
            for i in range(0, num_vertices):
                vertex = tuple(circle_center + radius * cos(theta * i) * u + radius * sin(theta * i) * v) 
                self.vertices.append(vertex)
                print("circle_vertex: " + str(vertex))
                circle_vertices.append(num_vertices_total + i) # save index in vertices list of the vertex being added
            wheel_dict["vertex_indices"]["outer_circle"] = circle_vertices
            
            # create edges between adjacent outer circle vertices
            circle_edges = [] # save indices of the edges creating the circle in the edges list
            total_edges = len(self.edges) # the total number of edges
            for i in range(0, num_vertices):
                if i < (num_vertices - 1): # all vertices besides the last one
                    self.edges.append([circle_vertices[i], circle_vertices[i + 1]]) # create an edge between this vertex and the next one
                    circle_edges.append(total_edges + i) # store what index that newly added edge is within the entire edges list
                else:
                    self.edges.append([circle_vertices[i], circle_vertices[0]]) # create an edge between the last vertex and the first one to close the circle
                    circle_edges.append(total_edges + i) # store what index that newly added edge is within the entire edges list
            wheel_dict["edge_indices"]["outer_circle"] = circle_edges
            
            
            # connect vertices on wheel to center of circle
            spoke_edges = [] # save indices of the edges creating the wheel spokes
            total_edges = len(self.edges) # the total number of edges
            center_vertex = wheel_dict["vertex_indices"]["center"]
            for i in range(0, len(circle_vertices)):
                self.edges.append([center_vertex, circle_vertices[i]])
                spoke_edges.append(total_edges + i) # store what index that newly added spoke edge is within the entire edges list
            wheel_dict["edge_indices"]["spoke"] = spoke_edges
            
    def connect_wheels(self):
        wheel1_circle_vertices = self.wheel1["vertex_indices"]["outer_circle"] # get list of numbers corresponding to the indices of the outer circle vertices in the vertices array
        wheel2_circle_vertices = self.wheel2["vertex_indices"]["outer_circle"] # get list of numbers corresponding to the indices of the outer circle vertices in the vertices array
        wheel_rings = [] # hold edges connecting the wheels (the rings the carts are on)
        total_edges = len(self.edges) # the total number of edges
        for i in range(self.num_carts):
            self.edges.append([wheel1_circle_vertices[i], wheel2_circle_vertices[i]]) # create an edge between the corresponding vertices on the two wheels
            wheel_rings.append(total_edges + i) # store what index that newly added edge is within the entire edges list
        self.wheel1["edge_indices"]["rings"] = wheel_rings
        self.wheel2["edge_indices"]["rings"] = wheel_rings
            
    def create_posts(self, wheel, wheel_num):
        # wheel_num is either 1 or 2, depending on which wheel of this ferris wheel is specified to create posts for
        # create the posts supporting this wheel
        wheel_center_ground = self.vertices[wheel["vertex_indices"]["center"]]
        wheel_center_ground = (wheel_center_ground[0], wheel_center_ground[1], -1.25 * self.size) # point on ground that has 1.25 * radius of the wheel
        wheel_post1_vertex = (wheel_center_ground[0] - (0.6 * self.size), wheel_center_ground[1], wheel_center_ground[2]) # the left post supporting the wheel, the point on ground that has 1.5 * radius of the wheel and that is to the left of the center vertex
        wheel_post2_vertex = (wheel_center_ground[0] + (0.6 * self.size), wheel_center_ground[1], wheel_center_ground[2]) # the right post supporting the wheel, the point on ground that has 1.5 * radius of the wheel and that is to the right of the center vertex
        # these vertices are support points for the wheel
        total_vertices = len(self.vertices)
        self.vertices.append(wheel_post1_vertex)
        self.vertices.append(wheel_post2_vertex)
        self.posts["vertex_indices"]["wheel" + str(wheel_num)] = [total_vertices, total_vertices + 1]
        
        # create edges between the center of the wheel and the post support points on the ground
        total_edges = len(self.edges)
        self.edges.append([wheel["vertex_indices"]["center"], total_vertices])
        self.edges.append([wheel["vertex_indices"]["center"], total_vertices + 1])
        self.posts["edge_indices"]["wheel" + str(wheel_num)] = [total_edges, total_edges + 1]
        
    def create_internal_support(self):
        # create bar that goes through the two center vertices
        wheel1_center_vi = self.wheel1["vertex_indices"]["center"]
        wheel2_center_vi = self.wheel2["vertex_indices"]["center"]
        self.edges.append([wheel1_center_vi, wheel2_center_vi]) # create an edge between the centers of the two wheels
        self.posts["edge_indices"]["center_bar"] = len(self.edges) - 1
        
        # create crossing bars 
        # determine number of crosses to make by getting half of self.size, rounding to the nearest integer
        num_crosses = int(self.size / 3)

        # iterate through every pair of spokes and create crosses between them
        num_spoke_pairs = self.num_carts
        wheel1_spokes = self.wheel1["edge_indices"]["spoke"]
        wheel2_spokes = self.wheel2["edge_indices"]["spoke"]
        
        self.wheel1["vertex_indices"]["crosses"] = {} # dictionary with keys as spoke_# and value as array of vertex indices for the cross points on that spoke
        self.wheel2["vertex_indices"]["crosses"] = {} # dictionary with keys as spoke_# and value as array of vertex indices for the cross points on that spoke
        
        self.wheel1["edge_indices"]["crosses"] = {} # dictionary with keys as spoke_# and value as array of edge indices that spoke is a part of
        self.wheel2["edge_indices"]["crosses"] = {} # dictionary with keys as spoke_# and value as array of edge indices that spoke is a part of
        
        # add vertices on the spokes for crossing bars
        for i in range(num_spoke_pairs):
            wheel1_spoke_v1_index = (self.edges[wheel1_spokes[i]])[0] 
            wheel1_spoke_v1 = self.vertices[wheel1_spoke_v1_index] # the center of the wheel, one end of the spoke
            
            wheel1_spoke_v2_index = (self.edges[wheel1_spokes[i]])[1]
            wheel1_spoke_v2 = self.vertices[wheel1_spoke_v2_index] # the end of the spoke near the cart
            
            wheel2_spoke_v1_index = (self.edges[wheel2_spokes[i]])[0] 
            wheel2_spoke_v1 = self.vertices[wheel2_spoke_v1_index] # the center of the wheel, one end of the spoke
            
            wheel2_spoke_v2_index = (self.edges[wheel2_spokes[i]])[1]
            wheel2_spoke_v2 = self.vertices[wheel2_spoke_v2_index] # the end of the spoke near the cart
            
            self.wheel1["vertex_indices"]["crosses"]["spoke_" + str(i)] = [wheel1_spoke_v1_index] # remember the end of the spoke near the center of the wheel
            self.wheel2["vertex_indices"]["crosses"]["spoke_" + str(i)] = [wheel2_spoke_v1_index] # remember the end of the spoke near the center of the wheel
            
            if num_crosses > 1:
                for j in range(1, num_crosses):
                    wheel1_cross_vertex = interpolate_edge_vertex(wheel1_spoke_v1, wheel1_spoke_v2, j / num_crosses)
                    wheel2_cross_vertex = interpolate_edge_vertex(wheel2_spoke_v1, wheel2_spoke_v2, j / num_crosses)
                    # these vertices are the cross points on the spokes
                    total_vertices = len(self.vertices)
                    self.vertices.append(wheel1_cross_vertex)
                    self.vertices.append(wheel2_cross_vertex)
                    
                    self.wheel1["vertex_indices"]["crosses"]["spoke_" + str(i)].append(total_vertices) # remember the index of the vertex on the spoke
                    self.wheel2["vertex_indices"]["crosses"]["spoke_" + str(i)].append(total_vertices + 1) # remember the index of the vertex on the spoke
            
            else:
                # create a cross vertex 0.3 in
                 wheel1_cross_vertex = interpolate_edge_vertex(wheel1_spoke_v1, wheel1_spoke_v2, 0.3)
                 wheel2_cross_vertex = interpolate_edge_vertex(wheel2_spoke_v1, wheel2_spoke_v2, 0.3)
                 # these vertices are the cross points on the spokes
                 total_vertices = len(self.vertices)
                 self.vertices.append(wheel1_cross_vertex)
                 self.vertices.append(wheel2_cross_vertex)
                 self.wheel1["vertex_indices"]["crosses"]["spoke_" + str(i)].append(total_vertices) # remember the index of the vertex on the spoke
                 self.wheel2["vertex_indices"]["crosses"]["spoke_" + str(i)].append(total_vertices + 1) # remember the index of the vertex on the spoke
                 
            final_num_crosses = len(self.wheel1["vertex_indices"]["crosses"]["spoke_" + str(i)]) - 1 # there are at least two vertices
            
            # create the crossing bar edges
            self.wheel1["edge_indices"]["crosses"]["spoke_" + str(i)] = [] # remember the edges added
            self.wheel2["edge_indices"]["crosses"]["spoke_" + str(i)] = [] # remember the edges added
            for j in range(final_num_crosses):
                # get the four vertices for the cross
                wheel1_cross_v1 = self.wheel1["vertex_indices"]["crosses"]["spoke_" + str(i)][j]
                wheel1_cross_v2 = self.wheel1["vertex_indices"]["crosses"]["spoke_" + str(i)][j + 1]
                wheel2_cross_v1 = self.wheel2["vertex_indices"]["crosses"]["spoke_" + str(i)][j]
                wheel2_cross_v2 = self.wheel2["vertex_indices"]["crosses"]["spoke_" + str(i)][j + 1]
                total_edges = len(self.edges)
                # create the X
                self.edges.append([wheel1_cross_v1, wheel2_cross_v2])
                self.edges.append([wheel2_cross_v1, wheel1_cross_v2])
                self.edges.append([wheel1_cross_v2, wheel2_cross_v2]) # connect the outer cross vertices
                self.wheel1["edge_indices"]["crosses"]["spoke_" + str(i)].extend([total_edges, total_edges + 1, total_edges + 2])
                self.wheel2["edge_indices"]["crosses"]["spoke_" + str(i)].extend([total_edges, total_edges + 1, total_edges + 2])
                
        # create circles connecting the cross vertices on each wheel for internal support
        num_cross_vertices = len(self.wheel1["vertex_indices"]["crosses"]["spoke_0"]) # every spoke has the same number of cross vertices, so just use the first spoke as reference   
        # connect corresponding cross vertices on the other spokes of that wheel with edges
        # iterate through the number of cross vertices on a spoke (there are at least two, since the vertex at the end of the spoke toward circle center also counts as one)
        self.wheel1["edge_indices"]["internal_circles"] = {}
        self.wheel2["edge_indices"]["internal_circles"] = {}
        for i in range(1, num_cross_vertices):
            wheel1_circle_edges = [] # list of edge indices forming this circle on wheel1
            wheel2_circle_edges = [] # list of edge indices forming this circle on wheel2
            for j in range(num_spoke_pairs):
                if j == (num_spoke_pairs - 1):
                    # the last spoke, connect the cross vertex on the last spoke to the corresponding cross vertex on the first spoke
                    wheel1_cross_spoke1_vertex = self.wheel1["vertex_indices"]["crosses"]["spoke_" + str(j)][i] # cross vertex on this spoke of wheel1
                    wheel1_cross_spoke2_vertex = self.wheel1["vertex_indices"]["crosses"]["spoke_" + str(0)][i] # corresponding cross vertex on the first spoke of wheel1 to close the loop
                    
                    wheel2_cross_spoke1_vertex = self.wheel2["vertex_indices"]["crosses"]["spoke_" + str(j)][i] # cross vertex on this spoke of wheel2
                    wheel2_cross_spoke2_vertex = self.wheel2["vertex_indices"]["crosses"]["spoke_" + str(0)][i] # corresponding cross vertex on the first spoke of wheel2 to close the loop
                    total_edges = len(self.edges)
                    self.edges.append([wheel1_cross_spoke1_vertex, wheel1_cross_spoke2_vertex])
                    self.edges.append([wheel2_cross_spoke1_vertex, wheel2_cross_spoke2_vertex])
                    wheel1_circle_edges.append(total_edges) # remember index of newly added edge
                    wheel2_circle_edges.append(total_edges + 1) # remember index of newly added edge   
                else:
                    # connect the cross vertex on this spoke with the corresponding cross vertex on the next spoke
                    wheel1_cross_spoke1_vertex = self.wheel1["vertex_indices"]["crosses"]["spoke_" + str(j)][i] # cross vertex on this spoke of wheel1
                    wheel1_cross_spoke2_vertex = self.wheel1["vertex_indices"]["crosses"]["spoke_" + str(j + 1)][i] # corresponding cross vertex on the next spoke of wheel1
                    
                    wheel2_cross_spoke1_vertex = self.wheel2["vertex_indices"]["crosses"]["spoke_" + str(j)][i] # cross vertex on this spoke of wheel2
                    wheel2_cross_spoke2_vertex = self.wheel2["vertex_indices"]["crosses"]["spoke_" + str(j + 1)][i] # corresponding cross vertex on the next spoke of wheel2
                    total_edges = len(self.edges)
                    self.edges.append([wheel1_cross_spoke1_vertex, wheel1_cross_spoke2_vertex])
                    self.edges.append([wheel2_cross_spoke1_vertex, wheel2_cross_spoke2_vertex])
                    wheel1_circle_edges.append(total_edges) # remember index of newly added edge
                    wheel2_circle_edges.append(total_edges + 1) # remember index of newly added edge
                
                # save the list of edge indices 
                self.wheel1["edge_indices"]["internal_circles"]["circle_" + str(i - 1)] = wheel1_circle_edges # i - 1 so that the key names start with circle_0
                self.wheel2["edge_indices"]["internal_circles"]["circle_" + str(i - 1)] = wheel2_circle_edges # i - 1 so that the key names start with circle_0
                          
        
    def create_carts(self):
        num_carts = self.num_carts # number of carts to create
    
        # error checking for height of cart
        wheel_center_vertex = self.vertices[self.wheel1["vertex_indices"]["center"]] # the vertex  for the center of the first wheel
        post_vertex = self.vertices[self.posts["vertex_indices"]["wheel1"][0]] # get vertex of post bottom on wheel1
        post_vertical_height = abs(wheel_center_vertex[2] - post_vertex[2])
        if (self.cart_height > (post_vertical_height - self.size)): # make sure the bottom cart doesn't go below the support structure posts when it's hanging down
            self.cart_height = 0.80 * (post_vertical_height - self.size) # adjust height of cart so that it doesn't go below the support structure posts   
        
        # create collection to hold carts
        cart_collection = bpy.data.collections.new('cart_collection')
        bpy.context.scene.collection.children.link(cart_collection)
        # add cart object to scene collection
        ring_edge_indices = self.wheel1["edge_indices"]["rings"] # get the indices of the edges of the rings
        print("ring_edge_indices" + str(ring_edge_indices))
        for i in range(num_carts): 
            # get center of ring
            ring_vertices = self.edges[ring_edge_indices[i]] # the indices of the two vertices forming the ring
            ring_v1 = self.vertices[ring_vertices[0]] # the first ring vertex
            ring_v2 = self.vertices[ring_vertices[1]] # the second ring vertex
            print("ring_v1" + str(ring_v1))
            print("ring_v2" + str(ring_v2))
            
            # compute center of cart top
            top_center = (ring_v1[0] + 0.5 * (ring_v2[0] - ring_v1[0]), ring_v1[1] + 0.5 * (ring_v2[1] - ring_v1[1]), ring_v1[2] + 0.5 * (ring_v2[2] - ring_v1[2]))
            new_cart = Cart(top_center, self.cart_width, self.cart_height)
            cart_collection.objects.link(new_cart.obj) # add the newly created cart object to the cart_collection
    
    def create_wheel_obj(self):
        # create wheel object from a list of vertices and edges
        wheel_cylinders = []
        post_cylinders = []
        
        # compute thickness of the bars forming the ferris wheel depending on user-specified wheel size
        wheel_size = self.size 
        radius = wheel_size * 0.01
        
        # every edge is a cylinder
        wheel1_edge_indices = get_all_values(self.wheel1["edge_indices"], []) # get all edge indices of the edges that make up wheel1
        wheel2_edge_indices = get_all_values(self.wheel2["edge_indices"], []) # get all edge indices of the edges that make up wheel2
        wheel_edge_indices = list(set(wheel1_edge_indices + wheel2_edge_indices)) # all edge indices of the wheel
        
        # iterate through each edge in the wheel and turn it into a cylinder
        for edge_index in wheel_edge_indices:
            edge = self.edges[edge_index] # get the edge
            v1 = self.vertices[edge[0]] # get vertex1 of the edge
            v2 = self.vertices[edge[1]] # get vertex2 of the edge
            new_cylinder = Cylinder(v1, v2, radius)
            wheel_cylinders.append(new_cylinder)
        
        post_edge_indices = get_all_values(self.posts["edge_indices"], []) # get all edge indices of the supporting posts of the wheel
        # iterate through each edge making up the posts and turn it into a cylinder
        for edge_index in post_edge_indices:
            edge = self.edges[edge_index] # get the edge
            v1 = self.vertices[edge[0]] # get vertex1 of the edge
            v2 = self.vertices[edge[1]] # get vertex2 of the edge
            new_cylinder = Cylinder(v1, v2, radius)
            post_cylinders.append(new_cylinder)

        # create a single object from all of these cylinders
        # make ferris wheel collection
        ferris_wheel_collection = bpy.data.collections.new('ferris_wheel_collection')
        bpy.context.scene.collection.children.link(ferris_wheel_collection)
        cylinders_to_obj(wheel_cylinders, ferris_wheel_collection) # create the ferris wheel
        cylinders_to_obj(post_cylinders, ferris_wheel_collection) # create the ferris wheel posts
        
 
# test
# new_wheel = Wheel((1.0, 2.0, 1.0), 15, 3) # smaller wheel
new_wheel = Wheel((1.0, 2.0, 1.0), 12, 10) # larger wheel

# flat_booth = Booth("rectangle", "flat", "flat", "food", 4, 8) # flat top booth
hex_booth = Booth("hexagon", "flat", "flat", "food", 4, 8) # flat top booth