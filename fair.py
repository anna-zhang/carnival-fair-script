
import bpy
from math import sqrt, pi, sin, cos
import numpy as np
from bpy.props import BoolProperty, FloatProperty, IntProperty, PointerProperty
from bpy.types import Operator, Panel, PropertyGroup
import bmesh

# returns the coordinates of the vertex along the edge connecting vertex_1 and vertex_2, according to a weight value (0.5 weight value gets halfway point)
def interpolate_edge_vertex(vertex_1, vertex_2, weight):
    vertex_x = vertex_1[0] + (vertex_2[0] - vertex_1[0]) * weight 
    vertex_y = vertex_1[1] + (vertex_2[1] - vertex_1[1]) * weight 
    vertex_z = vertex_1[2] + (vertex_2[2] - vertex_1[2]) * weight 
    return (vertex_x, vertex_y, vertex_z)

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
        cart_mesh = bpy.data.meshes.new('cart')
        cart_mesh.from_pydata(self.vertices, self.edges, self.faces)
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
        self.posts["edge_indices"]["rings"] = wheel_rings
            
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
                
    def create_carts(self):
        num_carts = self.num_carts # number of carts to create
        cart_collection = bpy.data.collections.new('cart_collection')
        bpy.context.scene.collection.children.link(cart_collection)
        # add cart object to scene collection
        ring_edge_indices = self.posts["edge_indices"]["rings"] # get the indices of the edges of the rings
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
        
 
# test
# new_wheel = Wheel((1.0, 2.0, 1.0), 15, 3) # smaller wheel
new_wheel = Wheel((1.0, 2.0, 1.0), 15, 10) # larger wheel

ferris_wheel_mesh = bpy.data.meshes.new('ferris_wheel')
ferris_wheel_mesh.from_pydata(new_wheel.vertices, new_wheel.edges, new_wheel.faces)
ferris_wheel_mesh.update()
# make ferris wheel object from ferris wheel mesh
ferris_wheel_object = bpy.data.objects.new('ferris_wheel_object', ferris_wheel_mesh)
# make ferris wheel collection
ferris_wheel_collection = bpy.data.collections.new('ferris_wheel_collection')
bpy.context.scene.collection.children.link(ferris_wheel_collection)
# add ferris wheel object to scene collection
ferris_wheel_collection.objects.link(ferris_wheel_object)