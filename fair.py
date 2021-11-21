
import bpy
from math import sqrt, pi, sin, cos
import numpy as np
from bpy.props import BoolProperty, FloatProperty, IntProperty, PointerProperty
from bpy.types import Operator, Panel, PropertyGroup
import bmesh


class Wheel:
    def __init__(self, center, num_carts, size):
        self.center = center # center of wheel, user-specified
        self.num_carts = 15 # number of carts, user-specified
        self.size = 3.0 # size of wheel, user-specified
        
        self.vertices = [] # vertices of the ferris wheel
        self.edges = [] # edges of the ferris wheel
        self.faces = [] # faces of the ferris wheel
        self.wheel1 = {"vertex_indices": {}, "edge_indices": {}}
        self.wheel2 = {"vertex_indices": {}, "edge_indices": {}}
        self.posts = {"vertex_indices": {}, "edge_indices": {}}
        
        self.create_wheel((self.center[0], self.center[1] - (0.15 * self.size), self.center[2]), self.size, self.num_carts, self.wheel1) # 2.0 is the radius of a ring holding a cart
        self.create_wheel((self.center[0], self.center[1] + (0.15 * self.size), self.center[2]), self.size, self.num_carts, self.wheel2)
        self.connect_wheels()
        self.create_posts(self.wheel1, 1)
        self.create_posts(self.wheel2, 2)
        self.create_internal_support()
    
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
            
                
 
# test
new_wheel = Wheel((1.0, 2.0, 1.0), 3, 8)

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