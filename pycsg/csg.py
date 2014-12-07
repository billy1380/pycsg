#
#  csg.py
#  pycsg
#
#  Created by William Shakour (billy1380) on 6 Dec 2014.
#  Copyright © 2014 SPACEHOPPER STUDIOS Ltd. All rights reserved.
#
import math

#Constructive Solid Geometry (CSG) is a modeling technique that uses Boolean
#operations like union and intersection to combine 3D solids. This library
#implements CSG operations on meshes elegantly and concisely using BSP trees,
#and is meant to serve as an easily understandable implementation of the
#algorithm. All edge cases involving overlapping coplanar polygons in both
#solids are correctly handled.
#
#Example usage:
#
#  var cube = CSG.cube();
#  var sphere = CSG.sphere({ radius: 1.3 });
#  var polygons = cube.subtract(sphere).toPolygons();
#
### Implementation Details
#
#All CSG operations are implemented in terms of two functions, `clipTo()` and
#`invert()`, which remove parts of a BSP tree inside another BSP tree and swap
#solid and empty space, respectively. To find the union of `a` and `b`, we
#want to remove everything in `a` inside `b` and everything in `b` inside `a`,
#then combine polygons from `a` and `b` into one solid:
#
#  a.clipTo(b);
#  b.clipTo(a);
#  a.build(b.allPolygons());
#
#The only tricky part is handling overlapping coplanar polygons in both trees.
#The code above keeps both copies, but we need to keep them in one tree and
#remove them in the other tree. To remove them from `b` we can clip the
#inverse of `b` against `a`. The code for union now looks like this:
#
#  a.clipTo(b);
#  b.clipTo(a);
#  b.invert();
#  b.clipTo(a);
#  b.invert();
#  a.build(b.allPolygons());
#
#Subtraction and intersection naturally follow from set operations. If
#union is `A | B`, subtraction is `A - B = ~(~A | B)` and intersection is
#`A & B = ~(~A | ~B)` where `~` is the complement operator.
#
### License
#
#Copyright (c) 2011 Evan Wallace (http:#madebyevan.com/), under the MIT license.

#Holds a binary space partition tree representing a 3D solid. Two solids can
#be combined using the `union()`, `subtract()`, and `intersect()` methods.
class Csg:
	polygons = []

	# Construct a CSG solid from a list of `CSG.Polygon` instances.
	@staticmethod
	def fromPolygons(polygons):
		csg = Csg()
		csg.polygons = polygons
		return csg

	def clone(self):
		csg = Csg()

		for p in self.polygons:
			csg.polygons.append(p.clone())

		return csg

	def toPolygons(self):
		return self.polygons

	#Return a CSG solid representing space in either this solid or in the
	# solid `csg`. Neither this solid nor the solid `csg` are modified.
	# 
	#     A.union(B)
	# 
	#     +-------+            +-------+
	#     |       |            |       |
	#     |   A   |            |       |
	#     |    +--+----+   =   |       +----+
	#     +----+--+    |       +----+       |
	#          |   B   |            |       |
	#          |       |            |       |
	#          +-------+            +-------+
	# 
	def union(self, csg):
		a = Node(self.clone().polygons)
		b = Node(csg.clone().polygons)
		a.clipTo(b)
		b.clipTo(a)
		b.invert()
		b.clipTo(a)
		b.invert()
		a.build(b.allPolygons())
		return Csg.fromPolygons(a.allPolygons())

	#Return a CSG solid representing space in this solid but not in the
	# solid `csg`. Neither this solid nor the solid `csg` are modified.
	# 
	#     A.subtract(B)
	# 
	#     +-------+            +-------+
	#     |       |            |       |
	#     |   A   |            |       |
	#     |    +--+----+   =   |    +--+
	#     +----+--+    |       +----+
	#          |   B   |
	#          |       |
	#          +-------+
	# 
	def subtract(self, csg):
		a = Node(self.clone().polygons)
		b = Node(csg.clone().polygons)
		a.invert()
		a.clipTo(b)
		b.clipTo(a)
		b.invert()
		b.clipTo(a)
		b.invert()
		a.build(b.allPolygons())
		a.invert();
		return Csg.fromPolygons(a.allPolygons())

	# Return a CSG solid representing space both this solid and in the
	# solid `csg`. Neither this solid nor the solid `csg` are modified.
	# 
	#     A.intersect(B)
	# 
	#     +-------+
	#     |       |
	#     |   A   |
	#     |    +--+----+   =   +--+
	#     +----+--+    |       +--+
	#          |   B   |
	#          |       |
	#          +-------+
	# 
	def intersect(self, csg):
		a = Node(self.clone().polygons)
		b = Node(csg.clone().polygons)
		a.invert()
		b.clipTo(a)
		b.invert()
		a.clipTo(b)
		b.clipTo(a)
		a.build(b.allPolygons())
		a.invert()
		return Csg.fromPolygons(a.allPolygons())

	# Return a CSG solid with solid and empty space switched. This solid is
	# not modified.
	def inverse(self):
		csg = self.clone()

		for polygon in csg.polygons:
			polygon.flip()

		return csg

	# Construct an axis-aligned solid cuboid. Optional parameters are `center` and
	# `radius`, which default to `[0, 0, 0]` and `[1, 1, 1]`. The radius can be
	# specified using a single number or a list of three numbers, one for each axis.
	# 
	# Example code:
	# 
	#	     var cube = CSG.cube({
	#	       center: [0, 0, 0],
	#	       radius: 1
	#	     });
	@staticmethod
	def cube(options = None):
		options = {} if options == None else options
		
		c = Vector([0, 0, 0] if 'center' not in options else options['center'])
		r = [1, 1, 1] if 'radius' not in options else options['radius'] if type(options['radius']) is list else [options['radius'], options['radius'], options['radius']]

		one = [[[0, 4, 6, 2], [-1, 0, 0]], [[1, 3, 7, 5], [+1, 0, 0]],[[0, 1, 5, 4], [0, -1, 0]], [[2, 6, 7, 3], [0, +1, 0]],[[0, 2, 3, 1], [0, 0, -1]], [[4, 5, 7, 6], [0, 0, +1]]]

		polygons = []

		for two in one:
			vertices = []
			for three in two[0]:
				pos = Vector(c.x + r[0] * (2 * (1 if (three & 1) != 0 else 0) - 1), c.y + r[1]
						* (2 * (1 if (three & 2) != 0 else 0) - 1), c.z + r[2] * (2 * (1 if (three & 4) != 0 else 0) - 1))
				v = Vertex(pos, Vector(two[1]))
				vertices.append(v)

			polygons.append(Polygon(vertices))

		return Csg.fromPolygons(polygons)


	@staticmethod
	def vertex(theta, phi, r, c, vertices):
		theta = theta * math.pi * 2
		phi = phi * math.pi
		direction = Vector(math.cos(theta) * math.sin(phi), math.cos(phi), math.sin(theta) * math.sin(phi))
		vertices.append(Vertex(c.plus(direction.times(r)), direction))

	# Construct a solid sphere. Optional parameters are `center`, `radius`,
	# `slices`, and `stacks`, which default to `[0, 0, 0]`, `1`, `16`, and `8`.
	# The `slices` and `stacks` parameters control the tessellation along the
	# longitude and latitude directions.
	# 
	# Example usage:
	# 
	#	     var sphere = CSG.sphere({
	#	       center: [0, 0, 0],
	#	       radius: 1,
	#	       slices: 16,
	#	       stacks: 8
	#	     });
	@staticmethod
	def sphere(options = None):
		options = {} if options == None else options

		c = Vector([0, 0, 0] if 'center' not in options else options['center'])
		r = 1 if 'radius' not in options else options['radius'] 
		slices = 16 if 'slices' not in options else options['slices'] 
		stacks = 8 if 'stacks' not in options else options['stacks']
		polygons = []

		for i in range(0, slices):
			for j in range(0, stacks):
				vertices = []
				Csg.vertex(i / slices, j / stacks, r, c, vertices)
				if j > 0:
					Csg.vertex((i + 1) / slices, j / stacks, r, c, vertices)
				if j < stacks - 1:
					Csg.vertex((i + 1) / slices, (j + 1) / stacks, r, c, vertices)
				Csg.vertex(i / slices, (j + 1) / stacks, r, c, vertices)
				polygons.append(Polygon(vertices))
		
		
		return Csg.fromPolygons(polygons)

	@staticmethod
	def point(stack, seg, normalBlend, s, axisX, axisY, axisZ, ray, r):
		angle = seg * math.pi * 2
		out = axisX.times(math.cos(angle)).plus(axisY.times(math.sin(angle)));
		pos = s.plus(ray.times(stack)).plus(out.times(r));
		normal = out.times(1 - math.fabs(normalBlend)).plus(axisZ.times(normalBlend));
		return Vertex(pos, normal)


	# Construct a solid cylinder. Optional parameters are `start`, `end`,
	# `radius`, and `slices`, which default to `[0, -1, 0]`, `[0, 1, 0]`, `1`, and
	# `16`. The `slices` parameter controls the tessellation.
	# 
	# Example usage:
	# 
	#	     var cylinder = CSG.cylinder({
	#	       start: [0, -1, 0],
	#	       end: [0, 1, 0],
	#	       radius: 1,
	#	       slices: 16
	#	     });
	@staticmethod
	def cylinder(options = None):
		options = {} if options == None else options
		
		s = Vector([0, -1, 0] if 'start' not in options else options['start'])
		e = Vector([0, 1, 0] if 'end' not in options else options['end'])
		ray = e.minus(s)
		r = 1 if 'radius' not in options else options['radius']
		slices = 16 if 'slices' not in options else options['slices']

		axisZ = ray.unit()
		isY = (math.fabs(axisZ.y) > 0.5)
		axisX = Vector(1 if isY else 0, 0 if isY else 1, 0).cross(axisZ).unit()
		axisY = axisX.cross(axisZ).unit()
		start = Vertex(s, axisZ.negated())
		end = Vertex(e, axisZ.unit())
		polygons = []

		for i in range(0, slices):
			t0 = i / slices
			t1 = (i + 1) / slices
			polygons.append(Polygon(start, Csg.point(0, t0, -1, s, axisX, axisY, axisZ, ray, r), Csg.point(0, t1, -1, s, axisX, axisY, axisZ, ray, r)))
			polygons.append(Polygon(Csg.point(0, t1, 0, s, axisX, axisY, axisZ, ray, r), Csg.point(0, t0, 0, s, axisX, axisY, axisZ, ray, r), Csg.point(1, t0, 0, s, axisX, axisY, axisZ, ray, r), Csg.point(1, t1, 0, s, axisX, axisY, axisZ, ray, r)))
			polygons.append(Polygon(end, Csg.point(1, t1, 1, s, axisX, axisY, axisZ, ray, r), Csg.point(1, t0, 1, s, axisX, axisY, axisZ, ray, r)))

		return Csg.fromPolygons(polygons)

#
#  Node.java
#  javacsg
#
#  Created by William Shakour (billy1380) on 6 Dec 2014.
#  Copyright © 2014 SPACEHOPPER STUDIOS Ltd. All rights reserved.
#

#Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
#by picking a polygon to split along. That polygon (and all other coplanar
#polygons) are added directly to that node and the other polygons are added to
#the front and/or back subtrees. This is not a leafy BSP tree since there is
#no distinction between internal and leaf nodes.
class Node:

	plane = None
	front = None
	back = None
	polygons = None
	
	def __init__(self, polygons = None):
		self.plane = None
		self.front = None
		self.back = None
		self.polygons = []

		if (polygons != None):
			self.build(polygons)

	def clone(self):
		node = Node()
		node.plane = None if self.plane == None else self.plane.clone()
		node.front = None if self.front == None else self.front.clone()
		node.back = None if self.back == None else self.back.clone()

		if self.polygons != None:
			node.polygons = []

			for polygon in self.polygons:
				node.polygons.append(polygon.clone())

		return node

	# Convert solid space to empty space and empty space to solid space.
	def invert(self):
		for i in range (0, len(self.polygons)):
			self.polygons[i].flip()

		self.plane.flip()

		if self.front != None:
			self.front.invert()
		if self.back != None:
			self.back.invert()

		temp = self.front
		self.front = self.back
		self.back = temp

	# Recursively remove all polygons in `polygons` that are inside this BSP
	# tree.
	def clipPolygons(self, polygons):
		if self.plane == None:
			return polygons[:]
		front = []
		back = []
		for i in range(0, len(polygons)):
			self.plane.splitPolygon(polygons[i], front, back, front, back)
			
		if self.front != None:
			front = self.front.clipPolygons(front)
		if self.back != None:
			back = self.back.clipPolygons(back)
		else:
			back = []
		
		return front + back

	# Remove all polygons in this BSP tree that are inside the other BSP tree
	# `bsp`.
	def clipTo(self, bsp):
		self.polygons = bsp.clipPolygons(self.polygons)
		
		if self.front != None:
			self.front.clipTo(bsp)

		if self.back != None:
			self.back.clipTo(bsp)

	# Return a list of all polygons in this BSP tree.
	def allPolygons(self):
		polygons = self.polygons[:]

		if self.front != None:
			polygons += self.front.allPolygons()

		if self.back != None:
			polygons += self.back.allPolygons()

		return polygons

	# Build a BSP tree out of `polygons`. When called on an existing tree, the
	# polygons are filtered down to the bottom of the tree and become new
	# nodes there. Each set of polygons is partitioned using the first polygon
	# (no heuristic is used to pick a good split).
	def build(self, polygons):
		if len(polygons) == 0:
			return
		
		if self.plane == None:
			self.plane = polygons[0].plane.clone()

		front = []
		back = []
		for i in range(0, len(polygons)):
			self.plane.splitPolygon(polygons[i], self.polygons, self.polygons, front, back)
		
		if len(front) > 0:
			if self.front == None:
				self.front = Node()
			
			self.front.build(front)
		
		if len(back) > 0:
			if self.back == None:
				self.back = Node()
			
			self.back.build(back)


#Represents a plane in 3D space.
class Plane:

	# `CSG.Plane.EPSILON` is the tolerance used by `splitPolygon()` to decide if a
	# point is on the plane.
	EPSILON = 1.0e-5
	
	COPLANAR = 0
	FRONT = 1
	BACK = 2
	SPANNING = 3

	def __init__(self, normal, w):
		self.normal = normal
		self.w = w

	@staticmethod
	def fromPoints(a, b, c):
		n = b.minus(a).cross(c.minus(a)).unit()
		return Plane(n, n.dot(a))

	def clone(self):
		return Plane(self.normal.clone(), self.w)

	def flip(self):
		self.normal = self.normal.negated()
		self.w = -self.w

	# Split `polygon` by this plane if needed, then put the polygon or polygon
	# fragments in the appropriate lists. Coplanar polygons go into either
	# `coplanarFront` or `coplanarBack` depending on their orientation with
	# respect to this plane. Polygons in front or in back of this plane go into
	# either `front` or `back`.
	def splitPolygon(self, polygon, coplanarFront, coplanarBack, front, back):
		# Classify each point as well as the entire polygon into one of the above
		# four classes.
		polygonType = 0
		types = []
		for i in range(0, len(polygon.vertices)):
			t = self.normal.dot(polygon.vertices[i].pos) - self.w
			currentType = Plane.BACK if (t < -Plane.EPSILON) else Plane.FRONT if (t > Plane.EPSILON) else Plane.COPLANAR
			polygonType |= currentType
			types.append(currentType)

		# Put the polygon in the correct list, splitting it when necessary.
		if (polygonType == Plane.COPLANAR):
			(coplanarFront if self.normal.dot(polygon.plane.normal) > 0 else coplanarBack).append(polygon)
		elif (polygonType == Plane.FRONT):
			front.append(polygon)
		if (polygonType == Plane.BACK):
			back.append(polygon)
		if (polygonType == Plane.SPANNING):
			f = []
			b = []
			for i in range(0, len(polygon.vertices)):
				j = (i + 1) % len(polygon.vertices)
				ti = types[i]
				tj = types[j]
				vi = polygon.vertices[i]
				vj = polygon.vertices[j]
				if (ti != Plane.BACK):
					f.append(vi)
				if (ti != Plane.FRONT):
					b.append(vi.clone() if ti != Plane.BACK else vi)
				if ((ti | tj) == Plane.SPANNING):
					t = (self.w - self.normal.dot(vi.pos)) / self.normal.dot(vj.pos.minus(vi.pos))
					v = vi.interpolate(vj, t)
					f.append(v)
					b.append(v.clone())
			
			if (len(f) >= 3):
				front.append(Polygon(f, polygon.shared))
			if (len(b) >= 3):
				back.append(Polygon(b, polygon.shared))
				
#Represents a convex polygon. The vertices used to initialize a polygon must
#be coplanar and form a convex loop. They do not have to be `CSG.Vertex`
#instances but they must behave similarly (duck typing can be used for
#customization).
#
#Each convex polygon has a `shared` property, which is shared between all
#polygons that are clones of each other or were split from the same polygon.
#This can be used to define per-polygon properties (such as surface color).
class Polygon:

	def __init__(self, vertices, shared = None):
		self.vertices = vertices
		self.shared = shared
		self.plane = Plane.fromPoints(vertices[0].pos, vertices[1].pos, vertices[2].pos)

	def clone(self):
		vertexClone = []

		for vertex in self.vertices:
			vertexClone.append(vertex)

		return Polygon(vertexClone, self.shared)

	def flip(self):
		self.vertices.reverse()

		for vertex in self.vertices:
			vertex.flip()

		self.plane.flip()

#Represents a 3D vector.
#
#Example usage:
#
#  CSG.Vector(1, 2, 3);
#  CSG.Vector([1, 2, 3]);
#  CSG.Vector({ x: 1, y: 2, z: 3 });
class Vector:

	x = 0.0
	y = 0.0
	z = 0.0

	def __init__(self, x, y = None, z = None):
		if y == None and z == None:
			self.x = x[0]
			self.y = x[1]
			self.z = x[2]
		else: 
			self.x = x
			self.y = y
			self.z = z
		
	def clone(self):
		return Vector(self.x, self.y, self.z)

	def times(self, a):
		return Vector(self.x * a, self.y * a, self.z * a)

	def dividedBy(self, a):
		return Vector(self.x / a, self.y / a, self.z / a)

	def plus(self, a):
		return Vector(self.x + a.x, self.y + a.y, self.z + a.z)

	def minus(self, a):
		return Vector(self.x - a.x, self.y - a.y, self.z - a.z)

	def unit(self):
		return self.dividedBy(self.length())

	def length(self):
		return math.sqrt(self.dot(self))

	def dot(self, a):
		return self.x * a.x + self.y * a.y + self.z * a.z

	def cross(self, a):
		return Vector(self.y * a.z - self.z * a.y, self.z * a.x - self.x
				* a.z, self.x * a.y - self.y * a.x)

	def negated(self):
		return Vector(-self.x, -self.y, -self.z)

	def lerp(self, a, t):
		return self.plus(a.minus(self).times(t))

#Represents a vertex of a polygon. Use your own vertex class instead of this
#one to provide additional features like texture coordinates and vertex
#colors. Custom vertex classes need to provide a `pos` property and `clone()`,
#`flip()`, and `interpolate()` methods that behave analogous to the ones
#defined by `CSG.Vertex`. This class provides `normal` so convenience
#functions like `CSG.sphere()` can return a smooth vertex normal, but `normal`
#is not used anywhere else.
class Vertex:

	pos = None
	normal = None

	def __init__(self, pos, normal):
		self.pos = pos.clone()
		self.normal = normal.clone()

	def clone(self):
		return Vertex(self.pos.clone(), self.normal.clone())

	# Invert all orientation-specific data (e.g. vertex normal). Called when the
	# orientation of a polygon is flipped.
	def flip(self):
		self.normal = self.normal.negated()

	# Create a vertex between this vertex and `other` by linearly
	# interpolating all properties using a parameter of `t`. Subclasses should
	# override this to interpolate additional properties.
	def interpolate(self, other, t):
		return Vertex(self.pos.lerp(other.pos, t), self.normal.lerp(other.normal, t));
