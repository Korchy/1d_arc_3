# Nikita Akimov
# interplanety@interplanety.org
#
# GitHub
#    https://github.com/Korchy/1d_arc_3
import math

import bmesh
import bpy
from bpy.props import BoolProperty, EnumProperty, FloatProperty, IntProperty
from bpy.types import Operator, Panel, Scene
from bpy.utils import register_class, unregister_class
from mathutils import Matrix, Vector
from math import isclose, log, pi

bl_info = {
    "name": "3 Points Arc",
    "description": "Creates arc from 2 edges (3 starting points).",
    "author": "Nikita Akimov, Paul Kotelevets",
    "version": (1, 1, 0),
    "blender": (2, 79, 0),
    "location": "View3D > Tool panel > 1D > 3 Points Arc",
    "doc_url": "https://github.com/Korchy/1d_arc_3",
    "tracker_url": "https://github.com/Korchy/1d_arc_3",
    "category": "All"
}


# MAIN CLASS

class Arc3:

    @classmethod
    def arc3(cls, context, obj, points=5, edge_length=1.25, base_points=32, base_points_rate=2.0,
             arc_mode='POINTS', invert_direction=False):
        # create arc from 2 edges (3 vertices)
        obj = obj if obj else context.active_object
        # current mode
        mode = obj.mode
        if obj.mode == 'EDIT':
            bpy.ops.object.mode_set(mode='OBJECT')
        # switch to vertex selection mode
        obj_world_matrix = obj.matrix_world.copy()
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        bm.verts.ensure_lookup_table()
        bm.edges.ensure_lookup_table()
        # process all vertex selection islands - try to form list of chunks, each chunk with 3 vertices
        selected_vertices = [vertex for vertex in bm.verts if vertex.select]
        _l = len(selected_vertices)
        _i = 0
        while selected_vertices:
            v0 = v1 = v2 = None
            v0 = selected_vertices[0]
            if len(v0.link_edges) == 1 and v0.link_edges[0].select == False:
                # selected vertex at the end of the edges loop (ex. after Step Extrude)
                v1 = v0.link_edges[0].other_vert(v0)
                # special case: if v0 and v1 placed on the same place (have the same coordinates)
                #   (after Step Extrude interruption)
                if cls._is_vectors_close(v0.co, v1.co):
                    bm.verts.remove(v0)
                    selected_vertices.remove(v0)
                    v0 = v1
                    v1 = v0.link_edges[0].other_vert(v0) if v0.link_edges else None
                if v1:
                    v2_edge = next((edge for edge in v1.link_edges if edge != v0.link_edges[0]), None)
                    v2 = v2_edge.other_vert(v1) if v2_edge else None
            elif len(v0.link_edges) == 2 and v0.link_edges[0].select and v0.link_edges[1].select:
                # vertex selected between two selected edges
                v1 = v0
                v0 = v1.link_edges[0].other_vert(v1)
                v2 = v1.link_edges[1].other_vert(v1)
            elif len(v0.link_edges) == 2 and (v0.link_edges[0].select or v0.link_edges[1].select):
                # corner vertex on two selected edges (left or right)
                selected_edge = next((_edge for _edge in v0.link_edges if _edge.select), None)
                v1 = selected_edge.other_vert(v0)
                selected_edge = next((_edge for _edge in v1.link_edges if _edge.select and _edge != selected_edge), None)
                if selected_edge:
                    v2 = selected_edge.other_vert(v1)
            # if all three vertices found - create 3-point arc by this 3 vertices
            if v0 and v1 and v2:
                cls._arc_by_3_vertices(
                    bm=bm,
                    v0=v0,
                    v1=v1,
                    v2=v2,
                    world_matrix=obj_world_matrix,
                    points=points,
                    edge_length=edge_length,
                    base_points=base_points,
                    base_points_rate=base_points_rate,
                    arc_mode=arc_mode,
                    invert_direction=invert_direction
                )
            # remove from selected vertices list
            for v in (_v for _v in (v0, v1, v2) if _v in selected_vertices):
                selected_vertices.remove(v)
            # alarm break
            _i += 1
            if _i > _l:
                print('processing selected vertices overflow err exit')
                break
        # save changed data to mesh
        bm.to_mesh(obj.data)
        bm.free()
        # return mode back
        bpy.ops.object.mode_set(mode=mode)

    @classmethod
    def butch_clean(cls, context, obj, points=5, edge_length=1.25, base_points=32, base_points_rate=2.0,
                    arc_mode='BASE_32'):
        # butch processing of raw (butch) arcs
        obj = obj if obj else context.active_object
        # current mode
        mode = obj.mode
        if obj.mode == 'EDIT':
            bpy.ops.object.mode_set(mode='OBJECT')
        # switch to vertex selection mode
        obj_world_matrix = obj.matrix_world.copy()
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        bm.verts.ensure_lookup_table()
        bm.edges.ensure_lookup_table()
        # select edges mode
        bm.select_mode = {'EDGE'}
        bm.select_flush_mode()
        # process all vertex selection islands - try to form list of butch arcs to process them
        # try to get bounding vertices - which have at least one deselected linked edge
        bounding_vertices = [vertex for vertex in bm.verts
                             if vertex.select and (len([_edge for _edge in vertex.link_edges if not _edge.select]) > 0
                                                   or len(vertex.link_edges) == 1)
                             ]
        # print([v.index for v in bounding_vertices])
        # from each bounding vertex try to walk to another bounding vertex to form a pair of bounding vertices for arc
        # pairs => list of pair bounding vertices and all other vertices belongs to this pair
        pairs = []  # [(bounding_vert_1, bounding_vert_2, [vertex, vertex, ...]), ...]
        checked_vertices = set()
        selected_vertices = [_vert for _vert in bm.verts if _vert.select]
        for vertex in bounding_vertices:
            if vertex not in checked_vertices:
                path = cls._possible_path(
                    start_point=vertex,
                    possible_points=selected_vertices,
                    possible_endpoints=bounding_vertices
                )
                if path:
                    # pair found
                    pairs.append((path[0], path[-1], []))
                    checked_vertices.update((path[0], path[-1]))
                else:
                    # get bounded vertex but can't find a pair to it
                    print('WARNING: can\'t find a pair for vertex:', vertex.index)
        # print(pairs)
        # for each non-bounding vertex check for what pair of bounding vertices it belongs
        selected_non_bounding_vertices = [_vert for _vert in bm.verts
                                          if _vert.select and _vert not in bounding_vertices]
        for vertex in selected_non_bounding_vertices:
            path = cls._possible_path(
                start_point=vertex,
                possible_points=selected_vertices,
                possible_endpoints=bounding_vertices
            )
            if path:
                # found path to one of bounding points - add to pair with this bounding point
                pair = next((_pair for _pair in pairs if path[-1] in {_pair[0], _pair[1]}), None)
                if pair:
                    pair[2].append(vertex)
            else:
                # can't find a path to on of bounding vertices
                print('WARNING: can\'t find a bounding vertex for vertex:', vertex.index)
        # for _p in pairs:
        #     print(_p)

        # remove pairs which has no belonging vertices between its bounding vertices
        broken_pairs = [_pair for _pair in pairs if not _pair[2]]
        if broken_pairs:
            print('WARNING: groups consists with only 2 vertices', broken_pairs)
        pairs = [_pair for _pair in pairs if _pair[2]]
        # for each pair process its belonging vertices to find the third point for building the circle by 3 points
        #   we need to save list of belonging vertices to remove them after creating arc
        #   the p1 will be removed automatically when creating arc, so the list shouldn't include it
        circles_points = []     # [(p0, p1, p2, [vertex, vertex, ...]), ...]
        for pair in pairs:
            p0 = pair[0]
            p2 = pair[1]
            # from vertices belonging to this pair find one which is mostly equidistant from bounding vertices
            #   length of the vector from first bounding point to this vertex should be close to the length of the
            #   vector from the second bounding point to this vertex
            p1 = min(((_vertex, abs((p0.co - _vertex.co).length - (p2.co - _vertex.co).length)) for _vertex in pair[2]),
                     key=lambda _item: _item[1])
            circles_points.append((p0, p1[0] if p1 else None, p2, set(pair[2])-{p1[0] if p1 else None}))
        # for _c in circles_points:
        #     print(_c)
        # now we have several blocks, each with 3 points, for creating circles
        for circle_data in circles_points:
            # create arc
            # print(circle_data)
            cls._arc_by_3_vertices(
                bm=bm,
                v0=circle_data[0],
                v1=circle_data[1],
                v2=circle_data[2],
                world_matrix=obj_world_matrix,
                points=points,
                edge_length=edge_length,
                base_points=base_points,
                base_points_rate=base_points_rate,
                arc_mode=arc_mode
            )
            # remove source vertices
            for vertex in circle_data[3]:
                bm.verts.remove(vertex)
        # save changed data to mesh
        bm.to_mesh(obj.data)
        bm.free()
        # return mode back
        bpy.ops.object.mode_set(mode=mode)

    @classmethod
    def _arc_by_3_vertices(cls, bm, v0, v1, v2, world_matrix, points=5, edge_length=1.25, base_points=32,
                           base_points_rate=2.0, arc_mode='POINTS', invert_direction=False):
        # create arc on three vertices v0 - v1 - v2
        if v0 and v1 and v2:
            world_matrix_i = world_matrix.copy()
            world_matrix_i.invert()
            # points coordinates in global space
            v0gco = world_matrix * v0.co
            v1gco = world_matrix * v1.co
            v2gco = world_matrix * v2.co
            # get circle data
            circle_center, circle_radius = cls._circle_by_3_points(
                v0=v0gco,
                v1=v1gco,
                v2=v2gco
            )
            # If vectors from circle center to v0 and v2 are close to collinear, cross of them when counting
            #   rotation matrix will be very inaccuracy. To solve we can use the middle vertex to count a normal
            #   (cross of v1-v0 and v1-v2 vectors), they are in the same plane, so it works.
            #   But if 2 source edges are close to collinear, that case founded cross will be very inaccuracy.
            #   To solve - test angle between vectors from circle center to v0 and v2 - if it is close to PI,
            #   counting axis using middle vertex
            axis = None
            if isclose((circle_center - v0gco).angle(circle_center - v2gco), pi, rel_tol=0.001):
                axis = (v1gco - v0gco).cross(v1gco - v2gco)
            # We need to set in what direction new vertices will be created. Matrices in common rotates by shortest
            #   path, but we may need or not this. We need to rotate in a direction in which middle vertex (v1) is
            #   placed. So - try to get a sign of rotation.
            invert_direction_auto = False
            invert_direction_auto = cls._invert_rotation_v21(
                v0=v0gco,
                v1=v1gco,
                v2=v2gco,
                axis=axis,
                circle_center=circle_center
            )
            # if we want to manually invert direction
            if invert_direction:
                invert_direction_auto = not invert_direction_auto
            # now we can create vertices on all point with more properly rotating direction
            # recalculate points amount if not "POINTS" mode is used
            if arc_mode == 'EDGE_LENGTH':
                # EDGES_LENGTH mode - recalculate points amount by edges_length input parameter
                v0_cc_v2_angle = (circle_center - v0gco).angle(circle_center - v2gco)
                if invert_direction_auto:
                    v0_cc_v2_angle = 2*pi - v0_cc_v2_angle
                chord_length = v0_cc_v2_angle * circle_radius
                # points = floor(chord_length / edge_length) - 1
                points = round(chord_length / edge_length) - 1
            elif arc_mode == 'BASE_32':
                # BASE_32 mode - recalculate points amount by base_points input parameter
                # get total points by base points: 32 * circle radius
                total_points_on_full_circle = int(base_points * circle_radius / base_points_rate)   # Paul + division on base_points_rate
                # get sector angle
                v0_cc_v2_angle = (circle_center - v0gco).angle(circle_center - v2gco)   # 0...2*pi
                if invert_direction_auto:
                    v0_cc_v2_angle = 2 * pi - v0_cc_v2_angle
                v0_cc_v2_ratio = v0_cc_v2_angle / (2 * pi)    # from 0...2*pi radians to 0...1 ratio
                # count points number for this sector
                points = int(total_points_on_full_circle * v0_cc_v2_ratio)
            elif arc_mode == 'BASE_32_LOG':
                # BASE_32_LOG mode - recalculate points amount by base_points and base_points_rate input parameter
                # get total points by base points: 32 * circle radius
                total_points_on_full_circle = int(base_points * log(circle_radius + 1, base_points_rate))
                # get sector angle
                v0_cc_v2_angle = (circle_center - v0gco).angle(circle_center - v2gco)   # 0...2*pi
                if invert_direction_auto:
                    v0_cc_v2_angle = 2 * pi - v0_cc_v2_angle
                v0_cc_v2_ratio = v0_cc_v2_angle / (2 * pi)    # from 0...2*pi radians to 0...1 ratio
                # count points number for this sector
                points = int(total_points_on_full_circle * v0_cc_v2_ratio)
            # create new vertices - count coordinates for all new vertices on the arc
            new_vertices_co = cls._arc_vertices_co(
                v0=v0gco,
                v2=v2gco,
                axis=axis,
                circle_center=circle_center,
                points=points,
                invert_direction=invert_direction_auto
            )
            # create new vertices by given list of coordinates
            new_vertices = []
            for vert in new_vertices_co:
                vert = world_matrix_i * vert    # get local coordinates from global
                new_vert = bm.verts.new(vert)
                new_vertices.append(new_vert)
            # create edges by new vertices
            # add v0 and v2 as first - last vertices to the new vertices list
            if invert_direction_auto:
                new_vertices.insert(0, v2)
                new_vertices.append(v0)
            else:
                new_vertices.insert(0, v0)
                new_vertices.append(v2)
            chunks = cls._chunks(new_vertices, n=2, offset=1)
            for chunk in (chunk for chunk in chunks if len(chunk) == 2):
                bm.edges.new(chunk)
            # remove v1 (clearing old geometry)
            bm.verts.remove(v1)

            # --- debug circle ---
            # bpy.context.scene.cursor_location = circle_center
            # # get normal to plane formed by all this 3 points
            # normal = (v1gco - v0gco).cross(v1gco - v2gco)
            # normal.normalize()
            # # transform matrix from normal (0, 0, 1) to normal of the 3 points plane
            # bpy.ops.mesh.primitive_circle_add(radius=circle_radius, location=circle_center)
            # m = cls._transform_matrix_from_normal_to_normal(
            #     src_normal=Vector((0.0, 0.0, 1.0)),
            #     dest_normal=normal
            # )
            # bpy.context.object.matrix_world *= m

    @staticmethod
    def _invert_rotation_v3(v0: Vector, v1: Vector, v2: Vector):
        # check if we need invert rotation direction
        #   EXPERIMENTAL, doesn't work right
        #   try to check if v1 lies left or right according to v0-v1 vector
        v1_v0 = (v0 - v1).normalized()
        v1_v2 = (v2 - v1).normalized()
        normal = v1_v0.cross(v1_v2)
        normal_v1 = v1_v0.cross(normal)
        dot = normal_v1.dot(v1_v0)
        return dot < 0.0

    @classmethod
    def _invert_rotation_v21(cls,  v0, v1, v2, axis, circle_center, points=25):
        # check if we need invert rotation direction
        #   More heavy version of Paul variant
        #   counting coordinates for 25 points in two directions, finding the closest point to v1
        #   in which way the closest point is found - use this variant
        single_vert_co = cls._arc_vertices_co(
            v0=v0,
            v2=v2,
            axis=axis,
            circle_center=circle_center,
            points=points,
            invert_direction=False
        )
        single_vert_co_inv_direction = cls._arc_vertices_co(
            v0=v0,
            v2=v2,
            axis=axis,
            circle_center=circle_center,
            points=points,
            invert_direction=True
        )
        # find the closes point for v1 in both lists
        v1_point_min_dist = min(map(lambda _v: (_v - v1).length, single_vert_co))
        v1_point_min_dist_inv = min(map(lambda _v: (_v - v1).length, single_vert_co_inv_direction))
        # if min dist for no inversion is less min dist for inverted variant - return no inversion
        #   False if v1_point_min_dist < v1_point_min_dist_inv else True
        return v1_point_min_dist > v1_point_min_dist_inv

    @classmethod
    def _invert_rotation_v2(cls,  v0, v1, v2, axis, circle_center):
        # check if we need invert rotation direction
        #   Paul variant
        #   Best variant but doesn't work in all cases
        #   built two arcs by 1 point, one arc with default rotation direction and another with opposite
        #   try to compare length of two vectors in these variants: v1 - center of created arc
        #   which variant length is shorter - this direction we will use
        single_vert_co = cls._arc_vertices_co(
            v0=v0,
            v2=v2,
            axis=axis,
            circle_center=circle_center,
            points=1,
            invert_direction=True
        )
        single_vert_co_inv_direction = cls._arc_vertices_co(
            v0=v0,
            v2=v2,
            axis=axis,
            circle_center=circle_center,
            points=1,
            invert_direction=False
        )
        # compare lengths
        l1 = (v1 - single_vert_co[0]).length
        l2 = (v1 - single_vert_co_inv_direction[0]).length
        return l1 < l2

    @staticmethod
    def _invert_rotation_v1(v0: Vector, v1: Vector, circle_center: Vector):
        # check if we need invert rotation direction
        #   EXPERIMENTAL, doesn't work right in all cases
        #   check angle between two vectors: v0-circle_center and v1-circle_center
        #   if angle is obtuse - invert direction
        #   If cc_v0 and cc_v1 vectors are close to collinear, angle_sign may be very inaccurate.
        #       if isclose(angle_sign, 0, abs_tol=0.001):   # collinear
        cc_v0g = (v0 - circle_center).normalized()
        cc_v1g = (v1 - circle_center).normalized()
        angle_sign = cc_v0g.dot(cc_v1g)
        # change rotating direction by sign
        return True if angle_sign < 0.0 else False

    @staticmethod
    def _is_vectors_close(v0: Vector, v1: Vector, tolerance=0.001):
        # check if two vectors are close by tolerance
        return True if isclose(v0.x, v1.x, rel_tol=tolerance) \
                       and isclose(v0.y, v1.y, rel_tol=tolerance) \
                       and isclose(v0.z, v1.z, rel_tol=tolerance) \
            else False

    @classmethod
    def _arc_vertices_co(cls, v0, v2, axis, circle_center, points, invert_direction):
        # get coordinates of new creating vertices on the arc in GLOBAL space
        # return [(x, y, z), (x, y, z), ...]
        new_vertices_co = []
        for i in range(points):
            # count rotating factor for current creating vertex
            f = (i + 1) / (points + 1)
            # create new vertex in v0 in global system
            new_vert_co = v0
            # we need to use angle v0 - v1 - v2 because we need to control direction of rotating
            #   to make rotation through v1 but not by shortest path (as commonly matrices works)
            v0_cc_v2_angle = (circle_center - v0).angle(circle_center - v2)
            angle = (v0_cc_v2_angle * f) if not invert_direction else (2 * pi - v0_cc_v2_angle) * f + v0_cc_v2_angle
            # create rotation matrix for new vertex to rotate it from v0 to its place by factor
            mat_rot = cls._rotation_matrix_from_vector_to_vector(
                src_vector=circle_center - v0,
                dest_vector=circle_center - v2,
                angle=angle,
                axis=axis
            )
            # rotate new vertex using counted matrices
            #   move it to the world origin (because matrices rotates around world origin), rotate, move back
            cc_v0_g = v0 - circle_center
            cc_v0_g = mat_rot * cc_v0_g
            cc_new_vert_g = circle_center + cc_v0_g
            # save it in list
            new_vertices_co.append(cc_new_vert_g)
        return new_vertices_co

    @staticmethod
    def _chunks(lst, n, offset=0):
        for i in range(0, len(lst), n - offset):
            yield lst[i:i + n]

    @staticmethod
    def _invert_vector_around_point(vector: Vector, point: Vector):
        # get coordinates of inverted vector around point
        vector_from_point = vector - point
        vector_from_point_inverted = vector_from_point * -1
        return point + vector_from_point_inverted

    @staticmethod
    def _circle_by_3_points(v0: Vector, v1: Vector, v2: Vector):
        # get circle from 3 points (coordinates in global space)
        #   return circle center coordinates and radius length
        vv1 = v1 - v0
        vv2 = v2 - v0
        def _dot(_v1, _v2):
            return _v1.x * _v2.x + _v1.y * _v2.y + _v1.z * _v2.z
        v1v1 = _dot(vv1, vv1)
        v2v2 = _dot(vv2, vv2)
        v1v2 = _dot(vv1, vv2)
        base = 0.5 / (v1v1 * v2v2 - v1v2 * v1v2)
        k1 = base * v2v2 * (v1v1 - v1v2)
        k2 = base * v1v1 * (v2v2 - v1v2)
        center = v0 + vv1 * k1 + vv2 * k2
        radius = (center - v0).length
        return center, radius

    @staticmethod
    def _rotation_matrix_from_vector_to_vector(src_vector: Vector, dest_vector: Vector, angle: float = None,
                                               axis: Vector = None):
        # get transform matrix to rotate one vector to another around axis
        # count full angle between two vectors, or we can count for custom angle if it is set in parameters
        angle = angle if angle else src_vector.angle(dest_vector)  # rad
        # manual setting axis may need if src and dest vectors are close to collinear (cross counting with big inaccuracy)
        if axis is None:
            # normal to plane formed with this two vectors
            axis = src_vector.cross(dest_vector)
            axis.normalize()
        # count rotation matrix
        return Matrix.Rotation(angle, 4, axis)

    @staticmethod
    def _transform_matrix_from_normal_to_normal(src_normal: Vector, dest_normal: Vector) -> Matrix:
        # get transform matrix from one normal to another
        v = dest_normal.cross(src_normal)
        c = dest_normal.dot(src_normal)
        if isclose(c, -1.0, rel_tol=0.0001):
            # vectors are collinear and opposite-directional - matrix just for inverse
            transform_matrix = Matrix(
                ((-1.0, 0.0, 0.0),
                 (0.0, -1.0, 0.0),
                 (0.0, 0.0, -1.0))
            )
        else:
            k = 1.0 / (1.0 + c)
            transform_matrix = Matrix(
                ((v.x * v.x * k + c, v.y * v.x * k - v.z, v.z * v.x * k + v.y),
                 (v.x * v.y * k + v.z, v.y * v.y * k + c, v.z * v.y * k - v.x),
                 (v.x * v.z * k - v.y, v.y * v.z * k + v.x, v.z * v.z * k + c))
            )
        # convert to 4x4 matrix for Blender transformations
        transform_matrix = transform_matrix.to_4x4()
        # convert from left-hand orientation to Blender right-hand orientation
        transform_matrix.transpose()
        # return final transformation matrix
        return transform_matrix

    @staticmethod
    def _possible_path(start_point, possible_points, possible_endpoints):
        # try to find a path starting from "start_point" which can follow only through "possible_points" (vertices and
        #   edges must be selected) and can end on any point from "possible_endpoints"
        path_found = None
        visited_vertices = {start_point, }
        possible_paths = [[start_point, ], ]
        _i = 0
        _l = len(possible_points)
        while not path_found and possible_paths:
            # explore first possible path
            path = possible_paths.pop(0)
            # get possible next vertices
            vertex = path[-1]
            next_vertices = [_edge.other_vert(vertex) for _edge in vertex.link_edges
                             if _edge.select    # edges must be selected to prevent case when two selected vertices could be connected with not selected edge
                             and _edge.other_vert(vertex) in possible_points
                             and _edge.other_vert(vertex) not in visited_vertices]
            # create new possible paths
            for vertex in next_vertices:
                possible_paths.insert(0, path[:] + [vertex,])   # new paths are prepending to the list
            # update visited vertices list
            visited_vertices.update(next_vertices)
            # check if successfully finished with one of possible_endpoints
            path_found = next((_path for _path in possible_paths if _path[-1] in possible_endpoints), None)
            # overflow control
            _i += 1
            if _i > _l:
                print('ERR: overflow when finding possible paths')
                break
        return path_found

    @staticmethod
    def ui(layout, context):
        # ui panel
        # 3 Points Arc
        box = layout.box()
        op = box.operator(
            operator='arc3.arc3',
            icon='PARTICLE_POINT'
        )
        op.points = context.scene.arc3_prop_points
        op.edge_length = context.scene.arc3_prop_edge_length
        op.base_points = context.scene.arc3_prop_base_points
        op.base_points_rate = context.scene.arc3_prop_base_points_rate
        op.mode = context.scene.arc3_prop_mode
        op.invert_direction = False
        # Butch Clean
        op = box.operator(
            operator='arc3.butch_clean',
            icon='FORCE_TURBULENCE',
            text='Arc Batch Clean'
        )
        op.base_points = context.scene.arc3_prop_base_points
        op.base_points_rate = context.scene.arc3_prop_base_points_rate
        op.mode = context.scene.arc3_prop_mode
        # PROPS
        if context.scene.arc3_prop_mode == 'EDGE_LENGTH':
            box.prop(
                data=context.scene,
                property='arc3_prop_edge_length'
            )
        elif context.scene.arc3_prop_mode == 'BASE_32':
            box.prop(
                data=context.scene,
                property='arc3_prop_base_points'
            )
            box.prop(
                data=context.scene,
                property='arc3_prop_base_points_rate',
                text='K'
            )
        elif context.scene.arc3_prop_mode == 'BASE_32_LOG':
            box.prop(
                data=context.scene,
                property='arc3_prop_base_points'
            )
            box.prop(
                data=context.scene,
                property='arc3_prop_base_points_rate'
            )
        else:
            box.prop(
                data=context.scene,
                property='arc3_prop_points'
            )
        row = box.row()
        row.prop(
            data=context.scene,
            property='arc3_prop_mode',
            expand=True
        )


# OPERATORS

class Arc3_OT_arc3(Operator):
    bl_idname = 'arc3.arc3'
    bl_label = '3 Points Arc'
    bl_options = {'REGISTER', 'UNDO'}

    points = IntProperty(
        name='Points Amount',
        default=5
    )
    edge_length = FloatProperty(
        name='Edge Length',
        default=1.25,
        min=0.0001
    )
    base_points = IntProperty(
        name='Points Base',
        default=32,
        min=1
    )
    base_points_rate = FloatProperty(
        name='K/LOG',
        default=2.0,
        min=0.0001
    )
    mode = EnumProperty(
        name='Mode',
        items=[
            ('POINTS', 'POINTS', 'POINTS AMOUNT', '', 0),
            ('EDGE_LENGTH', 'LENGTH', 'EDGE LENGTH', '', 1),
            ('BASE_32', 'LIN', 'LOG', '', 2),
            ('BASE_32_LOG', 'LOG', 'LOG', '', 3)
        ],
        default='POINTS'
    )
    invert_direction = BoolProperty(
        name='Invert Direction',
        default=False
    )

    def execute(self, context):
        Arc3.arc3(
            context=context,
            obj=context.active_object,
            points=self.points,
            edge_length=self.edge_length,
            base_points=self.base_points,
            base_points_rate=self.base_points_rate,
            arc_mode=self.mode,
            invert_direction=self.invert_direction
        )
        return {'FINISHED'}

class Arc3_OT_butch_clean(Operator):
    bl_idname = 'arc3.butch_clean'
    bl_label = 'Arc Batch Clean'
    bl_options = {'REGISTER', 'UNDO'}

    points = IntProperty(
        name='Points Amount',
        default=5
    )
    edge_length = FloatProperty(
        name='Edge Length',
        default=1.25,
        min=0.0001
    )
    base_points = IntProperty(
        name='Points Base',
        default=32,
        min=1
    )
    base_points_rate = FloatProperty(
        name='K/LOG',
        default=2.0,
        min=0.0001
    )
    mode = EnumProperty(
        name='Mode',
        items=[
            ('POINTS', 'POINTS', 'POINTS AMOUNT', '', 0),
            ('EDGE_LENGTH', 'LENGTH', 'EDGE LENGTH', '', 1),
            ('BASE_32', 'LIN', 'LIN', '', 2),
            ('BASE_32_LOG', 'LOG', 'LOG', '', 3)
        ],
        default='POINTS'
    )

    def execute(self, context):
        Arc3.butch_clean(
            context=context,
            obj=context.active_object,
            points=self.points,
            edge_length=self.edge_length,
            base_points=self.base_points,
            base_points_rate=self.base_points_rate,
            arc_mode=self.mode
        )
        return {'FINISHED'}


# PANELS

class Arc3_PT_panel(Panel):
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_label = '3 Points Arc'
    bl_category = '1D'

    def draw(self, context):
        Arc3.ui(
            layout=self.layout,
            context=context
        )


# REGISTER

def register(ui=True):
    # 3-arc
    Scene.arc3_prop_mode = EnumProperty(
        name='Mode',
        items=[
            ('POINTS', 'POINTS', 'POINTS AMOUNT', '', 0),
            ('EDGE_LENGTH', 'LENGTH', 'EDGE LENGTH', '', 1),
            ('BASE_32', 'LIN', 'LIN', '', 2),
            ('BASE_32_LOG', 'LOG', 'LOG', '', 3)
        ],
        default='POINTS'
    )
    Scene.arc3_prop_points = IntProperty(
        name='Points Amount',
        default=5,
        min=0
    )
    Scene.arc3_prop_edge_length = FloatProperty(
        name='Edge Length',
        default=1.25,
        min=0.0001
    )
    # Scene.arc3_prop_invert_direction = BoolProperty(
    #     name='Invert Direction',
    #     default=False
    # )
    register_class(Arc3_OT_arc3)
    # butch clean
    Scene.arc3_prop_base_points = IntProperty(
        name='Points Base',
        default=32,
        min=1
    )
    Scene.arc3_prop_base_points_rate = FloatProperty(
        name='Log',
        default=2.0,
        min=0.0001
    )
    register_class(Arc3_OT_butch_clean)
    if ui:
        register_class(Arc3_PT_panel)


def unregister(ui=True):
    if ui:
        unregister_class(Arc3_PT_panel)
    # butch clean
    unregister_class(Arc3_OT_butch_clean)
    del Scene.arc3_prop_base_points_rate
    del Scene.arc3_prop_base_points
    # 3-arc
    unregister_class(Arc3_OT_arc3)
    # del Scene.arc3_prop_invert_direction
    del Scene.arc3_prop_edge_length
    del Scene.arc3_prop_points
    del Scene.arc3_prop_mode


if __name__ == "__main__":
    register()
