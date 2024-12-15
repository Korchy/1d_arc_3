# Nikita Akimov
# interplanety@interplanety.org
#
# GitHub
#    https://github.com/Korchy/1d_arc_3

import bmesh
import bpy
from bpy.props import BoolProperty, IntProperty
from bpy.types import Operator, Panel, Scene
from bpy.utils import register_class, unregister_class
from mathutils import Matrix, Vector, Quaternion
from math import cos, sin, acos, isclose, pi

bl_info = {
    "name": "3 Points ARc",
    "description": "Creates arc from 2 edges (3 starting points).",
    "author": "Nikita Akimov, Paul Kotelevets",
    "version": (1, 0, 0),
    "blender": (2, 79, 0),
    "location": "View3D > Tool panel > 1D > 3 Points Arc",
    "doc_url": "https://github.com/Korchy/1d_arc_3",
    "tracker_url": "https://github.com/Korchy/1d_arc_3",
    "category": "All"
}


# MAIN CLASS

class Arc3:

    @classmethod
    def arc3(cls, context, obj, points=5, invert_direction=False):
        # create arc from 2 edges (3 vertices)
        obj = obj if obj else context.active_object
        # current mode
        mode = obj.mode
        if obj.mode == 'EDIT':
            bpy.ops.object.mode_set(mode='OBJECT')
        # switch to vertex selection mode
        # context.tool_settings.mesh_select_mode = (True, False, False)
        obj_world_matrix = obj.matrix_world.copy()
        obj_world_matrix_i = obj.matrix_world.copy()
        obj_world_matrix_i.invert()
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        bm.verts.ensure_lookup_table()
        bm.edges.ensure_lookup_table()
        # get 3 starting vertices
        v0 = v1 = v2 = None
        v0 = next((vertex for vertex in bm.verts if vertex.select), None)
        if v0:
            v1 = v0.link_edges[0].other_vert(v0) if v0.link_edges else None
            if v1:
                v2_edge = next((edge for edge in v1.link_edges if edge != v0.link_edges[0]), None)
                v2 = v2_edge.other_vert(v1) if v2_edge else None
        # if we have all three starting points
        if v0 and v1 and v2:
            # points coordinates in global space
            v0gco = obj_world_matrix * v0.co
            v1gco = obj_world_matrix * v1.co
            v2gco = obj_world_matrix * v2.co
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
            cc_v0g = (v0gco - circle_center).normalized()
            cc_v1g = (v1gco - circle_center).normalized()
            angle_sign = cc_v0g.dot(cc_v1g)
            # If cc_v0 and cc_v1 vectors are close to collinear, angle_sign may be very inaccurate. To fix - use
            #   the input parameter invert_direction
            # change rotating direction by sign
            if angle_sign < 0.0:
                invert_direction = not invert_direction
            # create new vertices - create in v0 and rotate it around circle center
            new_vertices = []
            for i in range(points):
                # count rotating factor for current creating vertex
                f = (i + 1) / (points + 1)
                # create new vertex in v0 in global system
                new_vert = bm.verts.new(v0gco)
                # save it for future creating edges
                new_vertices.append(new_vert)
                # we need to use angle v0 - v1 - v2 because we need to control direction of rotating
                #   to make rotation through v1 but not by shortest path (as commonly matrices works)
                v0_cc_v2_angle = (circle_center - v0gco).angle(circle_center - v2gco)
                angle = (v0_cc_v2_angle * f) if not invert_direction else (2*pi - v0_cc_v2_angle) * f + v0_cc_v2_angle
                # create rotation matrix for new vertex to rotate it from v0 to its place by factor
                mat_rot = cls._rotation_matrix_from_vector_to_vector(
                    src_vector=circle_center - v0gco,
                    dest_vector=circle_center - v2gco,
                    angle=angle,
                    axis=axis
                )
                # rotate new vertex using counted matrices
                #   move it to the world origin (because matrices rotates around world origin), rotate, move back
                cc_v0_g = v0gco - circle_center
                cc_v0_g = mat_rot * cc_v0_g
                cc_new_vert_g = circle_center + cc_v0_g
                cc_new_vert = obj_world_matrix_i * cc_new_vert_g
                new_vert.co = cc_new_vert
            # create edges by new vertices
            # add v0 and v2 as first - last vertices to the new vertices list
            if invert_direction:
                new_vertices.insert(0, v2)
                new_vertices.append(v0)
            else:
                new_vertices.insert(0, v0)
                new_vertices.append(v2)
            chunks = cls._chunks(new_vertices, n=2, offset=1)
            for chunk in (chunk for chunk in chunks if len(chunk) == 2):
                bm.edges.new(chunk)

            # --- debug circle ---
            # bpy.context.scene.cursor_location = circle_center
            # get normal to plane formed by all this 3 points
            # normal = (v1gco - v0gco).cross(v1gco - v2gco)
            # normal.normalize()
            # transform matrix from normal (0, 0, 1) to normal of the 3 points plane
            # bpy.ops.mesh.primitive_circle_add(radius=circle_radius, location=circle_center)
            # m = cls._transform_matrix_from_normal_to_normal(
            #     src_normal=Vector((0.0, 0.0, 1.0)),
            #     dest_normal=normal
            # )
            # bpy.context.object.matrix_world *= m

        # save changed data to mesh
        bm.to_mesh(obj.data)
        bm.free()
        # return mode back
        bpy.ops.object.mode_set(mode=mode)

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
        v1v1 = vv1.dot(vv1)
        v2v2 = vv2.dot(vv2)
        v1v2 = vv1.dot(vv2)
        base = (0.5 / (v1v1 * v2v2 - v1v2 * v1v2)) if ((v1v1 * v2v2 - v1v2 * v1v2) > 0.001) else 0
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
    def ui(layout, context):
        # ui panel
        # 3 Points Arc
        op = layout.operator(
            operator='arc3.arc3',
            icon='PARTICLE_POINT'
        )
        op.points = context.scene.arc3_prop_points
        # op.invert_direction = context.scene.arc3_prop_invert_direction
        op.invert_direction = False
        layout.prop(
            data=context.scene,
            property='arc3_prop_points'
        )
        layout.prop(
            data=context.scene,
            property='arc3_prop_invert_direction'
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
    invert_direction = BoolProperty(
        name='Invert Direction',
        default=False
    )

    def execute(self, context):
        Arc3.arc3(
            context=context,
            obj=context.active_object,
            points=self.points,
            invert_direction=self.invert_direction
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
    Scene.arc3_prop_points = IntProperty(
        name='Points Amount',
        default=5
    )
    # Scene.arc3_prop_invert_direction = BoolProperty(
    #     name='Invert Direction',
    #     default=False
    # )
    register_class(Arc3_OT_arc3)
    if ui:
        register_class(Arc3_PT_panel)


def unregister(ui=True):
    if ui:
        unregister_class(Arc3_PT_panel)
    unregister_class(Arc3_OT_arc3)
    # del Scene.arc3_prop_invert_direction
    del Scene.arc3_prop_points


if __name__ == "__main__":
    register()
