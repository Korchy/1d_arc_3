# Nikita Akimov
# interplanety@interplanety.org
#
# GitHub
#    https://github.com/Korchy/1d_arc_3

import bmesh
import bpy
from bpy.props import IntProperty
from bpy.types import Operator, Panel, Scene
from bpy.utils import register_class, unregister_class
from mathutils import Matrix, Vector
from math import cos, sin, acos

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
    def arc3(cls, context, obj, points=5):
        # create arc from 2 edges (3 vertices)
        obj = obj if obj else context.active_object
        # current mode
        mode = obj.mode
        if obj.mode == 'EDIT':
            bpy.ops.object.mode_set(mode='OBJECT')
        # switch to vertex selection mode
        # context.tool_settings.mesh_select_mode = (True, False, False)
        obj_world_matrix = obj.matrix_world.copy()
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

            # get normal to plane formed by all this 3 points
            normal = (v1gco - v0gco).cross(v1gco - v2gco)
            normal.normalize()

            bpy.context.scene.cursor_location = circle_center

            # for i in range(points):
            #     bm.verts.new(v)

            bpy.ops.mesh.primitive_circle_add(radius=circle_radius, location=circle_center)

            # transform matrix from normal (0, 0, 1) to normal of the 3 points plane
            m = cls._transform_matrix_from_normal_to_normal(
                src_normal=Vector((0.0, 0.0, 1.0)),
                dest_normal=normal
            )
            print('m', m)
            print('q', m.to_quaternion())

            bpy.context.object.matrix_world *= m

        # save changed data to mesh
        bm.to_mesh(obj.data)
        bm.free()
        # return mode back
        # bpy.ops.object.mode_set(mode=mode)

    @staticmethod
    def _circle_by_3_points(v0: Vector, v1: Vector, v2: Vector):
        # get circle from 3 points (coordinates in global space)
        #   return circle center coordinates and radius length
        vv1 = v1 - v0
        vv2 = v2 - v0
        v1v1 = vv1.dot(vv1)
        v2v2 = vv2.dot(vv2)
        v1v2 = vv1.dot(vv2)
        base = 0.5 / (v1v1 * v2v2 - v1v2 * v1v2)
        k1 = base * v2v2 * (v1v1 - v1v2)
        k2 = base * v1v1 * (v2v2 - v1v2)
        center = v0 + vv1 * k1 + vv2 * k2
        radius = (center - v0).length
        return center, radius

    @staticmethod
    def _transform_matrix_from_normal_to_normal(src_normal: Vector, dest_normal: Vector):
        # get transform matrix from one normal to another
        v = dest_normal.cross(src_normal)
        c = dest_normal.dot(src_normal)
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
        layout.prop(
            data=context.scene,
            property='arc3_prop_points'
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

    def execute(self, context):
        Arc3.arc3(
            context=context,
            obj=context.active_object,
            points=self.points
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
    register_class(Arc3_OT_arc3)
    if ui:
        register_class(Arc3_PT_panel)


def unregister(ui=True):
    if ui:
        unregister_class(Arc3_PT_panel)
    unregister_class(Arc3_OT_arc3)
    del Scene.arc3_prop_points


if __name__ == "__main__":
    register()
