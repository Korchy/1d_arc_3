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
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        bm.verts.ensure_lookup_table()
        bm.edges.ensure_lookup_table()

        obj_world_matrix = obj.matrix_world.copy()

        # get 3 starting vertices
        v0 = v1 = v2 = None
        v0 = next((vertex for vertex in bm.verts if vertex.select), None)
        if v0:
            v1 = v0.link_edges[0].other_vert(v0) if v0.link_edges else None
            if v1:
                v2_edge = next((edge for edge in v1.link_edges if edge != v0.link_edges[0]), None)
                v2 = v2_edge.other_vert(v1) if v2_edge else None

        # points_normal = (v0.co - v1.co).cross(v2.co - v1.co)
        points_normal = (v1.co * obj_world_matrix - v0.co * obj_world_matrix).cross(v1.co * obj_world_matrix - v2.co * obj_world_matrix)
        points_normal.normalize()

        if v0 and v1 and v2:

            # https://github.com/sergarrido/random/blob/master/circle3d/circle3d.cpp#L63
            # p1 = v0.co
            # p2 = v1.co
            # p3 = v2.co
            p1 = v0.co * obj_world_matrix
            p2 = v1.co * obj_world_matrix
            p3 = v2.co * obj_world_matrix

            v1 = p2 - p1
            v2 = p3 - p1

            v1v1 = v1.dot(v1)
            v2v2 = v2.dot(v2)
            v1v2 = v1.dot(v2)

            base = 0.5 / (v1v1 * v2v2 - v1v2 * v1v2)
            k1 = base * v2v2 * (v1v1 - v1v2)
            k2 = base * v1v1 * (v2v2 - v1v2)

            c = p1 + v1 * k1 + v2 * k2
            print('c', c)
            bpy.context.scene.cursor_location = c
            radius = (c - p1).length
            print('radius', radius)

            # for i in range(points):
            #     bm.verts.new(v)

            bpy.ops.mesh.primitive_circle_add(radius=radius, location=c)
            # bpy.ops.mesh.primitive_circle_add(radius=radius, location=(0, 0, 0))

            # transform matrix from normal (0, 0, 1) to normal of the 3 points plane

            m = cls._transform_from_normal_to_normal(
                src_normal=Vector((0.0, 0.0, 1.0)),
                dest_normal=points_normal
            )
            print(m)

            print(bpy.context.object)

            bpy.context.object.matrix_world *= m

        # save changed data to mesh
        bm.to_mesh(obj.data)
        bm.free()
        # return mode back
        # bpy.ops.object.mode_set(mode=mode)

    @staticmethod
    def _transform_from_normal_to_normal(src_normal, dest_normal):

        z = dest_normal
        d = src_normal

        # v = cross(z, d);
        v = z.cross(d)
        # c = dot(z, d);
        c = z.dot(d)
        k = 1.0 / (1.0 + c)

        m = Matrix(
            ((v.x * v.x * k + c, v.y * v.x * k - v.z, v.z * v.x * k + v.y),
             (v.x * v.y * k + v.z, v.y * v.y * k + c, v.z * v.y * k - v.x),
             (v.x * v.z * k - v.y, v.y * v.z * k + v.x, v.z * v.z * k + c))
        )
        m = m.to_4x4()

        m.transpose()
        
        return m

        # https://stackoverflow.com/questions/52189123/calculate-rotation-matrix-to-transform-one-vector-to-another

        # # V = normalize(cross(src_normal, dest_normal));
        # V = src_normal.cross(dest_normal)
        # V.normalize()
        #
        # # phi = acos(dot(src_normal, dest_normal));
        # phi = acos(src_normal.dot(dest_normal))
        #
        # # rcos = cos(phi);
        # rcos = cos(phi)
        #
        # # rsin = sin(phi);
        # rsin = sin(phi)
        #
        # M = Matrix()
        #
        # M[0][0] = rcos + V.x * V.x * (1.0 - rcos)
        # M[1][0] = V.z * rsin + V.y * V.x * (1.0 - rcos)
        # M[2][0] = -V.y * rsin + V.z * V.x * (1.0 - rcos)
        # M[3][0] = 0.0
        #
        # M[0][1] = -V.z * rsin + V.x * V.y * (1.0 - rcos)
        # M[1][1] = rcos + V.y * V.y * (1.0 - rcos)
        # M[2][1] = -V.x * rsin + V.z * V.y * (1.0 - rcos)
        # M[3][1] = 0.0
        #
        # M[0][2] = V.y * rsin + V.x * V.z * (1.0 - rcos)
        # M[1][2] = -V.x * rsin + V.y * V.z * (1.0 - rcos)
        # M[2][2] = rcos + V.z * V.z * (1.0 - rcos)
        # M[3][2] = 0.0
        #
        # M[0][3] = 0.0
        # M[1][3] = 0.0
        # M[2][3] = 0.0
        # M[3][3] = 1.0
        #
        # return M

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
