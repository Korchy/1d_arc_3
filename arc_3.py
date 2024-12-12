# Nikita Akimov
# interplanety@interplanety.org
#
# GitHub
#    https://github.com/Korchy/1d_arc_3

import bpy
from bpy.props import IntProperty
from bpy.types import Operator, Panel, Scene
from bpy.utils import register_class, unregister_class

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

    @staticmethod
    def arc3(context, obj, points=5):
        # create arc from 2 edges (3 vertices)
        print('arc3')

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
