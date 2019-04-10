from io import BytesIO
import win32clipboard
from PIL import Image
import sys
from direct.showbase.ShowBase import ShowBase
from direct.gui.DirectButton import DirectButton, DirectFrame
from direct.gui.DirectLabel import DirectLabel
from panda3d.core import VBase4, AmbientLight, GeomVertexFormat
from panda3d.core import GeomVertexData, GeomVertexWriter, Geom, GeomLines
from panda3d.core import GeomNode, GeomTriangles, TextNode, AntialiasAttrib
from panda3d.core import PointLight, CollisionNode, CollisionRay
from panda3d.core import CollisionTraverser, CollisionHandlerQueue
from panda3d.core import CollisionTube, CollisionSphere
from edsFEM.renderutils.eventhandlers import ClickHandler
from edsFEM.renderutils.rendernodes import RenderNode
import numpy as np


class Renderer(ShowBase):
    """
    The Renderer class allows the visual representation of the calculations
    performed in the System class.
    The Renderer must be set up in the script calling the Renderer, including
    the buttons to be rendered on-screen and the initial state of the render.
    """

    def __init__(self, curr_file=None, scale='mm'):
        """
        Initializing the Renderer class sets up the render and prepares the
        collision detection for use.

        The renderer is designed for use on the Windows platform and is
        untested in different environments.

        The renderer uses the Panda3d open source framework for rendering.
        https://www.panda3d.org/
        """

        ShowBase.__init__(self)

        self.s = self.scale_check(scale)

        # Set the background color of the render
        self.setBackgroundColor(0.985, 0.985, 0.985, 0.4)

        # Enable fast exit
        self.accept("escape", sys.exit)

        self.window = self.winList[0]

        # Create variables
        self.curr_render = None
        self.tooltipText = None

        # Lists containing the paths to rendered nodes and beams RenderNodes
        self.Nodes = []
        self.Beams = []

        # Lists containing the paths to the various rendered objects' nodepaths
        self.sup_NodePaths = []
        self.releaseNodePaths = []
        self.displacedbeamNodePaths = []
        self.beammomentNodePaths = []
        self.beamnforcesNodePaths = []
        self.normalstressNodePaths = []
        self.shearstressNodePaths = []
        self.textNodePaths = []
        self.beamprofNodePaths = []
        self.q_loadPaths = []
        self.p_loadPaths = []
        self.m_loadPaths = []
        self.framesPaths = []
        self.plights = []
        self.lightNodePaths = []

        # Add an ambient light
        alight = AmbientLight('alight')
        alight.setColor(VBase4(0.6, 0.6, 0.6, 1))
        alnp = self.render.attachNewNode(alight)
        self.render.setLight(alnp)
        self.render.setShaderAuto()
        self.render.setAntialias(AntialiasAttrib.MAuto)

        # add point light
        plight = PointLight('plight')
        plight.setColor(VBase4(0.7, 0.7, 0.7, 1))
        plight.setAttenuation((1, 0, 0))
        plnp = self.render.attachNewNode(plight)
        plnp.setPos(0, 0, 1000)
        self.render.setLight(plnp)

        # Add a button to make screenshots
        self.make_button('screenshot')
        self.accept('control-x', self.cmd_screenshot)

        # Create the collision handling to make object clickable
        self.traverser = CollisionTraverser('traverser')
        self.cTrav = self.traverser

        self.pickerNode = CollisionNode('mouseRay')
        self.pickerNP = self.camera.attachNewNode(self.pickerNode)
        self.pickerNode.setFromCollideMask(GeomNode.getDefaultCollideMask())
        self.pickerRay = CollisionRay()
        self.pickerNode.addSolid(self.pickerRay)

        self.myHandler = CollisionHandlerQueue()

        self.click_handler = ClickHandler(self, self.s)

        self.traverser.addCollider(self.pickerNP, self.myHandler)

    def render_tooltip(self, text):
        """
        This function renders a tooltip containing the given text at the
        current mouse position.
        """

        if self.mouseWatcherNode.hasMouse():

            # get the mouse position
            x = self.mouseWatcherNode.getMouseX()
            y = self.mouseWatcherNode.getMouseY()

            wp = self.win.getProperties()
            aspX = 1.0
            aspY = 1.0
            wpXSize = wp.getXSize()
            wpYSize = wp.getYSize()
            if wpXSize > wpYSize:
                aspX = wpXSize / float(wpYSize)
            else:
                aspY = wpYSize / float(wpXSize)

            self.tooltipText = DirectLabel(text=text,
                                           text_fg=(1, 1, 1, 1),
                                           text_scale=0.05,
                                           text_align=TextNode.ALeft,
                                           frameColor=(0.15, 0.15, 0.15, 0.85),
                                           pad=(0.02, 0.02),
                                           pos=((x*aspX), 0, (y*aspY)))
            self.tooltipText.setTransparency(True)

            self.tooltipText.show()

    def hide_tooltip(self):
        """
        This function hides any rendered tooltip contained in the
        self.tooltipText variable
        """

        if self.tooltipText is not None:
            self.tooltipText.hide()

    def render_origin(self, scale='mm'):
        """
        This function renders the origin by way of rendering the three
        axes of the Cartesian axis system.
        """
        s = self.scale_check(scale)

        fmt = GeomVertexFormat.getV3n3c4()
        vertexData = GeomVertexData('origin', fmt, Geom.UHStatic)

        vertexData.setNumRows(4)
        vertices = GeomVertexWriter(vertexData, 'vertex')
        normals = GeomVertexWriter(vertexData, 'normal')
        colors = GeomVertexWriter(vertexData, 'color')

        vertices.addData3f(0, 0, 0)
        normals.addData3f(0, 0, 1)
        colors.addData4f(0, 255, 0, 1)
        vertices.addData3f(400/s, 0, 0)
        normals.addData3f(0, 0, 1)
        colors.addData4f(255, 0, 0, 1)
        vertices.addData3f(0, 400/s, 0)
        normals.addData3f(0, 0, 1)
        colors.addData4f(0, 255, 0, 1)
        vertices.addData3f(0, 0, 400/s)
        normals.addData3f(0, 0, 1)
        colors.addData4f(0, 0, 255, 1)

        primitive = GeomLines(Geom.UHStatic)

        primitive.addVertices(0, 1)
        primitive.addVertices(0, 2)
        primitive.addVertices(0, 3)

        geom = Geom(vertexData)
        geom.addPrimitive(primitive)

        node = GeomNode('Origin')
        node.addGeom(geom)

        self.originPath = self.render.attachNewNode(node)
        self.originPath.setRenderModeThickness(2)
        self.originPath.setRenderModePerspective(True)

    def render_p_loads(self, node, scale='mm', dscale=250):
        """
        This function renders all point loads on the structural system as
        arrows pointing in the load direction at the scale given in the
        dscale argument.
        """
        s = self.scale_check(scale)

        fmt = GeomVertexFormat.getV3n3c4()
        vertexData = GeomVertexData('p_load', fmt, Geom.UHStatic)

        p = 5
        n = 0
        Fx = 0
        Fy = 0
        Fz = 0
        if node.Fx != 0:
            n += 1
            Fx = node.Fx
        if node.Fy != 0:
            n += 1
            Fy = node.Fy
        if node.Fz != 0:
            n += 1
            Fz = node.Fz

        x = node.point.x
        y = node.point.y
        z = node.point.z

        if n > 0:

            n *= 6
            vertexData.setNumRows(n)

            vertices = GeomVertexWriter(vertexData, 'vertex')
            normals = GeomVertexWriter(vertexData, 'normal')
            colors = GeomVertexWriter(vertexData, 'color')

            prim_lines = GeomLines(Geom.UHStatic)
            prim_triangles = GeomTriangles(Geom.UHStatic)

            g = 0

            if Fx != 0:

                vertices.addData3f((x)/s, (y)/s, (z)/s)
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - Fx / dscale), (y / s), (z / s))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Fx / dscale) / p)),
                                   (y / s - ((Fx / dscale) / (p*2))),
                                   (z / s - ((Fx / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Fx / dscale) / p)),
                                   (y / s - ((Fx / dscale) / (p*2))),
                                   (z / s + ((Fx / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Fx / dscale) / p)),
                                   (y / s + ((Fx / dscale) / (p*2))),
                                   (z / s + ((Fx / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Fx / dscale) / p)),
                                   (y / s + ((Fx / dscale) / (p*2))),
                                   (z / s - ((Fx / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                prim_lines.addVertices(g, g+1)

                prim_triangles.addVertices(g, g+2, g+3)
                prim_triangles.addVertices(g+3, g+2, g)

                prim_triangles.addVertices(g, g+3, g+4)
                prim_triangles.addVertices(g+4, g+3, g)

                prim_triangles.addVertices(g, g+4, g+5)
                prim_triangles.addVertices(g+5, g+4, g)

                prim_triangles.addVertices(g, g+5, g+2)
                prim_triangles.addVertices(g+2, g+5, g)

                g += 6

            if Fy != 0:

                vertices.addData3f((x)/s, (y)/s, (z)/s)
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s), (y / s - Fy / dscale), (z / s))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Fy / dscale) / (p*2))),
                                   (y / s - ((Fy / dscale) / p)),
                                   (z / s - ((Fy / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Fy / dscale) / (p*2))),
                                   (y / s - ((Fy / dscale) / p)),
                                   (z / s + ((Fy / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s + ((Fy / dscale) / (p*2))),
                                   (y / s - ((Fy / dscale) / p)),
                                   (z / s + ((Fy / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s + ((Fy / dscale) / (p*2))),
                                   (y / s - ((Fy / dscale) / p)),
                                   (z / s - ((Fy / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                prim_lines.addVertices(g, g+1)

                prim_triangles.addVertices(g, g+2, g+3)
                prim_triangles.addVertices(g+3, g+2, g)

                prim_triangles.addVertices(g, g+3, g+4)
                prim_triangles.addVertices(g+4, g+3, g)

                prim_triangles.addVertices(g, g+4, g+5)
                prim_triangles.addVertices(g+5, g+4, g)

                prim_triangles.addVertices(g, g+5, g+2)
                prim_triangles.addVertices(g+2, g+5, g)

                g += 6

            if Fz != 0:

                vertices.addData3f((x)/s, (y)/s, (z)/s)
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s), (y / s), (z / s - Fz / dscale))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Fz / dscale) / (p*2))),
                                   (y / s - ((Fz / dscale) / (p*2))),
                                   (z / s - ((Fz / dscale) / p)))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Fz / dscale) / (p*2))),
                                   (y / s + ((Fz / dscale) / (p*2))),
                                   (z / s - ((Fz / dscale) / p)))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s + ((Fz / dscale) / (p*2))),
                                   (y / s + ((Fz / dscale) / (p*2))),
                                   (z / s - ((Fz / dscale) / p)))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s + ((Fz / dscale) / (p*2))),
                                   (y / s - ((Fz / dscale) / (p*2))),
                                   (z / s - ((Fz / dscale) / p)))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                prim_lines.addVertices(g, g+1)

                prim_triangles.addVertices(g, g+2, g+3)
                prim_triangles.addVertices(g+3, g+2, g)

                prim_triangles.addVertices(g, g+3, g+4)
                prim_triangles.addVertices(g+4, g+3, g)

                prim_triangles.addVertices(g, g+4, g+5)
                prim_triangles.addVertices(g+5, g+4, g)

                prim_triangles.addVertices(g, g+5, g+2)
                prim_triangles.addVertices(g+2, g+5, g)

                g += 6

            geo_lines = Geom(vertexData)
            geo_lines.addPrimitive(prim_lines)
            geo_triangles = Geom(vertexData)
            geo_triangles.addPrimitive(prim_triangles)

            n_lines = GeomNode('PLoad_lines')
            n_lines.addGeom(geo_lines)

            n_triangles = GeomNode('PLoad_triangles')
            n_triangles.addGeom(geo_triangles)

            self.p_loadPaths.append(self.render.attachNewNode(n_lines))
            self.p_loadPaths[-1].setRenderModeThickness(1)
            self.p_loadPaths[-1].setRenderModePerspective(True)

            self.p_loadPaths.append(self.render.attachNewNode(n_triangles))
            self.p_loadPaths[-1].setRenderModeThickness(1)
            self.p_loadPaths[-1].setRenderModePerspective(True)

    def render_m_loads(self, node, scale='mm', dscale=250):
        """
        This function renders all moment loads on the system using double
        arrows pointing in the direction of the rotation axis (rotations are
        clockwise around the rotation axis when facing axis direction)
        at the scale specified in the dscale argument.
        """
        s = self.scale_check(scale)

        fmt = GeomVertexFormat.getV3n3c4()
        vertexData = GeomVertexData('m_load', fmt, Geom.UHStatic)

        p = 5
        n = 0
        Tx = 0
        Ty = 0
        Tz = 0
        if node.Tx != 0:
            n += 1
            Tx = node.Tx
        if node.Ty != 0:
            n += 1
            Ty = node.Ty
        if node.Tz != 0:
            n += 1
            Tz = node.Tz

        x = node.point.x
        y = node.point.y
        z = node.point.z

        if n > 0:

            n *= 11
            vertexData.setNumRows(n)

            vertices = GeomVertexWriter(vertexData, 'vertex')
            normals = GeomVertexWriter(vertexData, 'normal')
            colors = GeomVertexWriter(vertexData, 'color')

            prim_lines = GeomLines(Geom.UHStatic)
            prim_triangles = GeomTriangles(Geom.UHStatic)

            g = 0

            if Tx != 0:

                vertices.addData3f((x)/s, (y)/s, (z)/s)
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Tx / dscale) / p)),
                                   (y / s - ((Tx / dscale) / (p*2))),
                                   (z / s - ((Tx / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Tx / dscale) / p)),
                                   (y / s - ((Tx / dscale) / (p*2))),
                                   (z / s + ((Tx / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Tx / dscale) / p)),
                                   (y / s + ((Tx / dscale) / (p*2))),
                                   (z / s + ((Tx / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Tx / dscale) / p)),
                                   (y / s + ((Tx / dscale) / (p*2))),
                                   (z / s - ((Tx / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Tx / dscale) / p)),
                                   (y / s), (z / s))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - 2*((Tx / dscale) / p)),
                                   (y / s - ((Tx / dscale) / (p*2))),
                                   (z / s - ((Tx / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - 2*((Tx / dscale) / p)),
                                   (y / s - ((Tx / dscale) / (p*2))),
                                   (z / s + ((Tx / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - 2*((Tx / dscale) / p)),
                                   (y / s + ((Tx / dscale) / (p*2))),
                                   (z / s + ((Tx / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - 2*((Tx / dscale) / p)),
                                   (y / s + ((Tx / dscale) / (p*2))),
                                   (z / s - ((Tx / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - Tx / dscale), (y / s), (z / s))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                prim_lines.addVertices(g, g+5)
                prim_lines.addVertices(g+5, g+10)

                prim_triangles.addVertices(g, g+1, g+2)
                prim_triangles.addVertices(g+2, g+1, g)

                prim_triangles.addVertices(g, g+2, g+3)
                prim_triangles.addVertices(g+3, g+2, g)

                prim_triangles.addVertices(g, g+3, g+4)
                prim_triangles.addVertices(g+4, g+3, g)

                prim_triangles.addVertices(g, g+4, g+1)
                prim_triangles.addVertices(g+1, g+4, g)

                prim_triangles.addVertices(g+5, g+6, g+7)
                prim_triangles.addVertices(g+7, g+6, g+5)

                prim_triangles.addVertices(g+5, g+7, g+8)
                prim_triangles.addVertices(g+8, g+7, g+5)

                prim_triangles.addVertices(g+5, g+8, g+9)
                prim_triangles.addVertices(g+9, g+8, g+5)

                prim_triangles.addVertices(g+5, g+9, g+6)
                prim_triangles.addVertices(g+6, g+9, g+5)

                g += 11

            if Ty != 0:

                vertices.addData3f((x)/s, (y)/s, (z)/s)
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Ty / dscale) / (p*2))),
                                   (y / s - ((Ty / dscale) / p)),
                                   (z / s - ((Ty / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Ty / dscale) / (p*2))),
                                   (y / s - ((Ty / dscale) / p)),
                                   (z / s + ((Ty / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s + ((Ty / dscale) / (p*2))),
                                   (y / s - ((Ty / dscale) / p)),
                                   (z / s + ((Ty / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s + ((Ty / dscale) / (p*2))),
                                   (y / s - ((Ty / dscale) / p)),
                                   (z / s - ((Ty / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x)/s, (y / s - ((Ty / dscale) / p)), (z)/s)
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Ty / dscale) / (p*2))),
                                   (y / s - 2*((Ty / dscale) / p)),
                                   (z / s - ((Ty / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Ty / dscale) / (p*2))),
                                   (y / s - 2*((Ty / dscale) / p)),
                                   (z / s + ((Ty / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s + ((Ty / dscale) / (p*2))),
                                   (y / s - 2*((Ty / dscale) / p)),
                                   (z / s + ((Ty / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s + ((Ty / dscale) / (p*2))),
                                   (y / s - 2*((Ty / dscale) / p)),
                                   (z / s - ((Ty / dscale) / (p*2))))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s), (y / s - Ty / dscale), (z / s))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                prim_lines.addVertices(g, g+5)
                prim_lines.addVertices(g+5, g+10)

                prim_triangles.addVertices(g, g+1, g+2)
                prim_triangles.addVertices(g+2, g+1, g)

                prim_triangles.addVertices(g, g+2, g+3)
                prim_triangles.addVertices(g+3, g+2, g)

                prim_triangles.addVertices(g, g+3, g+4)
                prim_triangles.addVertices(g+4, g+3, g)

                prim_triangles.addVertices(g, g+4, g+1)
                prim_triangles.addVertices(g+1, g+4, g)

                prim_triangles.addVertices(g+5, g+6, g+7)
                prim_triangles.addVertices(g+7, g+6, g+5)

                prim_triangles.addVertices(g+5, g+7, g+8)
                prim_triangles.addVertices(g+8, g+7, g+5)

                prim_triangles.addVertices(g+5, g+8, g+9)
                prim_triangles.addVertices(g+9, g+8, g+5)

                prim_triangles.addVertices(g+5, g+9, g+6)
                prim_triangles.addVertices(g+6, g+9, g+5)

                g += 11

            if Tz != 0:

                vertices.addData3f((x)/s, (y)/s, (z)/s)
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Tz / dscale) / (p*2))),
                                   (y / s - ((Tz / dscale) / (p*2))),
                                   (z / s - ((Tz / dscale) / p)))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Tz / dscale) / (p*2))),
                                   (y / s + ((Tz / dscale) / (p*2))),
                                   (z / s - ((Tz / dscale) / p)))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s + ((Tz / dscale) / (p*2))),
                                   (y / s + ((Tz / dscale) / (p*2))),
                                   (z / s - ((Tz / dscale) / p)))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s + ((Tz / dscale) / (p*2))),
                                   (y / s - ((Tz / dscale) / (p*2))),
                                   (z / s - ((Tz / dscale) / p)))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x)/s, (y)/s, (z / s - ((Tz / dscale) / p)))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Tz / dscale) / (p*2))),
                                   (y / s - ((Tz / dscale) / (p*2))),
                                   (z / s - 2*((Tz / dscale) / p)))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s - ((Tz / dscale) / (p*2))),
                                   (y / s + ((Tz / dscale) / (p*2))),
                                   (z / s - 2*((Tz / dscale) / p)))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s + ((Tz / dscale) / (p*2))),
                                   (y / s + ((Tz / dscale) / (p*2))),
                                   (z / s - 2*((Tz / dscale) / p)))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s + ((Tz / dscale) / (p*2))),
                                   (y / s - ((Tz / dscale) / (p*2))),
                                   (z / s - 2*((Tz / dscale) / p)))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                vertices.addData3f((x / s), (y / s), (z / s - Tz / dscale))
                normals.addData3f(0, 0, 1)
                colors.addData4f(255, 0, 0, 1)

                prim_lines.addVertices(g, g+5)
                prim_lines.addVertices(g+5, g+10)

                prim_triangles.addVertices(g, g+1, g+2)
                prim_triangles.addVertices(g+2, g+1, g)

                prim_triangles.addVertices(g, g+2, g+3)
                prim_triangles.addVertices(g+3, g+2, g)

                prim_triangles.addVertices(g, g+3, g+4)
                prim_triangles.addVertices(g+4, g+3, g)

                prim_triangles.addVertices(g, g+4, g+1)
                prim_triangles.addVertices(g+1, g+4, g)

                prim_triangles.addVertices(g+5, g+6, g+7)
                prim_triangles.addVertices(g+7, g+6, g+5)

                prim_triangles.addVertices(g+5, g+7, g+8)
                prim_triangles.addVertices(g+8, g+7, g+5)

                prim_triangles.addVertices(g+5, g+8, g+9)
                prim_triangles.addVertices(g+9, g+8, g+5)

                prim_triangles.addVertices(g+5, g+9, g+6)
                prim_triangles.addVertices(g+6, g+9, g+5)

                g += 11

            geo_lines = Geom(vertexData)
            geo_lines.addPrimitive(prim_lines)
            geo_triangles = Geom(vertexData)
            geo_triangles.addPrimitive(prim_triangles)

            n_lines = GeomNode('MLoad_lines')
            n_lines.addGeom(geo_lines)

            n_triangles = GeomNode('MLoad_triangles')
            n_triangles.addGeom(geo_triangles)

            self.m_loadPaths.append(self.render.attachNewNode(n_lines))
            self.m_loadPaths[-1].setRenderModeThickness(1)
            self.m_loadPaths[-1].setRenderModePerspective(True)

            self.m_loadPaths.append(self.render.attachNewNode(n_triangles))
            self.m_loadPaths[-1].setRenderModeThickness(1)
            self.m_loadPaths[-1].setRenderModePerspective(True)

    def render_q_load(self, load, scale='mm', v=4, dscale=None):
        """
        This function should render distributed loads, however it is in need
        of fixing.
        """
        # TODO: Fix q loads
        s = self.scale_check(scale)

        fmt = GeomVertexFormat.getV3n3c4()
        vertexData = GeomVertexData('q_load', fmt, Geom.UHStatic)

        n = ((v+1)*2)*2
        vertexData.setNumRows(n)

        vertices = GeomVertexWriter(vertexData, 'vertex')
        normals = GeomVertexWriter(vertexData, 'normal')
        colors = GeomVertexWriter(vertexData, 'color')

        p1 = load.start_point
        p2 = load.end_point
        x = p1.x
        y = p1.y
        z = p1.z
        dx = p2.x-p1.x
        dy = p2.y-p1.y
        dz = p2.z-p1.z
        bdir = np.array([dx, dy, dz])
        bdir = bdir / np.linalg.norm(bdir)
        L = np.sqrt(dx**2 + dy**2 + dz**2)

        if dscale is None:
            dscale = L / (6 * load.start_val)

        d = -load.direction
        q = load.start_val
        p = (1/6)

        for i in range(int(n/2)):

            vertices.addData3f(((x + (dx/v)*i))/s,
                               ((y + (dy/v)*i))/s,
                               ((z + (dz/v)*i))/s)
            normals.addData3f(0, 0, 1)
            colors.addData4f(255, 0, 0, 1)

            vertices.addData3f(((x + (dx/v)*i) + q*d[0] * dscale)/s,
                               ((y + (dy/v)*i) + q*d[1] * dscale)/s,
                               ((z + (dz/v)*i) + q*d[2] * dscale)/s)
            normals.addData3f(0, 0, 1)
            colors.addData4f(255, 0, 0, 1)

            vertices.addData3f(((x + (dx/v)*i) - (q*d[0]*bdir[0]*p*dscale))/s,
                               ((y + (dy/v)*i) - (q*d[1]*bdir[1]*p*dscale))/s,
                               ((z + (dz/v)*i) - (q*d[2]*bdir[2]*p*dscale))/s)
            normals.addData3f(0, 0, 1)
            colors.addData4f(255, 0, 0, 1)

            vertices.addData3f(((x + (dx/v)*i) + (q*d[0]*bdir[0]*p*dscale))/s,
                               ((y + (dy/v)*i) + (q*d[1]*bdir[1]*p*dscale))/s,
                               ((z + (dz/v)*i) + (q*d[2]*bdir[2]*p*dscale))/s)
            normals.addData3f(0, 0, 1)
            colors.addData4f(255, 0, 0, 1)

        primitive = GeomLines(Geom.UHStatic)

        for i in range(4, n):
            if i % 4 == 0:
                primitive.addVertices(i-4, i-3)
                primitive.addVertices(i-4, i-2)
                primitive.addVertices(i-4, i-1)
                primitive.addVertices(i-3, i+1)
            elif i % n-1 == 0:
                primitive.addVertices(i-3, i)
                primitive.addVertices(i-3, i-1)
                primitive.addVertices(i-3, i-2)

        geom = Geom(vertexData)
        geom.addPrimitive(primitive)

        node = GeomNode('Beam')
        node.addGeom(geom)

        self.q_loadPaths.append(self.render.attachNewNode(node))
        self.q_loadPaths[-1].setRenderModeThickness(1)
        self.q_loadPaths[-1].setRenderModePerspective(True)

    def render_beam(self, beam, scale='mm', axis=None):
        """
        This function renders the beam passed to it as a line.
        """
        s = self.scale_check(scale)

        fmt = GeomVertexFormat.getV3n3c4()
        vertexData = GeomVertexData('beam', fmt, Geom.UHStatic)

        if axis is True:
            vertexData.setNumRows(5)
        else:
            vertexData.setNumRows(2)

        vertices = GeomVertexWriter(vertexData, 'vertex')
        normals = GeomVertexWriter(vertexData, 'normal')
        colors = GeomVertexWriter(vertexData, 'color')

        # 0
        vertices.addData3f(beam.point_1.x/s, beam.point_1.y/s,
                           beam.point_1.z/s)
        normals.addData3f(0, 0, 1)
        colors.addData4f(0, 255, 0, 1)

        # 1
        vertices.addData3f(beam.point_2.x/s, beam.point_2.y/s,
                           beam.point_2.z/s)
        normals.addData3f(0, 0, 1)
        colors.addData4f(0, 255, 0, 1)

        if axis is True:
            x, y, z = beam.elements[0].local_cs[1]
            dscale = beam.L / (4 * np.linalg.norm(x))
            vertices.addData3f((beam.point_1.x + x[0] * dscale)/s,
                               (beam.point_1.y + x[1] * dscale)/s,
                               (beam.point_1.z + x[2] * dscale)/s)
            normals.addData3f(0, 0, 1)
            colors.addData4f(255, 0, 0, 1)
            vertices.addData3f((beam.point_1.x + y[0] * dscale)/s,
                               (beam.point_1.y + y[1] * dscale)/s,
                               (beam.point_1.z + y[2] * dscale)/s)
            normals.addData3f(0, 0, 1)
            colors.addData4f(0, 255, 0, 1)
            vertices.addData3f((beam.point_1.x + z[0] * dscale)/s,
                               (beam.point_1.y + z[1] * dscale)/s,
                               (beam.point_1.z + z[2] * dscale)/s)
            normals.addData3f(0, 0, 1)
            colors.addData4f(0, 0, 255, 1)

        primitive = GeomLines(Geom.UHStatic)

        primitive.addVertex(0)
        primitive.addVertex(1)

        if axis is True:
            primitive.addVertices(0, 2)
            primitive.addVertices(0, 3)
            primitive.addVertices(0, 4)

        geom = Geom(vertexData)
        geom.addPrimitive(primitive)

        node = GeomNode('Beam')
        node.addGeom(geom)

        r_beam_path = self.render.attachNewNode(node)
        r_beam_path.setRenderModeThickness(2)
        r_beam_path.setRenderModePerspective(True)

        c_tube = CollisionTube(beam.point_1.x/s,
                               beam.point_1.y/s,
                               beam.point_1.z/s,
                               beam.point_2.x/s,
                               beam.point_2.y/s,
                               beam.point_2.z/s,
                               25/s)

        c_tube_node = CollisionNode('Beam')
        c_tube_node.addSolid(c_tube)
        c_tube_node.setIntoCollideMask(GeomNode.getDefaultCollideMask())

        c_beam_path = self.render.attachNewNode(c_tube_node)

        self.Beams.append(RenderNode(beam, r_beam_path, c_beam_path))

        """
        self.beamNodePaths.append(self.render.attachNewNode(node))
        self.beamNodePaths[-1].setRenderModeThickness(2)
        self.beamNodePaths[-1].setRenderModePerspective(True)

        c_tube = CollisionTube(beam.point_1.x/s,
                               beam.point_1.y/s,
                               beam.point_1.z/s,
                               beam.point_2.x/s,
                               beam.point_2.y/s,
                               beam.point_2.z/s,
                               25/s)

        c_tube_node = CollisionNode('Beam')

        c_tube_node.addSolid(c_tube)

        c_tube_node.setIntoCollideMask(GeomNode.getDefaultCollideMask())

        self.beamCPaths.append(self.render.attachNewNode(c_tube_node))
        """

    def render_beam_releases(self, beam, scale='mm', dscale=1):
        """
        This method renders beam rotational end releases. Translational end
        releases are not rendered by this method at this time.
        Freed rotational releases are displayed as circles in the plane
        perpendicular to the freed rotation axis.
        Fixed rotational releases are displayed as line perpendicular to the
        rotation axis.
        """
        s = self.scale_check(scale)
        ss = dscale

        color = np.array([255, 102, 0, 255]) / 255

        if beam.release is not None:

            rel = beam.release

            cs = beam.elements[0].local_cs[1]

            l_0 = int(beam.x.shape[0]/6)
            b_0 = np.array([beam.x[l_0], beam.y[l_0], beam.z[l_0]]) / s

            l_1 = int(5*beam.x.shape[0]/6)
            b_l = np.array([beam.x[l_1], beam.y[l_1], beam.z[l_1]]) / s

            fmt = GeomVertexFormat.getV3n3c4()
            vertexData = GeomVertexData('releases', fmt, Geom.UHStatic)

            vertices = GeomVertexWriter(vertexData, 'vertex')
            normals = GeomVertexWriter(vertexData, 'normal')
            colors = GeomVertexWriter(vertexData, 'color')

            prim_lines = GeomLines(Geom.UHStatic)

            n = 0

            for b, t in zip([b_0, b_l], [0, 6]):

                for k, m, o in zip([3, 4, 5], [2, 0, 1], [1, 2, 0]):

                    if rel[k + t] == 0:

                        vertices.addData3f(*(b - cs[m] * 0.5 * ss))
                        normals.addData3f(0, 0, 1)
                        colors.addData4f(*color)

                        if n != 0:
                            n += 1

                        vertices.addData3f(*(b + cs[m] * 0.5 * ss))
                        normals.addData3f(0, 0, 1)
                        colors.addData4f(*color)

                        n += 1

                        prim_lines.add_vertices(n - 1, n)

                    elif rel[k + t] == 1:

                        vertices.addData3f(*(b + cs[m]*0.5*ss))
                        normals.addData3f(0, 0, 1)
                        colors.addData4f(*color)

                        n += 1

                        for j in range(25):

                            ang = j * (2 * np.pi / 24)

                            vertices.addData3f(*(b +
                                                 cs[o]*0.5*np.sin(ang)*ss +
                                                 cs[m]*0.5*np.cos(ang)*ss))
                            normals.addData3f(0, 0, 1)
                            colors.addData4f(*color)

                            n += 1

                            prim_lines.add_vertices(n - 1, n)

            geo_lines = Geom(vertexData)
            geo_lines.addPrimitive(prim_lines)

            n_lines = GeomNode('release_lines')
            n_lines.addGeom(geo_lines)

            self.releaseNodePaths.append(self.render.attachNewNode(n_lines))
            self.releaseNodePaths[-1].setRenderModeThickness(1)
            self.releaseNodePaths[-1].setRenderModePerspective(True)

    def render_beam_profile(self, beam, scale='mm'):
        """
        If a member section is assigned to the beam, this method renders the
        section and displays the sections designation and possible stiffeners.
        """
        s = self.scale_check(scale)

        if beam.profile is not None:
            prof = beam.profile
            text = TextNode('profile')
            txt = beam.profile.designation

            if prof.stiffeners is not None:
                for stiff in prof.stiffeners:
                    if stiff.thickness > 0:
                        if stiff.direction == 'hor':
                            txt += (' s. ' + str(int(stiff.thickness)) +
                                    'mm long.')
                        else:
                            txt += (' s. ' + str(int(stiff.thickness)) +
                                    'mm')

            text.setText(txt)
            text.setTextColor(0, 0, 0, 1.0)
            text.setAlign(text.ACenter)

            self.textNodePaths.append(self.render.attachNewNode(text))

            T = beam.elements[0].trans_matrix[:3, :3]
            b_av = (np.array([np.average(beam.x),
                              np.average(beam.y),
                              np.average(beam.z)]) / s)

            self.textNodePaths[-1].setPos(*(b_av + np.dot(T.T,
                                                          np.array([0, 0,
                                                                    (prof.h /
                                                                     2) * 1.25
                                                                    / s])
                                                          )))
            self.textNodePaths[-1].setScale(2)

            fmt = GeomVertexFormat.getV3n3c4()
            vertexData = GeomVertexData('something', fmt, Geom.UHStatic)

            vertexData.setNumRows(24)

            vertices = GeomVertexWriter(vertexData, 'vertex')
            normals = GeomVertexWriter(vertexData, 'normal')
            colors = GeomVertexWriter(vertexData, 'color')

            b_0 = np.array([beam.x[0], beam.y[0], beam.z[0]]) / s
            b_l = np.array([beam.x[-1], beam.y[-1], beam.z[-1]]) / s

            # 0
            vertices.addData3f(*(b_0 + np.dot(T.T,
                                              np.array([0,
                                                        -prof.b / 2,
                                                        -prof.h / 2])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 1
            vertices.addData3f(*(b_0 + np.dot(T.T,
                                              np.array([0,
                                                        prof.b / 2,
                                                        -prof.h / 2])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 2
            vertices.addData3f(*(b_0 + np.dot(T.T,
                                              np.array([0,
                                                        prof.b / 2,
                                                        -prof.h / 2 +
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 3
            vertices.addData3f(*(b_0 + np.dot(T.T,
                                              np.array([0,
                                                        prof.t_w/2,
                                                        -prof.h / 2 +
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
            # 4
            vertices.addData3f(*(b_0 + np.dot(T.T,
                                              np.array([0,
                                                        prof.t_w/2,
                                                        prof.h / 2 -
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
            # 5
            vertices.addData3f(*(b_0 + np.dot(T.T,
                                              np.array([0,
                                                        prof.b / 2,
                                                        prof.h / 2 -
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 6
            vertices.addData3f(*(b_0 + np.dot(T.T,
                                              np.array([0,
                                                        prof.b / 2,
                                                        prof.h / 2])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 7
            vertices.addData3f(*(b_0 + np.dot(T.T,
                                              np.array([0,
                                                        -prof.b / 2,
                                                        prof.h / 2])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 8
            vertices.addData3f(*(b_0 + np.dot(T.T,
                                              np.array([0,
                                                        -prof.b / 2,
                                                        prof.h / 2 -
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 9
            vertices.addData3f(*(b_0 + np.dot(T.T,
                                              np.array([0,
                                                        - prof.t_w/2,
                                                        prof.h / 2 -
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
            # 10
            vertices.addData3f(*(b_0 + np.dot(T.T,
                                              np.array([0,
                                                        - prof.t_w/2,
                                                        - prof.h/2 +
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
            # 11
            vertices.addData3f(*(b_0 + np.dot(T.T,
                                              np.array([0,
                                                        -prof.b / 2,
                                                        - prof.h/2 +
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 12
            vertices.addData3f(*(b_l + np.dot(T.T,
                                              np.array([0,
                                                        -prof.b / 2,
                                                        -prof.h / 2])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 13
            vertices.addData3f(*(b_l + np.dot(T.T,
                                              np.array([0,
                                                        prof.b / 2,
                                                        -prof.h / 2])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 14
            vertices.addData3f(*(b_l + np.dot(T.T,
                                              np.array([0,
                                                        prof.b / 2,
                                                        -prof.h / 2 +
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 15
            vertices.addData3f(*(b_l + np.dot(T.T,
                                              np.array([0,
                                                        prof.t_w/2,
                                                        -prof.h / 2 +
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
            # 16
            vertices.addData3f(*(b_l + np.dot(T.T,
                                              np.array([0,
                                                        prof.t_w/2,
                                                        prof.h / 2 -
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
            # 17
            vertices.addData3f(*(b_l + np.dot(T.T,
                                              np.array([0,
                                                        prof.b / 2,
                                                        prof.h / 2 -
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 18
            vertices.addData3f(*(b_l + np.dot(T.T,
                                              np.array([0,
                                                        prof.b / 2,
                                                        prof.h / 2])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 19
            vertices.addData3f(*(b_l + np.dot(T.T,
                                              np.array([0,
                                                        -prof.b / 2,
                                                        prof.h / 2])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 20
            vertices.addData3f(*(b_l + np.dot(T.T,
                                              np.array([0,
                                                        -prof.b / 2,
                                                        prof.h / 2 -
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)
            # 21
            vertices.addData3f(*(b_l + np.dot(T.T,
                                              np.array([0,
                                                        - prof.t_w/2,
                                                        prof.h / 2 -
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
            # 22
            vertices.addData3f(*(b_l + np.dot(T.T,
                                              np.array([0,
                                                        - prof.t_w/2,
                                                        - prof.h/2 +
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.19607843, 0.63529412, 0.75686275, 1)
            # 23
            vertices.addData3f(*(b_l + np.dot(T.T,
                                              np.array([0,
                                                        -prof.b / 2,
                                                        - prof.h/2 +
                                                        prof.t_f])
                                              / s)))
            normals.addData3f(0, 0, 1)
            colors.addData4f(0.25882353, 0.80784314, 0.95686275, 1)

            primitive = GeomTriangles(Geom.UHStatic)

            # front
            primitive.addVertices(0, 1, 2)
            primitive.addVertices(2, 1, 0)
            primitive.addVertices(0, 2, 11)
            primitive.addVertices(11, 2, 0)
            primitive.addVertices(10, 3, 4)
            primitive.addVertices(4, 3, 10)
            primitive.addVertices(10, 4, 9)
            primitive.addVertices(9, 4, 10)
            primitive.addVertices(6, 7, 5)
            primitive.addVertices(5, 7, 6)
            primitive.addVertices(7, 8, 5)
            primitive.addVertices(5, 8, 7)

            # left bottom
            primitive.addVertices(0, 12, 23)
            primitive.addVertices(23, 12, 0)
            primitive.addVertices(0, 23, 11)
            primitive.addVertices(11, 23, 0)

            # bottom
            primitive.addVertices(0, 1, 13)
            primitive.addVertices(13, 1, 0)
            primitive.addVertices(0, 13, 12)
            primitive.addVertices(12, 13, 0)

            # bottom top left
            primitive.addVertices(11, 23, 22)
            primitive.addVertices(22, 23, 11)
            primitive.addVertices(11, 22, 10)
            primitive.addVertices(10, 22, 11)

            # bottom top right
            primitive.addVertices(3, 15, 14)
            primitive.addVertices(14, 3, 15)
            primitive.addVertices(3, 14, 2)
            primitive.addVertices(2, 14, 3)

            # web left
            primitive.addVertices(10, 22, 21)
            primitive.addVertices(21, 22, 10)
            primitive.addVertices(10, 21, 9)
            primitive.addVertices(9, 21, 10)

            # web right
            primitive.addVertices(3, 16, 15)
            primitive.addVertices(15, 16, 3)
            primitive.addVertices(3, 4, 16)
            primitive.addVertices(16, 4, 3)

            # right bottom
            primitive.addVertices(1, 14, 13)
            primitive.addVertices(13, 14, 1)
            primitive.addVertices(1, 2, 14)
            primitive.addVertices(14, 2, 1)

            # left top
            primitive.addVertices(8, 20, 19)
            primitive.addVertices(19, 20, 8)
            primitive.addVertices(8, 19, 7)
            primitive.addVertices(7, 19, 8)

            # top
            primitive.addVertices(7, 18, 19)
            primitive.addVertices(19, 18, 7)
            primitive.addVertices(7, 18, 6)
            primitive.addVertices(6, 18, 7)

            # top bottom left
            primitive.addVertices(8, 20, 9)
            primitive.addVertices(9, 20, 8)
            primitive.addVertices(20, 9, 21)
            primitive.addVertices(21, 9, 20)

            # top bottom right
            primitive.addVertices(4, 16, 17)
            primitive.addVertices(17, 16, 4)
            primitive.addVertices(5, 16, 17)
            primitive.addVertices(17, 16, 5)
            primitive.addVertices(4, 17, 5)
            primitive.addVertices(5, 17, 4)
            primitive.addVertices(5, 4, 16)
            primitive.addVertices(16, 4, 5)

            # right top
            primitive.addVertices(6, 18, 17)
            primitive.addVertices(17, 18, 6)
            primitive.addVertices(6, 17, 5)
            primitive.addVertices(5, 17, 6)

            # back
            primitive.addVertices(12, 13, 14)
            primitive.addVertices(2+12, 1+12, 0+12)
            primitive.addVertices(0+12, 2+12, 11+12)
            primitive.addVertices(11+12, 2+12, 0+12)
            primitive.addVertices(10+12, 3+12, 4+12)
            primitive.addVertices(4+12, 3+12, 10+12)
            primitive.addVertices(10+12, 4+12, 9+12)
            primitive.addVertices(9+12, 4+12, 10+12)
            primitive.addVertices(6+12, 7+12, 5+12)
            primitive.addVertices(5+12, 7+12, 6+12)
            primitive.addVertices(7+12, 8+12, 5+12)
            primitive.addVertices(5+12, 8+12, 7+12)

            geom = Geom(vertexData)
            geom.addPrimitive(primitive)

            node = GeomNode('I_frame')
            node.addGeom(geom)

            self.beamprofNodePaths.append(self.render.attachNewNode(node))

    def render_beam_displaced(self, beam, scale='mm', dscale=25):
        """
        This method renders the displaced beam as a dashed line.
        The displacement is rendered at the scale specified in the dscale
        argument.
        """
        s = self.scale_check(scale)

        fmt = GeomVertexFormat.getV3n3c4()
        vertexData = GeomVertexData('displacement', fmt, Geom.UHStatic)

        n = beam.x.shape[0]
        vertexData.setNumRows(n)

        vertices = GeomVertexWriter(vertexData, 'vertex')
        normals = GeomVertexWriter(vertexData, 'normal')
        colors = GeomVertexWriter(vertexData, 'color')

        for i in range(n):

            x = (beam.x[i] + beam.dx[i] * dscale)/s
            y = (beam.y[i] + beam.dy[i] * dscale)/s
            z = (beam.z[i] + beam.dz[i] * dscale)/s

            vertices.addData3f(x, y, z)

            normals.addData3f(0, 0, 1)
            colors.addData4f(255, 0, 0, 1.0)

            if np.abs(beam.dx[i]) == np.max(np.abs(beam.dx)):
                if not np.isclose(beam.dx[i], 0, atol=1e-1):
                    text = TextNode('max dx')
                    text.setText('Max dx: ' +
                                 str(np.round(beam.dx[i], 1)) + ' mm')
                    text.setTextColor(255, 0, 0, 1.0)
                    text.setAlign(text.ACenter)

                    self.textNodePaths.append(self.render.attachNewNode(text))
                    self.textNodePaths[-1].setPos(beam.x[i]/s, beam.y[i]/s,
                                                  beam.z[i]/s + 1)
                    self.textNodePaths[-1].setScale(1.5)

            if np.abs(beam.dy[i]) == np.max(np.abs(beam.dy)):
                if not np.isclose(beam.dy[i], 0, atol=1e-2):
                    text = TextNode('max dy')
                    text.setText('Max dy: ' +
                                 str(np.round(beam.dy[i], 1)) + ' mm')
                    text.setTextColor(255, 0, 0, 1.0)
                    text.setAlign(text.ACenter)

                    self.textNodePaths.append(self.render.attachNewNode(text))
                    self.textNodePaths[-1].setPos(beam.x[i]/s, beam.y[i]/s,
                                                  beam.z[i]/s + 3)
                    self.textNodePaths[-1].setScale(1.5)

            if np.abs(beam.dz[i]) == np.max(np.abs(beam.dz)):
                if not np.isclose(beam.dz[i], 0, atol=1e-1):
                    text = TextNode('max dz')
                    text.setText('Max dz: ' +
                                 str(np.round(beam.dz[i], 1)) + ' mm')
                    text.setTextColor(255, 0, 0, 1.0)
                    text.setAlign(text.ACenter)

                    self.textNodePaths.append(self.render.attachNewNode(text))
                    self.textNodePaths[-1].setPos(beam.x[i]/s, beam.y[i]/s,
                                                  beam.z[i]/s + 5)
                    self.textNodePaths[-1].setScale(1.5)

        primitive = GeomLines(Geom.UHStatic)

        for i in range(n):
            primitive.addVertex(i)

        geom = Geom(vertexData)
        geom.addPrimitive(primitive)

        node = GeomNode('Displaced beam')
        node.addGeom(geom)

        self.displacedbeamNodePaths.append(self.render.attachNewNode(node))
        self.displacedbeamNodePaths[-1].setRenderModeThickness(1)
        self.displacedbeamNodePaths[-1].setRenderModePerspective(True)

    def render_beam_moment(self, beam, d='y', scale='mm', dscale=None):
        """
        This method renders the beam moment line for the local direction
        (x, y or z) specified in the 'd' argument for the given beam object.
        The scale of the moment line is specified in the dscale argument.
        """
        s = self.scale_check(scale)

        d_t = d

        if d == 'x':
            M = beam.Mx
            d = - beam.elements[0].local_cs[1][1]
        elif d == 'y':
            M = beam.My
            d = beam.elements[0].local_cs[1][2]
        elif d == 'z':
            M = beam.Mz
            d = beam.elements[0].local_cs[1][1]

        if (4*np.max(np.abs(M[:, 3]))) > 0:
            if dscale is None:
                dscale = beam.L / (4*np.max(np.abs(M[:, 3])))

            fmt = GeomVertexFormat.getV3n3c4()
            vertexData = GeomVertexData('moment ' + d_t, fmt, Geom.UHStatic)

            n = M.shape[0]
            vertexData.setNumRows(n+2)

            vertices = GeomVertexWriter(vertexData, 'vertex')
            normals = GeomVertexWriter(vertexData, 'normal')
            colors = GeomVertexWriter(vertexData, 'color')

            vertices.addData3f(M[0, 0]/s, M[0, 1]/s, M[0, 2]/s)
            normals.addData3f(0, 1, 0)
            colors.addData4f(255, 0, 0, 1.0)

            for i in range(n):

                Mr = M[i, :3] + np.dot(M[i, 3], d*dscale)

                try:
                    if beam.Mr.shape == (3, ):
                        beam.Mr = np.append(beam.Mr.reshape((1, 3)),
                                            Mr.reshape((1, 3)), axis=0)
                    elif beam.Mr.shape == (n, 3):
                        beam.Mr = Mr
                    else:
                        beam.Mr = np.append(beam.Mr,
                                            Mr.reshape((1, 3)), axis=0)
                except AttributeError:
                    beam.Mr = Mr

                if beam.Mr.shape == (3,):
                    x = beam.Mr[0]/s
                    y = beam.Mr[1]/s
                    z = beam.Mr[2]/s
                    vertices.addData3f(x, y, z)
                else:
                    x = beam.Mr[i, 0]/s
                    y = beam.Mr[i, 1]/s
                    z = beam.Mr[i, 2]/s
                    vertices.addData3f(x, y, z)

                normals.addData3f(0, 1, 0)
                colors.addData4f(255, 0, 0, 1.0)

                if M[i, 3] == np.max(M[:, 3]):
                    text = TextNode('max moment')
                    text.setText(str(np.round(np.abs(M[i, 3]) / 10**6, 2)) +
                                 ' kNm')
                    text.setTextColor(255, 0, 0, 1.0)
                    text.setAlign(text.ACenter)

                    self.textNodePaths.append(self.render.attachNewNode(text))
                    self.textNodePaths[-1].setPos(x, y, z)
                    self.textNodePaths[-1].setScale(2)

                if M[i, 3] == np.min(M[:, 3]):
                    text = TextNode('max moment')
                    text.setText(str(np.round(np.abs(M[i, 3]) / 10**6, 2)) +
                                 ' kNm')
                    text.setTextColor(255, 0, 0, 1.0)
                    text.setAlign(text.ACenter)

                    self.textNodePaths.append(self.render.attachNewNode(text))
                    self.textNodePaths[-1].setPos(x, y, z)
                    self.textNodePaths[-1].setScale(2)

            vertices.addData3f(M[-1, 0]/s, M[-1, 1]/s, M[-1, 2]/s)
            normals.addData3f(0, 1, 0)
            colors.addData4f(255, 0, 0, 1.0)

            primitive = GeomLines(Geom.UHStatic)

            for i in range(1, n+2):
                primitive.addVertices(i-1, i)

            geom = Geom(vertexData)
            geom.addPrimitive(primitive)
            node = GeomNode('Moment')
            node.addGeom(geom)

            self.beammomentNodePaths.append(self.render.attachNewNode(node))
            self.beammomentNodePaths[-1].setRenderModeThickness(1)
            self.beammomentNodePaths[-1].setRenderModePerspective(True)

    def render_beam_normal_forces(self, beam, scale='mm', dscale=None):
        """
        This method renders the given beam objects normal force line at the
        scale specified in the dscale argument.
        """
        s = self.scale_check(scale)

        N = beam.N
        d = - beam.elements[0].local_cs[1][1]

        if (4*np.max(np.abs(N[:, 3]))) > 0:
            if dscale is None:
                dscale = beam.L / (4*np.max(np.abs(N[:, 3])))

            fmt = GeomVertexFormat.getV3n3c4()
            vertexData = GeomVertexData('beam', fmt, Geom.UHStatic)

            n = N.shape[0]
            vertexData.setNumRows(n+2)
            r_1 = False
            r_2 = False

            vertices = GeomVertexWriter(vertexData, 'vertex')
            normals = GeomVertexWriter(vertexData, 'normal')
            colors = GeomVertexWriter(vertexData, 'color')

            vertices.addData3f(N[0, 0]/s, N[0, 1]/s, N[0, 2]/s)
            normals.addData3f(0, 1, 0)
            colors.addData4f(255, 0, 0, 1.0)

            for i in range(n):

                Nr = N[i, :3] + np.dot(N[i, 3], d*dscale)

                try:
                    if beam.Nr.shape == (3,):
                        beam.Nr = np.append(beam.Nr.reshape((1, 3)),
                                            Nr.reshape((1, 3)), axis=0)
                    else:
                        beam.Nr = np.append(beam.Nr,
                                            Nr.reshape((1, 3)), axis=0)
                except AttributeError:
                    beam.Nr = Nr

                if beam.Nr.shape == (3,):
                    x = beam.Nr[0]/s
                    y = beam.Nr[1]/s
                    z = beam.Nr[2]/s
                    vertices.addData3f(x, y, z)
                else:
                    x = beam.Nr[i, 0]/s
                    y = beam.Nr[i, 1]/s
                    z = beam.Nr[i, 2]/s
                    vertices.addData3f(x, y, z)

                normals.addData3f(0, 1, 0)
                colors.addData4f(255, 0, 0, 1.0)

                if N[i, 3] == np.max(N[:, 3]):
                    if r_1 is False:
                        text = TextNode('max normal force')
                        text.setText(str(np.round(N[i, 3] / 10**3, 2)) +
                                     ' kN')
                        text.setTextColor(255, 0, 0, 1.0)
                        text.setAlign(text.ACenter)

                        self.textNodePaths.append((self.render
                                                   ).attachNewNode(text))
                        self.textNodePaths[-1].setPos(x, y, z)
                        self.textNodePaths[-1].setScale(2)

                        r_1 = True

                if N[i, 3] == np.min(N[:, 3]):
                    if r_2 is False:
                        text = TextNode('min normal force')
                        text.setText(str(np.round(N[i, 3] / 10**3, 2)) +
                                     ' kN')
                        text.setTextColor(255, 0, 0, 1.0)
                        text.setAlign(text.ACenter)

                        self.textNodePaths.append((self.render
                                                   ).attachNewNode(text))
                        self.textNodePaths[-1].setPos(x, y, z)
                        self.textNodePaths[-1].setScale(2)

                        r_2 = True

            vertices.addData3f(N[-1, 0]/s, N[-1, 1]/s, N[-1, 2]/s)
            normals.addData3f(0, 1, 0)
            colors.addData4f(255, 0, 0, 1.0)

            primitive = GeomLines(Geom.UHStatic)

            for i in range(1, n+2):
                primitive.addVertices(i-1, i)

            geom = Geom(vertexData)
            geom.addPrimitive(primitive)
            node = GeomNode('normal forces')
            node.addGeom(geom)

            self.beamnforcesNodePaths.append(self.render.attachNewNode(node))
            self.beamnforcesNodePaths[-1].setRenderModeThickness(1)
            self.beamnforcesNodePaths[-1].setRenderModePerspective(True)

    def render_beam_shear_forces(self, beam, d='y', scale='mm', dscale=None):
        """
        This method renders the given beam objects shear force line in its
        local direction specified with the 'd' argument (y or z) at the
        scale specified in the dscale argument.
        """
        s = self.scale_check(scale)

        if d == 'y':
            V = beam.Vy
            d = beam.elements[0].local_cs[1][1]
        elif d == 'z':
            V = beam.Vz
            d = beam.elements[0].local_cs[1][2]

        if (4*np.max(np.abs(V[:, 3]))) > 0:
            if dscale is None:
                dscale = beam.L / (4*np.max(np.abs(V[:, 3])))

            fmt = GeomVertexFormat.getV3n3c4()
            vertexData = GeomVertexData('beam', fmt, Geom.UHStatic)

            n = V.shape[0]
            vertexData.setNumRows(n+2)
            r_1 = False
            r_2 = False

            vertices = GeomVertexWriter(vertexData, 'vertex')
            normals = GeomVertexWriter(vertexData, 'normal')
            colors = GeomVertexWriter(vertexData, 'color')

            vertices.addData3f(V[0, 0]/s, V[0, 1]/s, V[0, 2]/s)
            normals.addData3f(0, 1, 0)
            colors.addData4f(255, 0, 0, 1.0)

            for i in range(n):

                Vr = V[i, :3] + np.dot(V[i, 3], d*dscale)

                try:
                    if beam.Vr.shape == (3,):
                        beam.Vr = np.append(beam.Vr.reshape((1, 3)),
                                            Vr.reshape((1, 3)), axis=0)
                    elif beam.Vr.shape == (n, 3):
                        beam.Vr = Vr
                    else:
                        beam.Vr = np.append(beam.Vr,
                                            Vr.reshape((1, 3)), axis=0)
                except AttributeError:
                    beam.Vr = Vr

                if beam.Vr.shape == (3,):
                    x = beam.Vr[0]/s
                    y = beam.Vr[1]/s
                    z = beam.Vr[2]/s
                    vertices.addData3f(x, y, z)
                else:
                    x = beam.Vr[i, 0]/s
                    y = beam.Vr[i, 1]/s
                    z = beam.Vr[i, 2]/s
                    vertices.addData3f(x, y, z)

                normals.addData3f(0, 1, 0)
                colors.addData4f(255, 0, 0, 1.0)

                if V[i, 3] == np.max(V[:, 3]):
                    if r_1 is False:
                        text = TextNode('max shear force')
                        text.setText(str(np.round(V[i, 3] / 10**3, 2)) +
                                     ' kN')
                        text.setTextColor(255, 0, 0, 1.0)
                        text.setAlign(text.ACenter)

                        self.textNodePaths.append((self.render
                                                   ).attachNewNode(text))
                        self.textNodePaths[-1].setPos(x, y, z)
                        self.textNodePaths[-1].setScale(2)

                        r_1 = True

                if V[i, 3] == np.min(V[:, 3]):
                    if r_2 is False:
                        text = TextNode('min shear force')
                        text.setText(str(np.round(V[i, 3] / 10**3, 2)) +
                                     ' kN')
                        text.setTextColor(255, 0, 0, 1.0)
                        text.setAlign(text.ACenter)

                        self.textNodePaths.append((self.render
                                                   ).attachNewNode(text))
                        self.textNodePaths[-1].setPos(x, y, z)
                        self.textNodePaths[-1].setScale(2)

                        r_2 = True

            vertices.addData3f(V[-1, 0]/s, V[-1, 1]/s, V[-1, 2]/s)
            normals.addData3f(0, 1, 0)
            colors.addData4f(255, 0, 0, 1.0)

            primitive = GeomLines(Geom.UHStatic)

            for i in range(1, n+2):
                primitive.addVertices(i-1, i)

            geom = Geom(vertexData)
            geom.addPrimitive(primitive)
            node = GeomNode('shear forces')
            node.addGeom(geom)

            self.beamnforcesNodePaths.append(self.render.attachNewNode(node))
            self.beamnforcesNodePaths[-1].setRenderModeThickness(1)
            self.beamnforcesNodePaths[-1].setRenderModePerspective(True)

    def render_beam_normal_stresses(self, beam, scale='mm', s_max=235):
        """
        This method renders the given beams absolute value of local maximal
        normal stresses using a color scale based on the maximal rendered
        stress in the s_max argument.
        """

        s = self.scale_check(scale)
        fmt = GeomVertexFormat.getV3n3c4()

        vtx_data = []
        vertices = []
        normals = []
        colors = []
        primitives = []
        geoms = []
        nodes = []

        c_map = [[0, 0, 1, 1],
                 [0, 0.5, 1, 1],
                 [0, 1, 1, 1],
                 [0, 1, 0.5, 1],
                 [0, 1, 0, 1],
                 [0.5, 1, 0, 1],
                 [1, 1, 0, 1],
                 [1, 0.5, 0, 1],
                 [1, 0, 0, 1]]

        c_l = len(c_map)

        texts = []
        sc = 0.07

        for c, i in zip(reversed(c_map), range(c_l)):

            self.framesPaths.append(DirectFrame(frameColor=c,
                                                frameSize=(-sc, sc, -sc, sc)))
            self.framesPaths[-1].setPos(-1.3, 0, (1-sc) - (sc*2)*i)

            text = TextNode('stress')
            if i == 0:
                txt = str(int(s_max)) + ' < S    (N / mm^2)'
            else:
                s_1 = str(int(s_max - (s_max / (c_l-1)) * i))
                s_2 = str(int(s_max - (s_max / (c_l-1)) * (i-1)))
                txt = s_1 + ' < S < ' + s_2

            text.setText(txt)
            text.setTextColor(0, 0, 0, 1.0)

            texts.append(self.framesPaths[-1].attachNewNode(text))

            texts[-1].setPos(1.1*sc, 0, -0.25*sc)
            texts[-1].setScale(0.8 * sc)

        for el in beam.elements:

            vtx_data.append(GeomVertexData('element', fmt, Geom.UHStatic))

            vtx_data[-1].setNumRows(2)

            vertices.append(GeomVertexWriter(vtx_data[-1], 'vertex'))
            normals.append(GeomVertexWriter(vtx_data[-1], 'normal'))
            colors.append(GeomVertexWriter(vtx_data[-1], 'color'))

            sigma = np.max(el.sigma)

            if sigma > s_max:
                n = len(c_map) - 1
                # TODO: On-screen this
                print('exceeded max stress')
            else:
                n = int((sigma / s_max) * (c_l-1))
            color = c_map[n]

            vertices[-1].addData3f(el.point_1.x/s, el.point_1.y/s,
                                   el.point_1.z/s)
            normals[-1].addData3f(0, 0, 1)
            colors[-1].addData4f(color[0], color[1], color[2], color[3])

            vertices[-1].addData3f(el.point_2.x/s, el.point_2.y/s,
                                   el.point_2.z/s)
            normals[-1].addData3f(0, 0, 1)
            colors[-1].addData4f(color[0], color[1], color[2], color[3])

            primitives.append(GeomLines(Geom.UHStatic))

            primitives[-1].add_next_vertices(2)

            primitives[-1].close_primitive()

            geoms.append(Geom(vtx_data[-1]))
            geoms[-1].addPrimitive(primitives[-1])

            nodes.append(GeomNode('normal stresses'))
            nodes[-1].addGeom(geoms[-1])

            self.normalstressNodePaths.append((self.render
                                               ).attachNewNode(nodes[-1]))
            self.normalstressNodePaths[-1].setColor(*color)
            self.normalstressNodePaths[-1].setRenderModeThickness(2)
            self.normalstressNodePaths[-1].setRenderModePerspective(True)

    def render_beam_shear_stresses(self, beam, scale='mm', s_max=235):
        """
        This method renders the given beams value of shear stresses using a
        color scale based on the maximal rendered stress in the s_max argument.
        """

        s = self.scale_check(scale)
        fmt = GeomVertexFormat.getV3n3c4()

        s_max = s_max / np.sqrt(3)

        vtx_data = []
        vertices = []
        normals = []
        colors = []
        primitives = []
        geoms = []
        nodes = []

        c_map = [[0, 0, 1, 1],
                 [0, 0.5, 1, 1],
                 [0, 1, 1, 1],
                 [0, 1, 0.5, 1],
                 [0, 1, 0, 1],
                 [0.5, 1, 0, 1],
                 [1, 1, 0, 1],
                 [1, 0.5, 0, 1],
                 [1, 0, 0, 1]]

        c_l = len(c_map)

        texts = []
        sc = 0.07

        for c, i in zip(reversed(c_map), range(c_l)):

            self.framesPaths.append(DirectFrame(frameColor=c,
                                                frameSize=(-sc, sc, -sc, sc)))
            self.framesPaths[-1].setPos(-1.3, 0, (1-sc) - (sc*2)*i)

            text = TextNode('stress')
            if i == 0:
                txt = str(int(s_max)) + ' < S    (N / mm^2)'
            else:
                s_1 = str(int(s_max - (s_max / (c_l-1)) * i))
                s_2 = str(int(s_max - (s_max / (c_l-1)) * (i-1)))
                txt = s_1 + ' < S < ' + s_2

            text.setText(txt)
            text.setTextColor(0, 0, 0, 1.0)

            texts.append(self.framesPaths[-1].attachNewNode(text))

            texts[-1].setPos(1.1*sc, 0, -0.25*sc)
            texts[-1].setScale(0.8 * sc)

        for el in beam.elements:

            vtx_data.append(GeomVertexData('element', fmt, Geom.UHStatic))

            vtx_data[-1].setNumRows(2)

            vertices.append(GeomVertexWriter(vtx_data[-1], 'vertex'))
            normals.append(GeomVertexWriter(vtx_data[-1], 'normal'))
            colors.append(GeomVertexWriter(vtx_data[-1], 'color'))

            tau = np.max(el.tau)

            if tau > s_max:
                n = len(c_map) - 1
                # TODO: On-screen this
                print('exceeded max stress')
            else:
                n = int((tau / s_max) * (c_l-1))
            color = c_map[n]

            vertices[-1].addData3f(el.point_1.x/s, el.point_1.y/s,
                                   el.point_1.z/s)
            normals[-1].addData3f(0, 0, 1)
            colors[-1].addData4f(color[0], color[1], color[2], color[3])

            vertices[-1].addData3f(el.point_2.x/s, el.point_2.y/s,
                                   el.point_2.z/s)
            normals[-1].addData3f(0, 0, 1)
            colors[-1].addData4f(color[0], color[1], color[2], color[3])

            primitives.append(GeomLines(Geom.UHStatic))

            primitives[-1].add_next_vertices(2)

            primitives[-1].close_primitive()

            geoms.append(Geom(vtx_data[-1]))
            geoms[-1].addPrimitive(primitives[-1])

            nodes.append(GeomNode('shear stresses'))
            nodes[-1].addGeom(geoms[-1])

            self.shearstressNodePaths.append((self.render
                                              ).attachNewNode(nodes[-1]))
            self.shearstressNodePaths[-1].setColor(*color)
            self.shearstressNodePaths[-1].setRenderModeThickness(2)
            self.shearstressNodePaths[-1].setRenderModePerspective(True)

    def render_beam_von_mises(self, beam, scale='mm', s_max=235):
        """
        Renders a simplified von Mises stress : Which does not occur in reality
        in the nature calculated here. The maximum shear and bending stresses
        do not occur in the same location in the beams cross section as is
        assumed in this method.
        This is only to be used for very generalized calculations.

        The color scale is based on the maximum stress specified in the s_max
        argument.
        """

        s = self.scale_check(scale)
        fmt = GeomVertexFormat.getV3n3c4()

        s_max = s_max

        vtx_data = []
        vertices = []
        normals = []
        colors = []
        primitives = []
        geoms = []
        nodes = []

        c_map = [[0, 0, 1, 1],
                 [0, 0.5, 1, 1],
                 [0, 1, 1, 1],
                 [0, 1, 0.5, 1],
                 [0, 1, 0, 1],
                 [0.5, 1, 0, 1],
                 [1, 1, 0, 1],
                 [1, 0.5, 0, 1],
                 [1, 0, 0, 1]]

        c_l = len(c_map)

        texts = []
        sc = 0.07

        for c, i in zip(reversed(c_map), range(c_l)):

            self.framesPaths.append(DirectFrame(frameColor=c,
                                                frameSize=(-sc, sc, -sc, sc)))
            self.framesPaths[-1].setPos(-1.3, 0, (1-sc) - (sc*2)*i)

            text = TextNode('stress')
            if i == 0:
                txt = str(int(s_max)) + ' < S    (N / mm^2)'
            else:
                s_1 = str(int(s_max - (s_max / (c_l-1)) * i))
                s_2 = str(int(s_max - (s_max / (c_l-1)) * (i-1)))
                txt = s_1 + ' < S < ' + s_2

            text.setText(txt)
            text.setTextColor(0, 0, 0, 1.0)

            texts.append(self.framesPaths[-1].attachNewNode(text))

            texts[-1].setPos(1.1*sc, 0, -0.25*sc)
            texts[-1].setScale(0.8 * sc)

        for el in beam.elements:

            vtx_data.append(GeomVertexData('element', fmt, Geom.UHStatic))

            vtx_data[-1].setNumRows(2)

            vertices.append(GeomVertexWriter(vtx_data[-1], 'vertex'))
            normals.append(GeomVertexWriter(vtx_data[-1], 'normal'))
            colors.append(GeomVertexWriter(vtx_data[-1], 'color'))

            tau = np.max(el.tau)
            sigma = np.max(el.sigma)

            v_mis = sigma + (tau / np.sqrt(3))

            if v_mis > s_max:
                n = len(c_map) - 1
                # TODO: On-screen this
                print('exceeded max stress')
            else:
                n = int((v_mis / s_max) * (c_l-1))
            color = c_map[n]

            vertices[-1].addData3f(el.point_1.x/s, el.point_1.y/s,
                                   el.point_1.z/s)
            normals[-1].addData3f(0, 0, 1)
            colors[-1].addData4f(color[0], color[1], color[2], color[3])

            vertices[-1].addData3f(el.point_2.x/s, el.point_2.y/s,
                                   el.point_2.z/s)
            normals[-1].addData3f(0, 0, 1)
            colors[-1].addData4f(color[0], color[1], color[2], color[3])

            primitives.append(GeomLines(Geom.UHStatic))

            primitives[-1].add_next_vertices(2)

            primitives[-1].close_primitive()

            geoms.append(Geom(vtx_data[-1]))
            geoms[-1].addPrimitive(primitives[-1])

            nodes.append(GeomNode('von Mises stresses'))
            nodes[-1].addGeom(geoms[-1])

            self.shearstressNodePaths.append((self.render
                                              ).attachNewNode(nodes[-1]))
            self.shearstressNodePaths[-1].setColor(*color)
            self.shearstressNodePaths[-1].setRenderModeThickness(2)
            self.shearstressNodePaths[-1].setRenderModePerspective(True)

    def render_node(self, node, scale='mm', dscale=.4):
        """
        This method renders the given node at the scale specified in the dscale
        argument.
        """
        s = self.scale_check(scale)

        r_node_path = self.loader.load_model('models/misc/sphere')

        r_node_path.reparentTo(self.render)
        r_node_path.setScale(dscale, dscale, dscale)
        r_node_path.setPos(node.point.x/s, node.point.y/s,
                           node.point.z/s)
        r_node_path.setColor(0, 0, 1, 1)
        r_node_path.setName('Node')

        c_sphere = CollisionSphere(node.point.x/s,
                                   node.point.y/s,
                                   node.point.z/s,
                                   dscale + 0.02)

        c_sphere_node = CollisionNode('Node')

        c_sphere_node.addSolid(c_sphere)

        c_sphere_node.setIntoCollideMask(GeomNode.getDefaultCollideMask())

        c_node_path = self.render.attachNewNode(c_sphere_node)

        self.Nodes.append(RenderNode(node, r_node_path, c_node_path))

    def render_support(self, node, scale='mm', dscale=1):
        """
        This method renders the support attached to the given node at the scale
        specified in the dscale argument.
        """
        s = self.scale_check(scale)

        if node.support is not None:

            sup = node.support

            fmt = GeomVertexFormat.getV3n3c4()
            vertexData = GeomVertexData('support', fmt, Geom.UHStatic)

            vertices = GeomVertexWriter(vertexData, 'vertex')
            normals = GeomVertexWriter(vertexData, 'normal')
            colors = GeomVertexWriter(vertexData, 'color')

            prim_lines = GeomLines(Geom.UHStatic)

            x = node.point.x / s
            y = node.point.y / s
            z = node.point.z / s

            xyz = np.array([x, y, z])

            if sup.loc_cs is None:
                cs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

            color = np.array([255, 102, 0, 255]) / 255

            ss = dscale

            n = 0

            for k, m, o in zip([0, 1, 2], [2, 0, 1], [4, 5, 3]):

                vertices.addData3f(*(xyz + cs[k] * 0.5 * ss - cs[m] * ss))
                normals.addData3f(0, 0, 1)
                colors.addData4f(*color)

                if n != 0:
                    n += 1

                vertices.addData3f(*(xyz - cs[k] * 0.5 * ss - cs[m] * ss))
                normals.addData3f(0, 0, 1)
                colors.addData4f(*color)

                n += 1

                prim_lines.add_vertices(n - 1, n)

                if sup.v[k] == 0:

                    for i in range(11):

                        vertices.addData3f(*(xyz -
                                             cs[k] * (0.5 - i * 0.1) * ss -
                                             cs[m] * ss))
                        normals.addData3f(0, 0, 1)
                        colors.addData4f(*color)

                        n += 1

                        vertices.addData3f(*(xyz -
                                             cs[k] * (0.65 - i * 0.1) * ss -
                                             cs[m] * 1.25 * ss))
                        normals.addData3f(0, 0, 1)
                        colors.addData4f(*color)

                        n += 1

                        prim_lines.add_vertices(n - 1, n)

                elif sup.v[k] is False:

                    for i in range(5):

                        vertices.addData3f(*(xyz -
                                             cs[k] * (0.5 - i * 0.25) * ss -
                                             cs[m] * ss))
                        normals.addData3f(0, 0, 1)
                        colors.addData4f(*color)

                        n += 1

                        for j in range(25):

                            ang = j * (2 * np.pi / 24)

                            vertices.addData3f(*(xyz -
                                                 cs[k] * (0.5 - i * 0.25 -
                                                          0.125 * np.sin(ang))
                                                 * ss -
                                                 cs[m] * (1.125 -
                                                          0.125 * np.cos(ang))
                                                 * ss
                                                 ))
                            normals.addData3f(0, 0, 1)
                            colors.addData4f(*color)

                            n += 1

                            prim_lines.add_vertices(n - 1, n)

                if sup.v[o] == 0:

                    vertices.addData3f(*(xyz))
                    normals.addData3f(0, 0, 1)
                    colors.addData4f(*color)

                    n += 1

                    vertices.addData3f(*(xyz - cs[m] * ss))
                    normals.addData3f(0, 0, 1)
                    colors.addData4f(*color)

                    n += 1

                    prim_lines.add_vertices(n - 1, n)

                elif sup.v[o] is False:

                    vertices.addData3f(*(xyz + cs[k] * 0.5 * ss - cs[m] * ss))
                    normals.addData3f(0, 0, 1)
                    colors.addData4f(*color)

                    n += 1

                    vertices.addData3f(*(xyz - cs[k] * 0.5 * ss - cs[m] * ss))
                    normals.addData3f(0, 0, 1)
                    colors.addData4f(*color)

                    n += 1

                    vertices.addData3f(*(xyz))
                    normals.addData3f(0, 0, 1)
                    colors.addData4f(*color)

                    n += 1

                    prim_lines.add_vertices(n-2, n)
                    prim_lines.add_vertices(n-1, n)

            geo_lines = Geom(vertexData)
            geo_lines.addPrimitive(prim_lines)

            n_lines = GeomNode('Support_lines')
            n_lines.addGeom(geo_lines)

            self.sup_NodePaths.append(self.render.attachNewNode(n_lines))
            self.sup_NodePaths[-1].setRenderModeThickness(1)
            self.sup_NodePaths[-1].setRenderModePerspective(True)

    def make_button(self, to_do=None, items=[]):
        """
        This method creates all on-screen buttons for the renderer and attaches
        the relevant rendering methods to these buttons.
        The to_do argument specifies the button the be created and the items
        argument may contain a list additional arguments to be passed to the
        function creating the actual button.
        """
        try:
            self.buttons
        except AttributeError:
            self.buttons = []
        if to_do is None:
            return
        elif to_do == 'supports':
            butt = DirectButton(text=('Show supports',
                                      'Show supports',
                                      'Show supports',
                                      'Show supports'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_supports,
                                extraArgs=items,
                                pos=(1.05, 0, 0.925))
            self.buttons.append(butt)
        elif to_do == 'beams':
            butt = DirectButton(text=('Show beams',
                                      'Show beams',
                                      'Show beams',
                                      'Show beams'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_beams,
                                extraArgs=items,
                                pos=(1.05, 0, 0.85))
            self.buttons.append(butt)
        elif to_do == 'releases':
            butt = DirectButton(text=('Show beam releases',
                                      'Show beam releases',
                                      'Show beam releases',
                                      'Show beam releases'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_beam_releases,
                                extraArgs=items,
                                pos=(1.05, 0, 0.775))
            self.buttons.append(butt)
        elif to_do == 'displacements':
            butt = DirectButton(text=('Show displaced beams',
                                      'Show displaced beams',
                                      'Show displaced beams',
                                      'Show displaced beams'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_beams_displaced,
                                extraArgs=items,
                                pos=(1.05, 0, 0.7))
            self.buttons.append(butt)
        elif to_do == 'profiles':
            butt = DirectButton(text=('Show beam profiles',
                                      'Show beam profiles',
                                      'Show beam profiles',
                                      'Show beam profiles'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_beams_profile,
                                extraArgs=items,
                                pos=(1.05, 0, 0.625))
            self.buttons.append(butt)
        elif to_do == 'sigma':
            butt = DirectButton(text=('Show normal stresses',
                                      'Show normal stresses',
                                      'Show normal stresses',
                                      'Show normal stresses'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_beams_normal_stresses,
                                extraArgs=items,
                                pos=(1.05, 0, 0.55))
            self.buttons.append(butt)
        elif to_do == 'shear':
            butt = DirectButton(text=('Show shear stresses',
                                      'Show shear stresses',
                                      'Show shear stresses',
                                      'Show shear stresses'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_beams_shear_stresses,
                                extraArgs=items,
                                pos=(1.05, 0, 0.475))
            self.buttons.append(butt)
        elif to_do == 'v_mises':
            butt = DirectButton(text=('Show simp. Von Mises stresses',
                                      'Show simp. Von Mises stresses',
                                      'Show simp. Von Mises stresses',
                                      'Show simp. Von Mises stresses'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_beams_von_mises,
                                extraArgs=items,
                                pos=(1.05, 0, 0.4))
            self.buttons.append(butt)
        elif to_do == 'moments':
            butt = DirectButton(text=('Show moments, y',
                                      'Show moments, y',
                                      'Show moments, y',
                                      'Show moments, y'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_moments,
                                extraArgs=items + ['y'],
                                pos=(1.05, 0, 0.325))
            self.buttons.append(butt)

            butt = DirectButton(text=('Show moments, z',
                                      'Show moments, z',
                                      'Show moments, z',
                                      'Show moments, z'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_moments,
                                extraArgs=items + ['z'],
                                pos=(1.05, 0, 0.25))
            self.buttons.append(butt)

            butt = DirectButton(text=('Show moments, x',
                                      'Show moments, x',
                                      'Show moments, x',
                                      'Show moments, x'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_moments,
                                extraArgs=items + ['x'],
                                pos=(1.05, 0, 0.175))
            self.buttons.append(butt)

        elif to_do == 'n_forces':
            butt = DirectButton(text=('Show normal forces',
                                      'Show normal forces',
                                      'Show normal forces',
                                      'Show normal forces'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_normal_forces,
                                extraArgs=items,
                                pos=(1.05, 0, 0.1))
            self.buttons.append(butt)
        elif to_do == 'v_forces':
            butt = DirectButton(text=('Show shear forces, z',
                                      'Show shear forces, z',
                                      'Show shear forces, z',
                                      'Show shear forces, z'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_shear_forces,
                                extraArgs=items + ['z'],
                                pos=(1.05, 0, 0.025))
            self.buttons.append(butt)

            butt = DirectButton(text=('Show shear forces, y',
                                      'Show shear forces, y',
                                      'Show shear forces, y',
                                      'Show shear forces, y'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_shear_forces,
                                extraArgs=items + ['y'],
                                pos=(1.05, 0, -0.05))
            self.buttons.append(butt)

        elif to_do == 'q_loads':
            butt = DirectButton(text=('Show q loads',
                                      'Show q loads',
                                      'Show q loads',
                                      'Show q loads'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_q_loads,
                                extraArgs=items,
                                pos=(1.05, 0, -0.825))
            self.buttons.append(butt)
        elif to_do == 'm_loads':
            butt = DirectButton(text=('Show nodal moments',
                                      'Show nodal moments',
                                      'Show nodal moments',
                                      'Show nodal moments'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_m_loads,
                                extraArgs=items,
                                pos=(1.05, 0, -0.9))
            self.buttons.append(butt)
        elif to_do == 'p_loads':
            butt = DirectButton(text=('Show nodal loads',
                                      'Show nodal loads',
                                      'Show nodal loads',
                                      'Show nodal loads'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_render_p_loads,
                                extraArgs=items,
                                pos=(1.05, 0, -0.975))
            self.buttons.append(butt)
        elif to_do == 'screenshot':
            butt = DirectButton(text=('Make screenshot',
                                      'Make screenshot',
                                      'Make screenshot',
                                      'Make screenshot'),
                                text_fg=(0, 0, 0, 1),
                                text_bg=(0, 1, 0.6, 1), scale=0.05,
                                command=self.cmd_screenshot,
                                extraArgs=items,
                                pos=(1.05, 0, -0.75))
            self.buttons.append(butt)

    def cmd_render_supports(self, nodes, scale='mm', dscale=1):
        """
        This method is called when the render supports button is clicked and
        renders all supports on the given nodes using the render_support
        method.
        """
        print('render supports button clicked')

        self.curr_render = 'supports'
        self.empty_render(nodes=None, beams=None)

        for node in nodes:
            self.render_support(node, scale=scale, dscale=dscale)

    def cmd_render_beams(self, beams, scale='mm', axis=None):
        """
        This method is called when the render beams button is clicked and
        renders all given beams using the render_beam method.
        """
        print('render beams button clicked')

        self.curr_render = 'beams'
        self.empty_render(nodes=None)

        for beam in beams:
            self.render_beam(beam, scale=scale, axis=axis)

    def cmd_render_beam_releases(self, beams, scale='mm', dscale=1):
        """
        This method is called when the render beam releases button is clicked
        and renders all releases on the given beams using the
        render_beam_releases method.
        """
        print('render beam releases button clicked')

        self.curr_render = 'releases'
        self.empty_render(nodes=None)

        for beam in beams:
            self.render_beam(beam, scale=scale)
            self.render_beam_releases(beam, scale=scale, dscale=dscale)

    def cmd_render_beams_normal_stresses(self, beams, scale='mm', s_max=235):
        """
        This method is called when the render normal stresses button is clicked
        and renders the normal stresses in the given beams using the
        render_beam_normal_stresses method.
        """
        print('render normal stresses button clicked')

        self.curr_render = 'sigma'
        self.empty_render(nodes=None)

        for beam in beams:
            self.render_beam_normal_stresses(beam, scale=scale, s_max=s_max)

    def cmd_render_beams_shear_stresses(self, beams, scale='mm', s_max=235):
        """
        This method is called when the render shear stresses button is clicked
        and renders the shear stresses in the given beams using the
        render_beam_shear_stresses method.
        """
        print('render shear stresses button clicked')

        self.curr_render = 'shear'
        self.empty_render(nodes=None)

        for beam in beams:
            self.render_beam_shear_stresses(beam, scale=scale, s_max=s_max)

    def cmd_render_beams_von_mises(self, beams, scale='mm', s_max=235):
        """
        This method is called when the render von mises stresses button is
        clicked and renders the simplified von mises stresses in the given
        beams using the render_beam_von_mises method.
        """
        print('render von mises stresses button clicked')

        self.curr_render = 'v_mises'
        self.empty_render(nodes=None)

        for beam in beams:
            self.render_beam_von_mises(beam, scale=scale, s_max=s_max)

    def cmd_render_beams_displaced(self, beams, scale='mm', dscale=25):
        """
        This method is called when the render displaced beams button is
        clicked and renders the displacement for the given beams using the
        render_beam_displaced method.
        """
        print('render displaced beams button clicked')

        self.curr_render = 'displacements'
        self.empty_render(nodes=None)

        for beam in beams:
            self.render_beam(beam, scale=scale, axis=None)
            self.render_beam_displaced(beam, scale=scale, dscale=dscale)

    def cmd_render_beams_profile(self, beams, scale='mm'):
        """
        This method is called when the render beam profiles button is
        clicked and renders the profiles for the given beams using the
        render_beam_profile method.
        """
        print('render beam profiles button clicked')

        self.curr_render = 'profiles'
        self.empty_render(nodes=None)

        for beam in beams:
            self.render_beam(beam, scale=scale, axis=None)
            self.render_beam_profile(beam, scale=scale)

    def cmd_render_q_loads(self, beams, loads, scale='mm', v=4, dscale=None):
        """
        This method is called when the render q loads button is
        clicked and should render the variable loads on the given beams using
        the render_q_load method.
        """
        # TODO: Fix q loads
        print('render q_loads button clicked')

        self.curr_render = 'q_loads'
        self.empty_render(nodes=None)

        for beam in beams:
            self.render_beam(beam, scale=scale, axis=None)
        for load in loads:
            self.render_q_load(load, scale=scale, v=v, dscale=dscale)

    def cmd_render_p_loads(self, beams, scale='mm', dscale=250):
        """
        This method is called when the render p loads button is clicked and
        renders the point loads on the given beams using the render_p_loads
        method.
        """
        print('render p loads button clicked')

        self.curr_render = 'p_loads'
        self.empty_render(nodes=None)

        for beam in beams:
            self.render_beam(beam, scale=scale, axis=None)
            for el in beam.elements:
                self.render_p_loads(el.node_1, scale=scale, dscale=dscale)
                self.render_p_loads(el.node_2, scale=scale, dscale=dscale)

    def cmd_render_m_loads(self, beams, scale='mm', dscale=250):
        """
        This method is called when the render m loads button is clicked and
        renders the moment loads on the given beams using the render_m_loads
        method.
        """
        print('render m loads button clicked')

        self.curr_render = 'm_loads'
        self.empty_render(nodes=None)

        for beam in beams:
            self.render_beam(beam, scale=scale, axis=None)
            for el in beam.elements:
                self.render_m_loads(el.node_1, scale=scale, dscale=dscale)
                self.render_m_loads(el.node_2, scale=scale, dscale=dscale)

    def cmd_render_moments(self, beams, scale='mm', dscale=250, d='y'):
        """
        This method is called when the render moments button is clicked and
        renders the moments in the given beams using the render_beam_moment
        method.
        """
        print('render moments button clicked')

        self.curr_render = 'moments'
        self.empty_render(nodes=None)

        for beam in beams:
            self.render_beam(beam, scale=scale, axis=None)
            self.render_beam_moment(beam, d=d, scale=scale, dscale=dscale)

    def cmd_render_normal_forces(self, beams, scale='mm', dscale=250):
        """
        This method is called when the render normal forces button is clicked
        and renders the normal forces in the given beams using the
        render_beam_normal_forces method.
        """
        print('render normal forces button clicked')

        self.curr_render = 'normal_forces'
        self.empty_render(nodes=None)

        for beam in beams:
            self.render_beam(beam, scale=scale, axis=None)
            self.render_beam_normal_forces(beam, scale=scale, dscale=dscale)

    def cmd_render_shear_forces(self, beams, scale='mm', dscale=250, d='y'):
        """
        This method is called when the render shear forces button is clicked
        and renders the shear forces in the given beams using the
        render_beam_shear_forces method.
        """
        print('render shear forces button clicked')

        self.curr_render = 'shear_forces'
        self.empty_render(nodes=None)

        for beam in beams:
            self.render_beam(beam, scale=scale, axis=None)
            self.render_beam_shear_forces(beam, d=d, scale=scale,
                                          dscale=dscale)

    def cmd_screenshot(self):
        """
        This method is called when the make screenshot button is clicked, or
        when the ctrl-x keybind is pressed within the program.
        The screenshot is copied to the clipboard and may be pasted in other
        software.
        """
        print('make screenshot button clicked')

        for button in self.buttons:
            button.hide()

        self.graphicsEngine.renderFrame()

        ss = self.window.getScreenshot()

        ram_arr = ss.getRamImageAs("RGBA")

        vct = np.array(ram_arr)

        x = self.get_size()[1]
        y = self.get_size()[0]

        arr = vct.reshape(x, y, 4, order='C')

        im_2 = Image.fromarray(arr.astype('uint8'), mode='RGBA')

        im_2 = im_2.transpose(Image.FLIP_TOP_BOTTOM)

        out = BytesIO()

        im_2.save(out, "BMP")

        data = out.getvalue()[14:]

        out.close()

        def send_to_clipboard(clip_type, data):
            win32clipboard.OpenClipboard()
            win32clipboard.EmptyClipboard()
            win32clipboard.SetClipboardData(clip_type, data)
            win32clipboard.CloseClipboard()

        send_to_clipboard(win32clipboard.CF_DIB, data)

        for button in self.buttons:
            button.show()

        if ss is None:
            print('Screenshot not taken')

    def empty_render(self, nodes=None, beams=1):
        """
        This method empties the render between the different items chosen to be
        rendered.
        """

        path_list = [self.normalstressNodePaths,
                     self.releaseNodePaths,
                     self.textNodePaths,
                     self.displacedbeamNodePaths,
                     self.beammomentNodePaths,
                     self.beamnforcesNodePaths,
                     self.normalstressNodePaths,
                     self.q_loadPaths,
                     self.p_loadPaths,
                     self.m_loadPaths,
                     self.beamprofNodePaths,
                     self.sup_NodePaths,
                     self.shearstressNodePaths,
                     self.framesPaths]

        if nodes is not None:
            for node in self.Nodes:
                node.r_node_path.hide()
                node.c_node_path.hide()

        if beams is not None:
            for node in self.Beams:
                node.r_node_path.hide()
                node.c_node_path.hide()

        path_list = [k for m in path_list for k in m]

        for node in path_list:
            node.removeNode()

    def scale_check(self, scale):
        """
        This method checks if the scale of the project is given in mm.
        If it is not, the scaling is set to 1.
        """
        if scale == 'mm':
            return 100
        else:
            return 1
