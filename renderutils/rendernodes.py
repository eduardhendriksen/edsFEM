class RenderNode():
    """
    The RenderNode class connects the structural element from the System with
    the rendernode used in panda3d to render it and the collisionnode used in
    panda3d to detect if the rendered node is clicked.
    """

    def __init__(self, structural_node, r_node_path, c_node_path):
        """
        Initializes the RenderNode class.

        structural_node: Structural element from the System class (Beam,
                         StructureNode)
        r_node_path: RenderNode path in panda3d
        c_node_path: CollisionNode path in panda3d
        """
        self.structural_node = structural_node
        self.r_node_path = r_node_path
        self.c_node_path = c_node_path
