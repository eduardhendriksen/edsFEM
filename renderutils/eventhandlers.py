from direct.showbase import DirectObject


class ClickHandler(DirectObject.DirectObject):
    """
    The ClickHandler class allows rendered structural elements to be clickable
    in the render view.
    """

    def __init__(self, showbase, scale=None):
        """
        Initializes the ClickHandler class.

        showbase: The ShowBase instance in which the render is taking place

        Also sets the left mouse key to perform the self.print_obj method and
        the right mouse key to perform the self.tooltip_obj method.
        """
        self.base = showbase
        self.accept('mouse1', self.print_obj)
        self.accept('mouse3', self.tooltip_obj)

        if scale is not None:
            self.scale = scale
        else:
            self.scale = 1

    def print_obj(self):
        """
        This method prints the properties of the structural element clicked
        in the python window.
        """
        if self.base.mouseWatcherNode.hasMouse():
            mpos = self.base.mouseWatcherNode.getMouse()

            self.base.pickerRay.setFromLens(self.base.camNode,
                                            mpos.getX(),
                                            mpos.getY())

            self.base.traverser.traverse(self.base.render)

            if self.base.myHandler.getNumEntries() > 0:

                self.base.myHandler.sortEntries()

                picked_obj = self.base.myHandler.getEntry(0).getIntoNodePath()

                if not picked_obj.isEmpty():

                    obj = [t for t in self.base.Beams + self.base.Nodes if
                           t.c_node_path == picked_obj]

                    if obj:

                        obj = obj[0]

                        print(obj.structural_node)

    def tooltip_obj(self):
        """
        This method displays the properties of the structural element clicked
        in a tooltip in the render.
        """

        if self.base.tooltipText is not None:

            if self.base.tooltipText.isHidden():

                if self.base.mouseWatcherNode.hasMouse():
                    mpos = self.base.mouseWatcherNode.getMouse()

                    self.base.pickerRay.setFromLens(self.base.camNode,
                                                    mpos.getX(),
                                                    mpos.getY())

                    self.base.traverser.traverse(self.base.render)

                    if self.base.myHandler.getNumEntries() > 0:

                        self.base.myHandler.sortEntries()

                        picked_obj = (self.base
                                      ).myHandler.getEntry(0).getIntoNodePath()

                        if not picked_obj.isEmpty():

                            obj = [t for t in self.base.Beams + self.base.Nodes
                                   if
                                   t.c_node_path == picked_obj]

                            if obj:

                                obj = obj[0]

                                self.base.render_tooltip((obj.structural_node
                                                          ).__repr__())
            else:

                self.base.hide_tooltip()
        else:

            if self.base.mouseWatcherNode.hasMouse():

                mpos = self.base.mouseWatcherNode.getMouse()

                self.base.pickerRay.setFromLens(self.base.camNode,
                                                mpos.getX(),
                                                mpos.getY())

                self.base.traverser.traverse(self.base.render)

                if self.base.myHandler.getNumEntries() > 0:

                    self.base.myHandler.sortEntries()

                    picked_obj = (self.base
                                  ).myHandler.getEntry(0).getIntoNodePath()

                    if not picked_obj.isEmpty():

                        obj = [t for t in self.base.Beams + self.base.Nodes if
                               t.c_node_path == picked_obj]

                        if obj:

                            obj = obj[0]

                            self.base.render_tooltip((obj.structural_node
                                                      ).__repr__())
