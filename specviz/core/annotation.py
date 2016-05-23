from pyqtgraph import TextItem, ArrowItem

from ..third_party.qtpy.QtCore import QPointF


orientations = {
    'horizontal': {'anchor': (0.5, 1), 'angle': 0},
    'vertical':   {'anchor': (0, 0.5), 'angle': -90}
}

class LineIDMarker(TextItem):
    ''' This class handles the drawing of a composite pyqtgraph.GraphicsItem
        object comprised of a TextItem and an ArrowItem. These are used to
        generate spectral line ID markers on the plot surface.
    '''
    def __init__(self, text, plot_item, color=(0,0,0), orientation='horizontal'):

        self._plot_item = plot_item

        anchor = orientations[orientation]['anchor']
        angle = orientations[orientation]['angle']

        super(LineIDMarker, self).__init__(text=text, color=color, anchor=anchor, angle=angle)

        # some tweaking with the brush, fill, and line width parameters
        # is still needed; as is, the text looks somewhat fuzzy. Note
        # that this is drawing an arrow that points upwards, contrary
        # to what one would expect on a line ID marker. This was made
        # on purpose. It causes the arrow anchor point to be close to
        # the text item, thus making the set text-arrow to behave
        # appropriately when zoomed.

        self.arrow = ArrowItem(angle=90,  headLen=0, tailWidth=1,
                               brush=(0,0,0), tailLen=20)

    def setPos(self, *args):

        # Text coordinates must be converted to screen (pixel) coordinates in
        # order to place the arrow marker at an arbitrary but fixed distance
        # *in pixels* from the text. Note that this is necessary only in the
        # case we want to tweak the relative position of the arrow wrt the text.
        # In case no offset is desired, the arrow can be plotted directly at
        # the position defined by *args.

        text_coord = self._plot_item.vb.mapViewToScene(QPointF(args[0], args[1]))

        arrow_pos_x = text_coord.x()
        arrow_pos_y = text_coord.y() + 5

        arrow_pos = self._plot_item.vb.mapSceneToView(QPointF(arrow_pos_x, arrow_pos_y))

        # text positioning is handled by the base class directly with the
        # input coordinates. Arrow positioning uses the shifted coordinates.

        super(LineIDMarker, self).setPos(QPointF(args[0], args[1]))

        self.arrow.setPos(arrow_pos.x(), arrow_pos.y())

