import pyqtgraph as pg

orientations = {
    'horizontal': {'anchor': (0.5, 1), 'angle': 0},
    'vertical':   {'anchor': (0, 0.5), 'angle': -90}
}

class Annotation(pg.TextItem):
    ''' This class handles the drawing of a composite pyqtgraph.GraphicsItem object
        comprised of a TextItem and an ArrowItem. These are meant to generate
        spectral line ID markers on the plot surface.

    '''
    def __init__(self, text, plot_item, color=(0,0,0), orientation='horizontal'):

        self._plot_item = plot_item
        self._arrow_length = 20

        anchor = orientations[orientation]['anchor']
        angle = orientations[orientation]['angle']

        super(Annotation, self).__init__(text=text, color=color, anchor=anchor, angle=angle)

        from pyqtgraph import ArrowItem

        # some tweaking with the brush, fill, and line width parameters
        # is still needed; as is, the text looks somewhat fuzzy.
        self.arrow = ArrowItem(angle=-90,  headLen=0, tailWidth=1,
                                brush=(0,0,0), tailLen=self._arrow_length)

    def setPos(self, *args):

        # need to get rid of this reference to PyQt4
        from PyQt4 import QtCore

        # The text coordinates must be converted to screen (pixel) coordinates
        # in order to place the arrow marker a fixed distance *in pixels* from
        # the text.

        text_coord = self._plot_item.vb.mapViewToScene( QtCore.QPointF(args[0], args[1]))

        arrow_pos_x = text_coord.x()
        arrow_pos_y = text_coord.y() + self._arrow_length
        arrow_pos = self._plot_item.vb.mapSceneToView(QtCore.QPointF(arrow_pos_x, arrow_pos_y))

        # text positioning is handled by the base class directly with the
        # input coordinates. Arrow positioning uses the shifted coordinates.
        super(Annotation, self).setPos(QtCore.QPointF(args[0], args[1]))
        self.arrow.setPos(arrow_pos.x(), arrow_pos.y())

