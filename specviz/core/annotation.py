from pyqtgraph import TextItem, ArrowItem

from ..third_party.qtpy.QtCore import Qt, QPointF
from ..third_party.qtpy.QtGui import QPolygonF, QPen, QColor


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


    def paint(self, p, *args):

        super(LineIDMarker, self).paint(p, args)

        bounding_rect = self.boundingRect()

        print ('@@@@@@     line: 70  - ', bounding_rect)

        print ('@@@@@@     line: 75  - ', bounding_rect.x())
        print ('@@@@@@     line: 76  - ', bounding_rect.y())
        print ('@@@@@@     line: 77  - ', bounding_rect.width())
        print ('@@@@@@     line: 78  - ', bounding_rect.height())


        # points = [QPointF(x_, y_) for (x_, y_) in zip(x, y)]
        # points.append(QPointF(x[-1], 0))
        # points.append(QPointF(x[0], 0))

        points = []

        # points.append(QPointF(0, 0))
        # points.append(QPointF(24, 40))
        # points.append(QPointF(bounding_rect.x(), bounding_rect.y()))
        # points.append(QPointF(bounding_rect.x() + bounding_rect.width(),
        #                       bounding_rect.y() + bounding_rect.height()))

        x = bounding_rect.x()
        y = bounding_rect.y() + bounding_rect.height() / 2.

        points.append(QPointF(x, y))
        points.append(QPointF(x, y - 20))

        polygon = QPolygonF(points)

        print ('@@@@@@     line: 95  - ', polygon)

        pen = QPen(QColor(0,0,0))

        p.setPen(pen)
        # p.setBrush(self.brush)

        p.drawPolygon(polygon)


