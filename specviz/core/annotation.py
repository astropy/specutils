from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from pyqtgraph import functions, TextItem

from qtpy.QtCore import QPointF
from qtpy.QtGui import QPolygonF, QPen, QColor


orientations = {
    'horizontal': {'anchor': (0.5, 1), 'angle': 0},
    'vertical':   {'anchor': (0, 0.5), 'angle': -90}
}

class LineIDMarker(TextItem):
    ''' This class handles the drawing of a modified TextItem that's
        augmented with a linear vertical marker. These items are used
        to generate spectral line ID markers on the plot surface.

        Note the convoluted handling of the 'color' parameter. This is
        due to a bug in pyqtgraph's function 'functions.mkColor', which
        bombs when presented with an argument of type Qt.GlobalColor.
    '''
    def __init__(self, text, plot_item, color=(0,0,0), orientation='horizontal'):

        self._plot_item = plot_item
        self._orientation = orientation
        self._color = functions.mkColor(color)

        anchor = orientations[orientation]['anchor']
        angle = orientations[orientation]['angle']

        super(LineIDMarker, self).__init__(text=text, color=color, anchor=anchor, angle=angle)

    def paint(self, p, *args):
        ''' Overrides the default implementation so as
            to draw a vertical marker.
        '''
        # draw the text
        super(LineIDMarker, self).paint(p, args)

        # Add marker. Geometry depends on the
        # text being vertical or horizontal.
        points = []
        bounding_rect = self.boundingRect()

        if self._orientation == 'vertical':
            x = bounding_rect.x()
            y = bounding_rect.y() + bounding_rect.height() / 2.

            points.append(QPointF(x, y))
            points.append(QPointF(x - 20, y))
        else:
            x = bounding_rect.x() + bounding_rect.width() / 2.
            y = bounding_rect.y() + bounding_rect.height() * 2.

            points.append(QPointF(x, y))
            points.append(QPointF(x, y - 20))

        polygon = QPolygonF(points)

        pen = QPen(QColor(functions.mkColor(self._color)))
        p.setPen(pen)
        p.drawPolygon(polygon)


