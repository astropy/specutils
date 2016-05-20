import pyqtgraph as pg

orientations = {
    'horizontal': {'anchor': (0.5, 1), 'angle': 0},
    'vertical':   {'anchor': (0, 0.5), 'angle': -90}
}

class Annotation(pg.TextItem):

    def __init__(self, text, color=(0,0,0), orientation='horizontal'):

        anchor = orientations[orientation]['anchor']
        angle = orientations[orientation]['angle']

        super(Annotation, self).__init__(text=text, color=color, anchor=anchor, angle=angle)


