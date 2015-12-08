from __future__ import absolute_import, division, print_function


class Controller(object):
    def __init__(self, viewer):
        self._viewer = viewer

    def create_sub_window(self):
        from bokeh.plotting import figure, output_file, show
        from bokeh.resources import CDN
        from bokeh.embed import file_html

        # output to static HTML file
        output_file("line.html")

        p = figure(plot_width=400, plot_height=400, responsive=True)

        # add a circle renderer with a size, color, and alpha
        p.circle([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], size=20, color="navy", alpha=0.5)

        html = file_html(p, CDN, "my plot")

        new_sub_window, web_view = self._viewer.add_sub_window()
        web_view.setHtml(html)
