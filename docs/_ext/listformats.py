from docutils import nodes
from docutils.parsers.rst import Directive
from specutils import Spectrum1D

class ListFormats(Directive):

    def run(self):
        with open('read_formats.txt', 'w') as f:
            Spectrum1D.read.list_formats(out=f)

        file = open('read_formats.txt', 'r')

        formats = file.readlines()
        formats_as_str = ""

        for line in formats:
            formats_as_str += line

        literal_block_node = nodes.literal_block(text=formats_as_str)
        file.close()
        return [literal_block_node]

def setup(app):
    app.add_directive('listformats', ListFormats)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
