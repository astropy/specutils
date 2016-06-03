from ..third_party.qtpy.QtCore import QThread, pyqtSignal
import os
import logging

from ..core.data import Data


class FileLoadThread(QThread):
    status = pyqtSignal(str, int)
    result = pyqtSignal(Data)

    def __init__(self, parent=None, ):
        super(FileLoadThread, self).__init__(parent)
        self.file_name = ""
        self.file_filter = ""

    def __call__(self, file_name, file_filter):

        self.file_name = file_name
        self.file_filter = file_filter

    def run(self):
        self.status.emit("Loading file...", 0)
        data = self.read_file(self.file_name, self.file_filter)
        self.status.emit("File loaded successfully!", 1000)
        self.result.emit(data)

    def read_file(self, file_name, file_filter):
        """
        Convenience method that directly reads a spectrum from a file.
        This exists mostly to facilitate development workflow. In time it
        could be augmented to support fancier features such as wildcards,
        file lists, mixed file types, and the like.
        Note that the filter string is hard coded here; its details might
        depend on the intrincacies of the registries, loaders, and data
        classes. In other words, this is brittle code.
        """
        file_name = str(file_name)
        file_ext = os.path.splitext(file_name)[-1]

        if file_filter is None:
            if file_ext in ('.txt', '.dat'):
                file_filter = 'ASCII (*.txt *.dat)'
            else:
                file_filter = 'Generic Fits (*.fits *.mits)'

        try:
            data = Data.read(file_name, file_filter)
            return data
        except:
            logging.error("Incompatible loader for selected data.")
            return
