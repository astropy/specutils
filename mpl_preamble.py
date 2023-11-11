import matplotlib
from matplotlib import rcParams

if __name__ == '__main__':
    bfsize = 18

    # matplotlib backend render fancy
    # doesn't seem to work though
    matplotlib.use("pgf")

    # dark background !!
    matplotlib.pyplot.style.use("dark_background")

    # LaTeX settings
    rcParams['text.usetex'] = True
    rcParams["pgf.texsystem"] = "pdflatex"
    rcParams["font.family"] = "serif"
    rcParams["font.serif"] = ""

    # other stuff
    rcParams['figure.figsize'] = (7 / 2.54, 7 / 2.54)
    rcParams['font.size'] = 1.0 * bfsize
    rcParams['font.weight'] = 700
    rcParams['axes.linewidth'] = 2
    rcParams['axes.labelsize'] = 0.95 * bfsize
    rcParams['axes.labelweight'] = 700
    rcParams['axes.labelpad'] = 5
    rcParams['axes.grid'] = False

    rcParams['xtick.major.width'] = 1.4
    rcParams['ytick.major.width'] = 1.4
    rcParams['xtick.major.size'] = 0.25 * bfsize
    rcParams['ytick.major.size'] = 0.25 * bfsize
    rcParams['xtick.minor.width'] = 1.0
    rcParams['ytick.minor.width'] = 1.0
    rcParams['xtick.minor.size'] = 0.18 * bfsize
    rcParams['ytick.minor.size'] = 0.18 * bfsize
    rcParams['xtick.labelsize'] = 0.95 * bfsize
    rcParams['ytick.labelsize'] = 0.95 * bfsize
    rcParams['xtick.major.pad'] = '8'
    rcParams['ytick.major.pad'] = '8'
    rcParams['ytick.right'] = 'True'

    rcParams['legend.fancybox'] = 'True'
    rcParams['legend.fontsize'] = 0.85 * bfsize
    rcParams['legend.labelspacing'] = 0.2
    rcParams['legend.handletextpad'] = 0.3
    rcParams['legend.borderpad'] = 0.4
    rcParams['legend.borderaxespad'] = 0.35
    rcParams['legend.columnspacing'] = 0.2
