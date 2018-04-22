import matplotlib.pyplot as plt


class PlotBase(object):

    def __init__(self):
        pass

    def set_figure(self, figure_prefix, **args):
        fig = plt.figure(num=1,
                         dpi=600,
                         figsize=(10, 10),
                         edgecolor='black',
                         frameon=False,
                         )
        fig.suptitle(figure_prefix,
                     verticalalignment='top',
                     horizontalalignment='right',
                     fontsize=24
                     )
        return fig