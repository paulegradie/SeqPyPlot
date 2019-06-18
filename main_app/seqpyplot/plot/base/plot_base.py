import matplotlib.pyplot as plt
import os
import matplotlib.lines as mlines

# plt.style.use('bmh')
plt.style.use('seaborn-deep')


class PlotBase(object):

    def __init__(self):
        pass

    def set_figure(self, figure_prefix, **args):
        fig = plt.figure(num=1,
                         dpi=600,
                         figsize=(7, 7),
                         edgecolor='black',
                         frameon=False,
                         )
        fig.suptitle(figure_prefix,
                     verticalalignment='top',
                     horizontalalignment='right',
                     fontsize=24
                     )
        return fig

    # def create_output_directory(self):
    #     #TODO Fix this to make it work
    #     dir_name = 'default_dir'
    #     try:
    #         dir_name = self.config_obj.get('data_directory', 'output')
    #     except:
    #         pass

    #     if not os.path.exists(dir_name):
    #         os.mkdir(dir_name)
    #     else:
    #         pass

    #     return dir_name

    def set_line(self, kwargs={'color': 'white'}):
        return mlines.Line2D([], [], **kwargs)