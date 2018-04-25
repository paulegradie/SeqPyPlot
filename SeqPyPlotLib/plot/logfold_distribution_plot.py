   from base.plot_base import PlotBase




   class LogfoldDistPlot(PlotBase):

       def __init__(self, config_obj, container, analyzer):

           self.analyzer = analyzer
           self.container = container
           self.config = config_obj


        def 
       
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
    # def single_log_plots(self):

    #     data_list = self.analyzer.foldchange_map
    #     histogram_list = dict()
    #     setup = True
    #     for time in enumerate(self.args.time):
    #         timeidx = time[0]
    #         timekey = time[1]

    #         for gene in data_list.values():
    #             if setup is True:
    #                 histogram_list[timekey] = []
    #                 setup = False
    #             else:
    #                 if isinstance(gene[timeidx], int) or isinstance(gene[timeidx], float):
    #                     if abs(float(gene[timeidx])) < self.args.log:
    #                         pass
    #                     else:
    #                         histogram_list[timekey].append(gene[timeidx])
    #         setup = True

    #     counter = 0
    #     sublist = []
    #     figure_labels = []
    #     for name in self.args.time:
    #         counter += 1
    #         sublist.append(name)
    #         if counter == 4:
    #             figure_labels.append(sublist)
    #             counter = 0
    #             sublist = []
    #     if len(sublist) != 0:
    #         figure_labels.append(sublist)

    #     filecnt = 1
    #     fig_pos = 0

    #     for figure in figure_labels:

    #         n_bins = 100
    #         color = 'black'

    #         fig, axes = plt.subplots(nrows=2, ncols=2)
    #         ax0, ax1, ax2, ax3 = axes.flatten()
    #         try:
    #             ax0.hist(sorted(histogram_list[figure[0]]),
    #                      n_bins,
    #                      color=color,
    #                      range=(min(histogram_list[figure[0]]),
    #                             max(histogram_list[figure[0]])))
    #         except ValueError:
    #             pass

    #         ax0.set_title(figure_labels[fig_pos][0])

    #         if len(figure) > 1:
    #             try:
    #                 ax1.hist(sorted(histogram_list[figure[1]]),
    #                          n_bins,
    #                          color=color,
    #                          range=(min(histogram_list[figure[1]]),
    #                                 max(histogram_list[figure[1]])))
    #             except ValueError:
    #                 pass
    #             ax1.set_title(figure_labels[fig_pos][1])

    #         if len(figure) > 2:
    #             try:
    #                 ax2.hist(sorted(histogram_list[figure[2]]),
    #                          n_bins,
    #                          color=color,
    #                          range=(min(histogram_list[figure[2]]),
    #                                 max(histogram_list[figure[2]])))
    #             except ValueError:
    #                 pass
    #             ax2.set_title(figure_labels[fig_pos][2])

    #         if len(figure) > 3:
    #             try:
    #                 ax3.hist(sorted(histogram_list[figure[3]]),
    #                          n_bins,
    #                          color=color,
    #                          range=(min(histogram_list[figure[3]]),
    #                                 max(histogram_list[figure[3]])))
    #             except ValueError:
    #                 pass

    #             ax3.set_title(figure_labels[fig_pos][3])

    #         for ax in axes.flatten():
    #             ax.spines['top'].set_visible(False)
    #             ax.spines['right'].set_visible(False)
    #             ax.get_xaxis().tick_bottom()
    #             ax.get_yaxis().tick_left()
    #             ax.set_xlabel("Log2Fold", fontsize=8)
    #             ax.set_ylabel("No. of Geens", fontsize=8)
    #             ax.set_xlim([-6,6])
    #         fig.tight_layout()
    #         path = os.path.join('.', self.args.out, self.args.prefix)
    #         plt.savefig("{}_{}_{}.png".format(path,
    #                                           str(filecnt),
    #                                           'sample_log2fold_histograms'),
    #                     format='png',
    #                     bbox_inches='tight')
    #         filecnt += 1
    #         fig_pos += 1
    #         plt.close()
    #         if fig_pos == len(figure_labels):
    #             return
