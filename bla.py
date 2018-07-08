   fig=mycorner.corner(samples[iBurnIn:,iParamMin:iParamMax],
                       labels=Names[iParamMin:iParamMax],
                       range=Bounds[iParamMin:iParamMax],
                       truths=Truths[iParamMin:iParamMax],
                       truth_color='b', color='#E10014', quantiles=[0.16, 0.5, 0.84], bins=Bins[iParamMin:iParamMax],
                       smooth=2.,
                       smooth1d=1., hist_kwargs={'lw':3},
                       show_titles=True, title_fmt='.3f', title_kwargs={'size': 22},
                       label_kwargs={'fontsize': 24},
                       #                     levels=[0.68, 0.95], **{'plot_datapoints':False, 'plot_density':False, 'no_fill_contours':True})
                       levels=[0.68, 0.95], **{'plot_datapoints':False, 'plot_density':False, 'no_fill_contours':True})
