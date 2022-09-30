import numpy as np
nar = np.array
def scatter_fit(plt,xin,yin,use_scatter = 'scatter',add_fit_label=False,label=None,c='b',fit_range=None,
                text_loc=None, verbose=True, plot_points=True,plot_fit=True,plot_text=False,**kwargs):
    """Plots a scatter plot of *yin* vs *xin*.  Overplots a linear power law fit.
    *c* gives the color argument to both scatter and line.
    *label* gives the label to the scatter plot.  Line plot gets no label
    *use_scatter* = scatter, lscatter, or None"""
    #pdb.set_trace()
    xin = nar(xin)
    yin = nar(yin)
    args = np.argsort(xin)
    x = xin[args]
    y = yin[args]

    if len(c) == len(x):
        cused = nar(c)[args]
        fit_color = 'k'
    else:
        cused = c
        fit_color = cused
    plot = None
    if False:
        if fit_range is not None:
            OK = np.logical_and( x > fit_range[0], x < fit_range[1] )
            OK = np.logical_and( OK, y != 0 )
            OK = np.logical_and( OK, x != 0 )

            TheX = x[OK]
            TheY = y[OK]
        else:
            OK = np.logical_and( x!=0, y != 0 )
            TheX = x[OK]
            TheY = y[OK]
            
        fit = np.polyfit(np.log10(TheX), np.log10(TheY), 1)
        #stat(TheX, "Stat, TheX, ScatterFit")
        index = fit[0]
        offset = fit[1]
        x_0 = TheX[0]
        x_1 = TheX[-1]
        y_0 = 10**(offset) * x_0 ** (index)
        y_1 = 10**(offset) * x_1 ** (index)
        if text_loc is not None:
            text_x = text_loc[0]
            text_y = text_loc[1]
        else:
            text_x,text_y = x_1,y_1
        if plt is not None and plot_text == True:
            plt.text(text_x,text_y,"m = %0.2f"%index)
    else:
        fit = np.polyfit(x, y, 1)
        m = fit[0]
        b = fit[1]
        x_0 = x[0]
        x_1 = x[-1]
        y_0 = m * x_0 + b 
        y_1 = m * x_1 + b 
        index=m
        offset=b
        if text_loc is not None:
            text_x = text_loc[0]
            text_y = text_loc[1]
        else:
            text_x,text_y = x_1,y_1
        if plt is not None and plot_text == True:
            plt.text(text_x,text_y,"m = %0.2f"%index)
    if verbose:
        #print "x: %0.3g, %0.3g y: %0.3g, %0.3g m: %0.3g b: %0.3g"%(x_0,x_1,y_0,y_1,index,offset)
        label_to_use = None
    if add_fit_label:
        if label != None:
            label_to_use = label  

        label_to_use  += r'$\ m=%0.2f\ b=%0.2f$'%(index,offset)
    if plt is not None and plot_points == True:
        if use_scatter == 'lscatter':
            plot = lscatter(x,y,label=label_to_use,c=cused)
        elif use_scatter == 'scatter':
            plot = plt.scatter(x,y, label=label_to_use,c=cused,**kwargs)

    if plt is not None and plot_fit == True:
        print( "FU", [x_0,x_1],[y_0,y_1])
        plt.plot([x_0,x_1],[y_0,y_1], c=fit_color)#,label=None)

    return {'plot':plot,'fit':fit, 'x':[x_0,x_1],'y':[y_0,y_1]}
