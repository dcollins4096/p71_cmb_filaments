from starter import *
import scatter_fit
import scipy
#-----------------
#Gaussian Function
#----------------
def gauss(x, a, b, c):
        return a * np.exp(-((x**2)/(2*(b**2)))) + c
def point_selector(x,y,x0,x1,y0,y1,dx,dy):
    #F(x,y) = (y1-y0)*x + (x0-x1)*y +(x1*y0-x0*y1)
    #F==0, (x,y) is on the line.  <0 above, >0 below
    #All above or below, then there's no intersection. 
    #Sign change means that two points are on opposite sides of the line.
    c1 = y1-y0
    c2 = x0-x1
    c3 = x1*y0-x0*y1
    a   = c1*(x      )+c2*(y      )+c3
    amm = c1*(x-dx/2.)+c2*(y-dy/2.)+c3
    amp = c1*(x-dx/2.)+c2*(y+dy/2.)+c3
    app = c1*(x+dx/2.)+c2*(y+dy/2.)+c3
    apm = c1*(x+dx/2.)+c2*(y-dy/2.)+c3
    keep = amm*amp <= 0
    keep = np.logical_or( keep, amp*app <= 0 )
    keep = np.logical_or( keep, app*apm <= 0 )
    keep = np.logical_or( keep, apm*amm <= 0 ) 
    #pdb.set_trace()
    return keep
class filament():
    def __init__(self,nfilament, x_fil=None,y_fil=None, x_coords=None,y_coords=None,parent=None):
        """filament objects take a set of points in x and y, finds the centroid.  Really just a container."""
        self.x_fil=x_fil
        self.y_fil=y_fil
#        if nfilament in [48]:
#            #some filaments we want to curate.  I don't have a particularly good way to do this right now.
#            keep = x_fil < 15./256
#            self.x_fil = self.x_fil[keep]
#            self.y_fil = self.y_fil[keep]
            
        """Check for periodic jumps."""
        if (x_fil < parent.dx[0]).any() or (y_fil < parent.dx[1]).any():
            print( "Filament", nfilament, "probably has periodic wrap, its against the edge of the domain.")
        self.x_centroid = self.x_fil.sum()/self.x_fil.size
        self.y_centroid = self.y_fil.sum()/self.y_fil.size
        self.n_points = self.x_fil.size
        self.parent=parent

class perpandicular():
    def __init__(self,fil, ind, width=0.3, n_points_to_fit = 4):
        #fit = scatter_fit.scatter_fit(None,x_fil,y_fil, plot_points=False)#x_fil,y_fil)
        #slope = fit['fit'][0]
        #offset = fit['fit'][1]
        #ind = x_fil.size/2 
        spine_point_x = fil.x_fil[ind]; spine_point_y = fil.y_fil[ind]
        """the transverse line"""
        r = (spine_point_x-fil.x_fil)**2 + (spine_point_y - fil.y_fil)**2
        r_args = np.argsort(r)
        x_closest = fil.x_fil[r_args][0:n_points_to_fit]
        y_closest = fil.y_fil[r_args][0:n_points_to_fit]
        #poly fit has a problem with vertical lines.
        self.vertical_problem = np.abs((x_closest-x_closest.mean())).max() < 0.5*fil.parent.dx[0]
        self.horizontal_problem = np.abs((y_closest-y_closest.mean())).max() < 0.5*fil.parent.dx[1]
        if self.vertical_problem:
            Deltax = width
            Deltay = 0
            #y0 = max([spine_point_y - 0.5*Deltay,0.5]); y1 = min([spine_point_y + 0.5*Deltay, y.max()])
            self.y0=self.y1=spine_point_y
            self.x0 = max([spine_point_x - 0.5*Deltax,fil.parent.dx[0]]); self.x1 = min([spine_point_x + 0.5*Deltax, fil.parent.x.max()])
        else:
            #fit = scatter_fit.scatter_fit(None,x_closest,y_closest, plot_points=False)#x_fil,y_fil)
            to_plot = None; #low_res_ax
            fit = scatter_fit.scatter_fit(to_plot,x_closest,y_closest, plot_points=False)#x_fil,y_fil)
            slope = fit['fit'][0]
            offset = fit['fit'][1]
            perp_slope = -1./slope
            Theta = np.arctan(perp_slope)
            if slope < 1e-5:
                self.horizontal_problem = True
            #Spatial extent in line in X and Y
            Deltax = width * np.cos(Theta) 
            Deltay = width * np.sin(Theta) 
            self.x0 = max([spine_point_x - 0.5*Deltax,fil.parent.dx[0]]); self.x1 = min([spine_point_x + 0.5*Deltax, fil.parent.x.max()])
            self.y0 = max([spine_point_y - 0.5*Deltay,fil.parent.dx[1]]); self.y1 = min([spine_point_y + 0.5*Deltay, fil.parent.y.max()])
        #for the bounding box, we want the ordered set.  The bounding box will be used to extract subsets of data
        self.yLeft  = min([self.y0,self.y1]);  
        self.yRight = max([self.y0,self.y1]);  
        self.xLeft  = min([self.x0,self.x1]);  
        self.xRight = max([self.x0,self.x1]);  
        self.spine_point_x = spine_point_x
        self.spine_point_y = spine_point_y


class filament_tool():
    """Main image container object."""
    def __init__( self, source_file, resolution=[512,512],box_size_pc = 4.6,extraction_width_pc=0.3):
        self.extraction_width_pc=extraction_width_pc
        self.box_size_pc=box_size_pc
        self.resolution = np.array(resolution)
        self.data_nobuf = pyfits.open(source_file)[0].data
        self.data = np.zeros(resolution)
        resolution_offset = [ (R-N)/2 for R,N in zip(resolution, self.data_nobuf.shape)]
        #print "RES OFF", resolution_offset
        if resolution_offset[0]*resolution_offset[1] > 0:
            self.data[resolution_offset[0]:-resolution_offset[0],resolution_offset[1]:-resolution_offset[1]] = self.data_nobuf
        else:
            self.data = self.data_nobuf
        self.dx = 1./self.resolution
        self.Right = np.ones(2)
        self.y,self.x = np.mgrid[0.5*self.dx[0]:self.Right[0]-0.5*self.dx[0]:self.resolution[0]*1j,
                                 0.5*self.dx[1]:self.Right[0]-0.5*self.dx[1]:self.resolution[1]*1j]
        self.all_coords = None
        self.all_density = None
    def make_fig(self):
        self.image_fig = plt.figure()
        self.image_ax = self.image_fig.add_subplot(111)
        self.profile_fig = plt.figure()
        self.profile_ax = self.profile_fig.add_subplot(111)
    def make_image(self,cmap='gray', log=False):
        if log:
            to_plot = np.log10(self.data)
        else:
            to_plot = self.data
        d1=self.image_ax.imshow(to_plot, interpolation='nearest',origin='lower',cmap=cmap,extent=[0,1,0,1]) #extents are left,right,bottom,top
        self.image_fig.colorbar(d1)
        self.profile_aggregator = None
    def image_save(self,outname):
        self.image_ax.set_xlim(0,1)
        self.image_ax.set_ylim(0,1)

        self.image_fig.savefig(outname)
        print( outname)
    def profile_save(self,outname, take_log=False):
        if take_log:
            self.profile_fig.set_yscale('log')
        self.profile_ax.set_ylabel('density')
        self.profile_ax.set_xlabel('x[pc]')

        self.profile_fig.savefig(outname)
        print(outname)
    def plot_filament(self,this_filament, c='b', plot_number=None):
        self.image_ax.scatter(this_filament.x_fil, this_filament.y_fil, marker='o', c=c, linewidths=0,s=1)
        if plot_number is not None:
            self.image_ax.text(this_filament.x_centroid,this_filament.y_centroid, "%d"%plot_number, fontsize=5, color=c)

    def get_filament_points(self,nfilament):
        """subclassing should start here."""
        filament_mask = self.data == nfilament
        x_fil = self.x[ filament_mask ]
        y_fil = self.y[ filament_mask ]
        fil = filament(nfilament,x_fil,y_fil, parent=self)

        return fil

    def extract_profile(self, perp, shift_peak=True):
        #self.image_ax.plot([perp.xLeft,perp.xRight], [perp.yLeft,perp.yRight], c='g')
        #self.image_ax.plot([perp.x0,perp.x1], [perp.y0,perp.y1], c='y')
        #self.image_ax.scatter([perp.spine_point_x],[perp.spine_point_y], c='g')
        slice_x = slice(np.where(self.x[0,:]<=perp.xLeft)[0][-1], np.where(self.x[0,:]>=perp.xRight)[0][0]+1)
        slice_y = slice(np.where(self.y[:,0]<=perp.yLeft)[0][-1], np.where(self.y[:,0]>=perp.yRight)[0][0]+1)
        x_sub=self.x[slice_y,slice_x]
        y_sub=self.y[slice_y,slice_x]
        keep = point_selector(x_sub,y_sub,perp.x0,perp.x1,perp.y0,perp.y1,self.dx[0],self.dx[1])
        x_a = x_sub[keep]
        y_a = y_sub[keep]
        density = self.data[slice_y,slice_x][keep] # fullsetnar([fullset[ix,iy] for ix,iy in zip(x_i,y_i)])
        new_center_x = perp.spine_point_x; new_center_y = perp.spine_point_y
        if shift_peak:
            """ Shift the curve to match the peak.  To avoid multiple maxima, restrict to a window of 0.1 pc """
            peak_hunt_subset = (x_a - new_center_x)**2+(y_a-new_center_y)**2 < (0.05/4.6)**2 #0.05 pc in pixels
            max_coord=int( np.where(density[peak_hunt_subset]==density[peak_hunt_subset].max())[0].mean() )
            new_center_x =  x_a[peak_hunt_subset][max_coord]
            new_center_y =  y_a[peak_hunt_subset][max_coord]
        if (new_center_x-perp.spine_point_x)**2+(new_center_y-perp.spine_point_y)**2 > self.extraction_width_pc**2:
            error = "Error: Peak for filament %d index %d shift by more than 0.25 width"%(nfil, ind)
            error_list.append(error)
            print(error)
        sign_of_line   = np.sign(x_a-new_center_x)
        sign_of_line_y = np.sign(y_a-new_center_y) #this funny sign juggle takes care of vertical points
        sign_of_line[sign_of_line==0] = sign_of_line_y[sign_of_line==0]
        coordinates = sign_of_line*np.sqrt((x_a-new_center_x)**2+(y_a-new_center_y)**2) #This centers the profile on the Filament.
        return coordinates, density
        


    #-------------
    # Profile  Fit
    # ------------
    # this is what bins the data and calls the fitting functions used

    def profile_fit(self, the_x, the_y, bins=None):

        if bins is None:
            bins = np.linspace( the_x.min(), the_x.max(), 16)
        bin_means, bin_edges, binnumber=scipy.stats.binned_statistic(the_x, the_y, statistic='mean', bins=bins)
            
        popt, pconv = scipy.optimize.curve_fit( gauss, the_x, the_y)
        plt.clf()
        plt.plot( the_x, gauss(the_x, *popt))
        plt.scatter(the_x, the_y)
        plt.savefig('/home/dcollins4096/PigPen/fart.png')
        #self.plummer_fit(outname,nfil1)

     #-----------------




    #-------------
    # Gaussian Fit
    #-------------    
    def gauss_fit(self,outname,nfil2):
        global low_gauss_widths
        global high_gauss_widths
        global low_gauss_vars
        global high_gauss_vars
        global low_gauss_widths_NS
        global low_gauss_vars_NS
        global rangemin
        global psize
        global filnum
        global true_var
        global avgback
        i=0
        j=0
        xdata0 = None
        ydata0 = None
        xdata = None
        ydata = None
        xdata0 = nar(list(self.average_profile['uni_coords']))
        ydata0 = nar(list(self.average_profile['avg_dens']))
        xdata = nar(list(xdata0))
        ydata = nar(list(ydata0))
        print(ydata)
        l=0
        var = 0.0
        a0 = np.max(ydata) - np.min(ydata) # initial guess for maximum
        c0 = np.min(ydata) # initial guess for offset parameter
        nn=np.sum(ydata)

        for l in range(len(xdata)):
                var += ((((ydata[l])/(nn))*((xdata[l])**2)))

        b0 = (np.sqrt(var)) # initial guess for width of gaussian
        print( 'ymax0, sigma0, offset0= ')
        print( a0, b0, c0)
        y = gauss(xdata, a0, b0, c0)
        pmin=(psize/(2.0 * np.sqrt(2.0 * np.log(2.0))))
        pmax=(4.6/(2.0 * np.sqrt(2.0 * np.log(2.0))))
        popt, pcov = curve_fit(gauss, xdata, ydata, p0=[ a0, b0, c0],  bounds=([0.,pmin,0.], [(1.5*np.max(ydata)), pmax, np.max(ydata)]), method='trf')
        WG = (2.0 * (popt[1]) * np.sqrt(2.0 * np.log(2.0)))
        VG = ((2.0 * np.sqrt(2.0 * np.log(2.0)))**2)*(pcov[1][1])
##        print( 'ymax, sigma, offset= ')
##        print( popt)
        l=0
        var=0
        ydata2=nar((ydata)-(bgfac))
        mm=np.sum(ydata2)
        print('ydata2 is = ')
        print(ydata2)
        for l in range(len(xdata)):
                var += ((((ydata2[l])/(mm))*((xdata[l])**2)))
        print('my var =')
        print(var)
#        if (var<0.01 and var>0.0):
#          true_var = np.append(true_var, var)
#          avgback = np.append(avgback, popt[2])
         ## turn this on if you want to see the covarience matrix of the convergence parameters
        print( 'width of gaussian= ')
        print(WG)
        if ((WG > np.sqrt(0.000224*float(rangemin))) and (WG < 4.59) and (WG > (psize+0.00001)) and (WG < (filament_map.get_filament_points(nfil2).n_points*0.002995))):   #text formula (original)
            low_gauss_widths = np.append(low_gauss_widths, WG)
            low_gauss_vars = np.append(low_gauss_vars, VG)
            filnum = np.append(filnum, nfil2)
            true_var = np.append(true_var, var)
            avgback = np.append(avgback, popt[2])
        xs = np.linspace(-0.1, 0.1, 10000)
        plt.plot(xdata, ydata, 'k.')
        WG1=np.power(((WG*WG)-(0.000224*float(rangemin))),0.5)
        plt.plot( xs, gauss(xs, popt[0], popt[1], popt[2]), 'r', linewidth=1.0, label=r'$\langle W_{G} \rangle=%s$'%np.round_(WG1, decimals=3))

    #------------
    # Plummer Fit
    #------------
    def plummer_fit(self, outname,nfil2):
        global low_plum_widths
        global high_plum_widths
        global low_plum_vars
        global high_plum_vars
        global low_plum_widths_NS
        global low_plum_vars_NS
        global rangemin
        global psize
        i=0
        j=0
        xdata0 = None
        ydata0 = None
        xdata = None
        ydata = None
        xdata0 = nar(list(self.average_profile['uni_coords']))
        ydata0 = nar(list(self.average_profile['avg_dens']))
        xdata = nar(list(xdata0))
        ydata = nar(list(ydata0)) 

        l=0
        var = 0.0
#        nn=len(xdata)
        nn = np.sum(ydata)
        a0 = np.max(ydata) - np.min(ydata) # this is the initial maximum 
        c0 = np.min(ydata) # this is the initial offset parameter 
        for l in range(len(xdata)):
                var += (((ydata[l])*(xdata[l])**2)/(nn))

        b0 = np.sqrt(var) # this is the initial guess at rflat
#        d0 = 3.0 # this is the exponent in the plummer fit
        # d is not varied here it is fixed to 2, an option will be added later for ease
#        print( 'ymax0, rflat0, offset0, exponent0 ')
#        print( a0, b0, c0, d0)
        y = plummer(xdata, a0, b0, c0)
#        y = plummer(xdata, a0, b0, c0, d0)
#        pmin=psize/3.0
#        pmax=4.6/3.0
        pmin=psize/3.0
        pmax=4.6/1.5
        popt, pcov = curve_fit(plummer, xdata, ydata, p0=[a0, b0, c0], bounds=([0.,pmin,0.], [(1.5*np.max(ydata)) , pmax, np.max(ydata)]), method='trf')
#        popt, pcov = curve_fit(plummer, xdata, ydata, p0=[a0, b0, c0, d0], bounds=([0.,pmin,0.,1.99], [(1.5*np.max(ydata)) , pmax, np.max(ydata), 4.01]), method='trf')
        WP = (3.0)*(popt[1])  # p=2
        VP = (3.0**2)*(pcov[1][1])
#        WP = (1.53)*(popt[1])  # p=4
#        VP = (1.53**2)*(pcov[1][1])
#        WP = (1.5)*(4./popt[3])*(popt[1])
#        VP = (((1.5)*(4./popt[3]))**2)*(pcov[1][1])
        print('ymax, rflat, offset')
        print(popt)
        print('covarience matrix')
        print(pcov[1][1])
        print('width of plummer profile= ')
        print(WP)
        if ((WP > np.sqrt(0.000224*float(rangemin))) and (WP > (psize+0.00001)) and (WP < 4.59) and (WP < (filament_map.get_filament_points(nfil2).n_points*0.002995))):    #text formula (original)
            low_plum_widths = np.append(low_plum_widths, WP)
            low_plum_vars = np.append(low_plum_vars, VP)
        xs = np.linspace(-0.1, 0.1, 10000)
        WP1=np.power(((WP*WP)-(0.000224*float(rangemin))),0.5)
        plt.plot( xs, plummer(xs, popt[0], popt[1], popt[2]), 'b', linewidth=1.0, label=r'$\langle W_{P} \rangle=%s$'%np.round_(WP1, decimals=3))

        plt.xlabel('Distance(pc)', fontsize=16)
        plt.ylabel(r'Column Density($\Sigma_{0}$)', fontsize=16)

    #------------
