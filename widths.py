from starter import *
import filament
reload(filament)
import scatter_fit
plot_dir="/home/dcollins4096/PigPen"
def distribution(filament_name, density_name, prefix='test', resolution=[512,512], filament_slice=None):


    filament_map = filament.filament_tool( source_file = filament_name,resolution = resolution)

    filament_map.make_fig()
    filament_map.make_image(cmap='jet')
    filament_map.image_save("%s/%s_filaments.png"%(plot_dir,prefix))

    low_map = filament.filament_tool( density_name,resolution = resolution, box_size_pc = 4.6, extraction_width_pc=0.3) #used for resolution variation
    low_map.make_fig()
    low_map.make_image()
    low_map.image_save("%s/%s_image.png"%(plot_dir,prefix))
    if filament_slice is None:
        filament_slice = slice(None)
    for nfil in np.unique(filament_map.data)[filament_slice]:
        print(nfil)
        if (filament_map.get_filament_points(nfil).n_points < 17):
            print('filament too short, fail')
            continue
        this_one = filament_map.get_filament_points(nfil)
        """For actual profiles along the filament"""
        all_r=np.array([])
        all_rho=np.array([])
        for ind in  range(this_one.n_points):
            this_perp = filament.perpandicular(this_one,ind,0.3/4.6,n_points_to_fit=-1)

            if (np.size(np.where(low_map.x[0,:]<=this_perp.xLeft)[0])==0):
                print('fail 1', nfil)
                break
            elif (np.size(np.where(low_map.x[0,:]>=this_perp.xRight)[0])==0):
                print('fail 2', nfil)
                break
            elif (np.size(np.where(low_map.y[:,0]<=this_perp.yLeft)[0])==0):
                print('fail 3', nfil)
                break
            elif (np.size(np.where(low_map.y[:,0]>=this_perp.yRight)[0])==0):
                print('fail 4', nfil)
                break

            coords, density=low_map.extract_profile(this_perp)
            all_r = np.concatenate([all_r, coords])
            all_rho = np.concatenate([all_rho, density])
        else:
            low_map.profile_fit( nar(all_r), nar(all_rho))

density_name="/scratch3/dcollins/Paper49d/frbs/1_1/DD0030_density_x.fits"
filament_name='DD0030_density_x.fits_c1.up.NDskl.BRK.ASMB.fits'
distribution(filament_name=filament_name,density_name=density_name, filament_slice=None)
