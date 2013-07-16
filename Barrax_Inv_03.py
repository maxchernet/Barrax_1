import numpy as np
import scipy.optimize
import scipy.ndimage.interpolation
from scipy import signal
import gdal
import gdalnumeric
#from struct import unpack
import eoldas
import os

def readGeoFile(filename):
    #Open a file
    ds = gdal.Open(filename)
    
    #get some dataset information
    print 'Driver: ', ds.GetDriver().ShortName,'/', ds.GetDriver().LongName
    print 'Size is ',ds.RasterXSize,'x',ds.RasterYSize, 'x',ds.RasterCount
    geotransform = ds.GetGeoTransform()
    if not geotransform is None:
        print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'
    img = np.zeros((ds.RasterXSize, ds.RasterYSize, ds.RasterCount), dtype=np.int16)
    r,c = np.mgrid[1:ds.RasterXSize+1:1,1:ds.RasterYSize+1:1]
    img = ds.ReadAsArray()
    return img, r, c

#Read parameters from *.params files
def read_params(f_param):
    #n = 76
    lai = np.zeros(13)   
    #for i in range(1,n):
    f = open(f_param, 'r')
    tmp_str = f.read()
    f.close()
    list_str = tmp_str.split('\n')
    for j in range(4,17):
        #lai[j-2] = -2*np.log( float(list_str[1].split()[j]) )
        lai[j-4] = float(list_str[24].split()[j])
    return lai

path_data = '/media/sf_Barrax/Landsat/2004-07-18/fields/'
path_out = '/home/max/Barrax/output_fields/'
#field = np.chararray(8)
field = ['2004-07-18_C1', \
'2004-07-18_C9',\
'2004-07-18_C10',\
'2004-07-18_G1',\
'2004-07-18_P1',\
'2004-07-18_P3',\
'2004-07-18_SF1']

conf_field = ['LS7_1_reg_c1.conf', \
'LS7_1_reg_c9.conf',\
'LS7_1_reg_c10.conf',\
'LS7_1_reg_g1.conf',\
'LS7_1_reg_p1.conf',\
'LS7_1_reg_p3.conf',\
'LS7_1_reg_sf1.conf']

gamma = [0.001, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 2, 3, 4, 8, 10]

for g in range(0,12):
    for point in range(0,7):
        f_img = path_data + field[point] + '.img'
        print f_img
        if os.path.isfile(f_img):
                (img, r, c) = readGeoFile(f_img)
                print img.shape
                # set the qa flags for each sample to 1 (good data)
                qa_flag_hires = np.ones_like (img).astype( np.float32 )
                qa_flag_hires[np.where( img <= 0.0 )] = 0.0
                img = img*0.0001
                f_state = path_data + field[point] + '.dat'
                f_param = path_out + field[point] + '_%.3f'%gamma[g] + '.params'
                f_fwd = path_out + field[point] + '_%.3f'%gamma[g] + '.fwd'
                f_prior = path_out + field[point] + '_%.3f'%gamma[g] + '.prior'
                f = open(f_state, "w")
                f.write("#PARAMETERS mask row col vza vaa sza saa 482.5 565.0 660.0 825.0 1650.0 2220.0 sd-482.5 sd-560.0 sd-662.0 sd-835.0 sd-1648.0 sd-2206.0\n")
                #print 'len=' 
                #print len(img[0,:,:].flatten())
                for i in range(len(img[0,:,:].flatten())):
                        f.write("%d %d %d %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n" % \
                        (qa_flag_hires[0,:,:].flatten()[i], r.flatten()[i], \
                        c.flatten()[i], 0.0, 0.0, 26.65, 121.4, img[0,:,:].flatten()[i], img[1,:,:].flatten()[i], img[2,:,:].flatten()[i], img[3,:,:].flatten()[i], \
                        img[4,:,:].flatten()[i], img[5,:,:].flatten()[i], 0.003, 0.004, 0.004, 0.015, 0.01, 0.006))
                f.close()
                cmd = 'eoldas --conf=config/eoldas_config.conf --conf=config/'+conf_field[point]+' --operator.obs.y.state=%s'%f_state + \
                ' --parameter.result.filename=%s'%f_param + \
                ' --operator.obs.y.result.filename=%s'%f_fwd + \
                ' --parameter.x.default='+str(gamma[g])+','+str(gamma[g])+',0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0'
                #+ \
                #' --operator.prior.y.result.filename=%s'%f_prior
                print cmd
                self = eoldas.eoldas(cmd)
                self.solve(write=True)

'''
f = open('/media/sf_Barrax/LAI_points/lai_points_reg_7x7_reg10.txt','w')
f.write("N xlai xhc rpl xkab scen xkw xkm xleafn xs1 xs2 xs3 xs4 lad\n")
for point in range(1,47):
        f_param = '/home/max/Barrax/output_reg_7x7/2004-07-18_LS7_LAI_point_%d'%point+'.params'
        if os.path.isfile(f_param):
                lai = read_params(f_param)
                print 'lai_point_%d'%point+'=', lai
                f.write('%d %.4f\n'%(point, (-2*np.log(lai[0]))))
f.close()
'''         
        
print 'Done!'
