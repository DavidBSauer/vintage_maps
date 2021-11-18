import logging
logger = logging.getLogger('model_maker')
logger.setLevel(logging.INFO)
handler = logging.FileHandler('model_maker.log')
handler.setLevel(logging.INFO)
logger.addHandler(handler)
import argparse
import mrcfile
import sys
import math
from skimage import measure
from PIL import Image, ImageOps
import numpy as np
from scipy import ndimage

parser = argparse.ArgumentParser(description='Create a stack of PNG images for etching into plates for a vintage style protein model')
req = parser.add_argument_group('Required information')
req.add_argument("-mp","--map_file",action='store', type=str, help="The MRC map file.",dest='map',default=None)
#parser.add_argument("-cl","--color",action='store', type=str, help="Color of opaque regions of the slices, as hex code.",dest='color',default='#000000')
parser.add_argument("-ax", "--slice_axis",action='store', type=str, help="The axis to be used for slicing.",dest='axis',default='Z')
## more sophisticated would be defining an arbitrary axis but would require resampling/interpolating
parser.add_argument("-mk", "--mask_file",action='store', type=str, help="The MRC mask file. Used for defining the region to illustrate.",dest='mask',default=None)
parser.add_argument("-sl", "--slice_num",action='store', type=int, help="The number of slices to generate.",dest='slice_num',default=10)
parser.add_argument("-lw", "--line_width",action='store', type=int, help="Line width in pixels.",dest='line_width',default=1)

parser.add_argument("-ct_min", "--contour_min",action='store', type=float, help="The minimum value contour line to show.",dest='ct_min',default=0)
parser.add_argument("-ct_max", "--contour_max",action='store', type=float, help="The maximum value contour line to show.",dest='ct_max',default=float('inf'))
parser.add_argument("-ct_num", "--contour_num",action='store', type=int, help="The number of contour lines to show.",dest='ct_num',default=10)

parser.add_argument("-pd_x", "--pad_x",action='store', type=float, help="Pad the X-dimension with empty space by this percentage.",dest='pd_x',default=10)
parser.add_argument("-pd_y", "--pad_y",action='store', type=float, help="Pad the Y-dimension with empty space by this percentage.",dest='pd_y',default=10)

parser.add_argument("-rm", "--resample_m",action='store', type=float, help="Factor for resampling the map. Computationally expensive but gives sharper images with aliasing.",dest='resample_m',default=1)
parser.add_argument("-n", "--name",action='store', type=str, help="Prefix of the output files.",dest='name',default='')

#parser.add_argument("-p", "--parallel",action='store_true', help="Run in parallel where possible.",dest='parallel',default=False)


args = parser.parse_args()
map = args.map
logger.info('Map file: '+str(map))
#color = args.color
#logger.info('Color: '+color)
mask = args.mask
logger.info('Mask file: '+str(mask))
axis = args.axis
logger.info('Axis: '+axis)
slice_num = args.slice_num
logger.info('Number of slices: '+str(slice_num))
line_width = args.line_width
logger.info('Line_width: '+str(line_width))
ct_min = args.ct_min
logger.info('Minimum contour: '+str(ct_min))
ct_max = args.ct_max
logger.info('Maximum contour: '+str(ct_max))
ct_num = args.ct_num
logger.info('Number of contour lines: '+str(ct_num))
pd_x = args.pd_x
logger.info('Padding of the X-dimension: '+str(pd_x))
pd_y = args.pd_y
logger.info('Padding of the Y-dimension: '+str(pd_y))
#parallel = args.parallel
#logger.info('Run in parallel: '+str(parallel))
resample_m = args.resample_m
logger.info('Factor for resampling map: '+str(resample_m))
name = args.name
logger.info('File name prefix: '+str(name))

#open MRC file
if not(map==None):
	try:
		map_file = mrcfile.open(map, mode='r')
	except:
		logger.info('Problem opening the map file. Quitting.')
		print('Problem opening the map file. Quitting.')
		sys.exit()
	else:
		pass

else:
	logger.info('No map file provided. Quitting.')
	print('No map file provided. Quitting.')
	sys.exit()
logger.info('Map shape: ('+', '.join([str(x) for x in map_file.data.shape])+')')
#get the voxel Size
voxel_size = (float(map_file.voxel_size.x),float(map_file.voxel_size.y),float(map_file.voxel_size.z))
logger.info('Voxel size for X dimension: '+str(voxel_size[0]))
logger.info('Voxel size for Y dimension: '+str(voxel_size[1]))
logger.info('Voxel size for Z dimension: '+str(voxel_size[2]))

## load the mask
if not(mask==None):
	try:
		mask_file = mrcfile.open(mask, mode='r')
	except:
		logger.info('Problem opening the mask file. Quitting.')
		print('Problem opening the mask file. Quitting.')
		sys.exit()
	else:
		pass
	'''
	## check that mask and map have the same dimensions
	try:
		assert (mask_file.shape == map_file.shape)
	except:
		logger.info('Mask and map must be the same shape. Quitting.')
		print('Mask and map file must be the same shape. Quitting.')
		sys.exit()
	else:
		pass
	'''

	## element-wise multiply the mask and map
	logger.info('Masking the map.')
	logger.info('This assumes the mask and model are on the same coordinate system!')
	map_file = np.multiply(map_file.data,mask_file.data)
	mask_file.close()

else:
	map_file = map_file.data

##save the modified map file
with mrcfile.new(name+'_masked_map.mrc') as mrc:
	mrc.set_data(map_file)
	vsize = mrc.voxel_size.copy()
	vsize.x = voxel_size[0]
	vsize.y = voxel_size[1]
	vsize.z = voxel_size[2]
	mrc.voxel_size = vsize

#resample the map
logger.info('Map shape: ('+', '.join([str(x) for x in map_file.shape])+')')

if resample_m >0:
	#rescale the image using splines
	resampled_map = ndimage.zoom(map_file, resample_m)

	with mrcfile.new(name+'_resampled_map.mrc') as mrc:
		mrc.set_data(resampled_map)
		#mrc.reset_header_stats()
		vsize = mrc.voxel_size.copy()
		vsize.x = voxel_size[0]/resample_m
		vsize.y = voxel_size[1]/resample_m
		vsize.z = voxel_size[2]/resample_m
		mrc.voxel_size = vsize

	voxel_size = [x/resample_m for x in voxel_size]
	logger.info('New voxel size for X dimension: '+str(vsize.x))
	logger.info('New voxel size for Y dimension: '+str(vsize.y))
	logger.info('New voxel size for Z dimension: '+str(vsize.z))
	logger.info('Resampled map shape: ('+', '.join([str(x) for x in resampled_map.shape])+')')

	map_file = resampled_map

else:
	logger.info('Rescale factor must be greater than 0. Quitting.')
	print('Rescale factor must be greater than 0. Quitting.')
	sys.exit()

## create slices along X, Y, or Z
map_file = np.array(map_file)
ct_max = np.max(map_file)

if axis.upper()=='X':
	thick_pix = map_file.shape[0]
	slices = [map_file[x,:,:] for x in range(0,map_file.shape[0],1)]
elif axis.upper()=='Y':
	thick_pix = map_file.data.shape[1]
	slices = [map_file[:,x,:] for x in range(0,map_file.shape[1],1)]
elif axis.upper()=='Z':
	thick_pix = map_file.data.shape[2]
	slices = [map_file[:,:,x] for x in range(0,map_file.shape[2],1)]
else:
	logger.info('Axis is not x, y, or z. Quitting.')
	print('Axis is not x, y, or z. Quitting.')

## for each slice calculate the contour lines
line_width = round(math.ceil((line_width-1)/2))
logger.info('Half line width: '+str(line_width))

##contour simple code
new_slices = []
for slice in slices:
	new_slice=np.zeros(slice.shape,dtype=map_file.dtype)
	to_place = list(set([(x,y) for level in [ct_min+q*(ct_max-ct_min)/ct_num for q in range(0,ct_num,1)] for contour in measure.find_contours(slice, level) for pos in contour for x in range(round(pos[0])-line_width,round(pos[0])+line_width+1,1) for y in range(round(pos[1])-line_width,round(pos[1])+line_width+1,1)]))
	for (x,y) in to_place:
		new_slice[x,y]=1
	new_slices.append(new_slice)

# find the first and last slices with a non-zero value, create sub-stack
left = np.min(np.nonzero(np.array([np.max(slice) for slice in new_slices])))
right = np.max(np.nonzero(np.array([np.max(slice) for slice in new_slices])))
new_slices = [new_slices[x] for x in range(left,right+1,1)]
slices = [slices[x] for x in range(left,right+1,1)]

#create new map file
map_file = np.array(slices)
with mrcfile.new(name+'_sub_map.mrc') as mrc:
	mrc.set_data(map_file)
	vsize = mrc.voxel_size.copy()
	vsize.x = voxel_size[0]
	vsize.y = voxel_size[1]
	vsize.z = voxel_size[2]
	mrc.voxel_size = vsize


#create new map file
map_file = np.array(new_slices)
with mrcfile.new(name+'_contoured_sub_map.mrc') as mrc:
	mrc.set_data(map_file)
	vsize = mrc.voxel_size.copy()
	vsize.x = voxel_size[0]
	vsize.y = voxel_size[1]
	vsize.z = voxel_size[2]
	mrc.voxel_size = vsize

#select slices for illustration
slice_positions = [int(map_file.shape[0]*x/(slice_num-1)) for x in range(0,slice_num-1,1)]+[int(map_file.shape[0]-1)]
logger.info('Slice positions: '+', '.join([str(x) for x in slice_positions]))
slices = [np.add((1-x%1)*map_file[math.floor(x),:,:],(x%1)*map_file[math.ceil(x),:,:]) for x in slice_positions]

## trim all slices to the extreme dimensions
#find left bound
lf_bound = np.min([np.min([x for x in np.nonzero(np.sum(slice,axis=1))]) for slice in slices])
logger.info('Left bound in pixels: '+str(lf_bound))

#find right bound
rt_bound = np.max([np.max([x for x in np.nonzero(np.sum(slice,axis=1))]) for slice in slices])
logger.info('Right bound in pixels: '+str(rt_bound))

#find upper bound
up_bound = np.min([np.min([x for x in np.nonzero(np.sum(slice,axis=0))]) for slice in slices])
logger.info('Upper bound in pixels: '+str(up_bound))

#find lower bound
lw_bound = np.max([np.max([x for x in np.nonzero(np.sum(slice,axis=0))]) for slice in slices])
logger.info('Lower bound in pixels: '+str(lw_bound))

#trim the slices
logger.info('Initial slice dimensions: ('+', '.join([str(x) for x in slices[0].shape])+')')
slices = [slice[lf_bound:rt_bound+1,:] for slice in slices]
slices = [slice[:,up_bound:lw_bound+1] for slice in slices]
logger.info('Trimmed slice dimensions: ('+', '.join([str(x) for x in slices[0].shape])+')')

#pad the X-dimension
if pd_x >0:
	logger.info('Padding X dimension by '+str(pd_x)+'%')
	pd_x = (pd_x/100)
	padding = np.zeros((round(slices[0].shape[0]*pd_x),slices[0].shape[1]))
	slices = [np.concatenate((padding,slice,padding),axis=0) for slice in slices]

#pad the Y-dimension
if pd_y > 0:
	logger.info('Padding Y dimension by '+str(pd_y)+'%')
	pd_y = (pd_y/100)
	padding = np.zeros((slices[0].shape[0],round(slices[0].shape[1]*pd_y)))
	slices = [np.concatenate((padding,slice,padding),axis=1) for slice in slices]
logger.info('Padded slice dimensions: ('+', '.join([str(x) for x in slices[0].shape])+')')

sliced_map = np.zeros((slice_positions[-1]+1,slices[0].shape[0],slices[0].shape[1]))
plates = np.zeros((slice_positions[-1]+1,slices[0].shape[0],slices[0].shape[1]))
for x in range(0,len(slice_positions),1):
		sliced_map[slice_positions[x]]=slices[x]
		#plates[slice_positions[x]]=1
		plates[slice_positions[x]]=1-slices[x]

map_file = np.array(sliced_map,dtype=map_file.dtype)
with mrcfile.new(name+'_slice_map.mrc') as mrc:
	mrc.set_data(map_file)
	vsize = mrc.voxel_size.copy()
	vsize.x = voxel_size[0]
	vsize.y = voxel_size[1]
	vsize.z = voxel_size[2]
	mrc.voxel_size = vsize

map_file = np.array(plates,dtype=map_file.dtype)
with mrcfile.new(name+'_plate_map.mrc') as mrc:
	mrc.set_data(map_file)
	vsize = mrc.voxel_size.copy()
	vsize.x = voxel_size[0]
	vsize.y = voxel_size[1]
	vsize.z = voxel_size[2]
	mrc.voxel_size = vsize
