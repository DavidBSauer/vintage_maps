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
from fpdf import FPDF
import os

parser = argparse.ArgumentParser(description='Create a stack of PNG images for etching into plates for a vintage style protein model')
req = parser.add_argument_group('Required information')
req.add_argument("-mp","--map_file",action='store', type=str, help="The MRC map file.",dest='map',default=None)
req.add_argument("-st","--slice_thick",action='store', type=float, help="The real-world thickness of each slice, in mm.",dest='thick',default=None)
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

parser.add_argument("-hm_x", "--hole_mark_x",action='store', type=float, help="Distance of the hole marks from the X-axis edges (in mm).",dest='hm_x',default=0)
parser.add_argument("-hm_y", "--hole_mark_y",action='store', type=float, help="Distance of the hole marks from the Y-axis edges (in mm).",dest='hm_y',default=0)
parser.add_argument("-hm_w", "--hole_mark_width",action='store', type=int, help="Width of hole marks (in mm).",dest='hm_w',default=0)

parser.add_argument("-rm", "--resample_m",action='store', type=float, help="Factor for resampling the map. Computationally expensive but gives sharper images with aliasing.",dest='resample_m',default=1)
parser.add_argument("-ri", "--resample_i",action='store', type=float, help="Factor for resampling the images. Fast and removes aliasing, but creates blurier images.",dest='resample_i',default=1)
parser.add_argument("-n", "--name",action='store', type=str, help="Prefix of the output files.",dest='name',default='')

parser.add_argument("-b_x", "--border_x",action='store', type=float, help="Border added to each plate in the X dimension, in mm.",dest='border_x',default=0)
parser.add_argument("-b_y", "--border_y",action='store', type=float, help="Border added to each plate in the Y dimension, in mm.",dest='border_y',default=0)

parser.add_argument("-p", "--print",action='store', help="Create a printable version. Give page format(A3, A4, A5, Letter, Legal).",dest='print',default='False')
parser.add_argument("-m", "--margin",action='store', type=float, help="Margin to use on the page, in mm.",dest='margin',default=10)


args = parser.parse_args()
map = args.map
logger.info('Map file: '+str(map))
thick = args.thick
logger.info('Slice real-world thickness: '+str(thick))
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
hm_x = args.hm_x
logger.info('Hole mark distance from the X-axis: '+str(hm_x))
hm_y = args.hm_y
logger.info('Hole mark distance from the Y-axis: '+str(hm_y))
hm_w = args.hm_w
logger.info('Width of the holemark: '+str(hm_w))
resample_m = args.resample_m
logger.info('Factor for resampling map: '+str(resample_m))
resample_i = args.resample_i
logger.info('Factor for resampling images: '+str(resample_i))
name = args.name
logger.info('File name prefix: '+str(name))
border_x = args.border_x
logger.info('X border size on each plates: '+str(border_x))
border_y = args.border_y
logger.info('Y border size on each plates: '+str(border_y))
printable = args.print
logger.info('Create a printable version: '+str(printable))
margin = args.margin
logger.info('Margin to use on printable version: '+str(margin))

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

#test that the slice thickness is valid
try:
	assert (thick > 0)
except:
	logger.info('Problem with the thickness value. Quitting.')
	print('Problem with the thickness value. Quitting.')
	sys.exit()
else:
	if not(thick > 0):
		logger.info('Thickness value must be greater than 0. Quitting.')
		print('Thickness value must be greater than 0. Quitting.')
		sys.exit()

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

logger.info('Initial pixel size of the map in each dimension: ('+', '.join([str(voxel_size[x]) for x in range(0,len(voxel_size),1)])+')')

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

	logger.info('Pixel size of the map in each dimension: ('+', '.join([str(voxel_size[x]) for x in range(0,len(voxel_size),1)])+')')
	map_file = resampled_map
else:
	logger.info('Rescale factor must be greater than 0. Quitting.')
	print('Rescale factor must be greater than 0. Quitting.')
	sys.exit()

## create slices along X, Y, or Z
## interpolate along slice axis if box size is not a multiple of the number of slices
map_file = np.array(map_file)
ct_max = np.max(map_file)

if axis.upper()=='X':
	slices = [map_file[x,:,:] for x in range(0,map_file.shape[0],1)]
	scale_x = voxel_size[1]/voxel_size[0]
	scale_y = voxel_size[2]/voxel_size[0]

elif axis.upper()=='Y':
	slices = [map_file[:,x,:] for x in range(0,map_file.shape[1],1)]
	scale_x = voxel_size[0]/voxel_size[1]
	scale_y = voxel_size[2]/voxel_size[1]

elif axis.upper()=='Z':
	slices = [map_file[:,:,x] for x in range(0,map_file.shape[2],1)]
	scale_x = voxel_size[0]/voxel_size[2]
	scale_y = voxel_size[1]/voxel_size[2]


else:
	logger.info('Axis is not x, y, or z. Quitting.')
	print('Axis is not x, y, or z. Quitting.')

logger.info('Scale factor in the X dimension: '+str(scale_x))
logger.info('Scale factor in the Y dimension: '+str(scale_y))

## for each slice calculate the contour lines
half_line_width = round(math.ceil((line_width-1)/2))
logger.info('Half line width: '+str(half_line_width))

##contour simple code
new_slices = []
for slice in slices:
	new_slice=np.zeros(slice.shape,dtype=map_file.dtype)
	to_place = list(set([(x,y) for level in [ct_min+q*(ct_max-ct_min)/ct_num for q in range(0,ct_num,1)] for contour in measure.find_contours(slice, level) for pos in contour for x in range(round(pos[0])-half_line_width,round(pos[0])+half_line_width+1,1) for y in range(round(pos[1])-half_line_width,round(pos[1])+half_line_width+1,1)]))
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


#border the X-dimension
px_per_mm_z = (slice_positions[-1]-slice_positions[0])/(slice_num*thick)
logger.info('Initial pixels per mm (Z dimension): '+str(px_per_mm_z))
px_per_mm_x = px_per_mm_z/scale_x
logger.info('Initial pixels per mm (X dimension): '+str(px_per_mm_x))
px_per_mm_y = px_per_mm_z/scale_y
logger.info('Initial pixels per mm (Y dimension): '+str(px_per_mm_y))

#object physical dimensions
logger.info('Initial object dimensions Z (in mm): '+str(thick*slice_num))
logger.info('Initial object dimensions X (in mm): '+str(slices[0].shape[0]/px_per_mm_x))
logger.info('Initial object dimensions Y (in mm): '+str(slices[0].shape[1]/px_per_mm_y))

#adjusting plate size for borders
if border_x >0:
	border_x = px_per_mm_x*border_x
	if border_x > slices[0].shape[0]:
		logger.info('Adjusting X dimension to include border')
		padding_x = int(round((border_x-slices[0].shape[0])/2))
		padding = np.zeros((padding_x,slices[0].shape[1]))
		slices = [np.concatenate((padding,slice,padding),axis=0) for slice in slices]
	if border_x < slices[0].shape[0]:
		print('WARNING, border smaller than map in the X dimension.')
		logger.info('WARNING, border smaller than map in the X dimension.')

if border_y >0:
	border_y = px_per_mm_y*border_y
	if border_y > slices[0].shape[1]:
		logger.info('Adjusting Y dimension to include border')
		padding_y = int(round((border_y-slices[0].shape[1])/2))
		padding = np.zeros((slices[0].shape[0],padding_y))
		slices = [np.concatenate((padding,slice,padding),axis=1) for slice in slices]
	if border_y < slices[0].shape[1]:
		print('WARNING, border smaller than map in the Y dimension.')
		logger.info('WARNING, border smaller than map in the Y dimension.')

logger.info('Bordered slice dimensions: ('+', '.join([str(x) for x in slices[0].shape])+')')

logger.info('Adding border markings to each plate.')
if ((border_x > 0) and (border_y > 0)):
	logger.info('Drawing centered borders.')
	(center_x,center_y) = (slices[0].shape[0],slices[0].shape[1])
	center_x = int(center_x/2)
	center_y = int(center_y/2)
	left_edge = max(0,int(center_x-border_x/2))
	right_edge = min(int(center_x+border_x/2),slices[0].shape[0])
	bottom_edge = max(0,int(center_y-border_y/2))
	top_edge = min(int(center_y+border_y/2),slices[0].shape[1])
	for slice in slices:
		for x in range(left_edge,right_edge,1):
			for y in range(bottom_edge,bottom_edge+line_width,1):
				slice[x,y]=1
			for y in range(top_edge-line_width,top_edge,1):
				slice[x,y]=1
		for y in range(bottom_edge,top_edge,1):
			for x in range(left_edge,left_edge+line_width,1):
				slice[x,y]=1
			for x in range(right_edge-line_width,right_edge,1):
				slice[x,y]=1

else:
	logger.info('At-least one border size is zero. Edge the images.')
	for slice in slices:
		for x in range(0,slice.shape[0],1):
			for y in range(0,line_width,1):
				slice[x,y]=1
			for y in range(slice.shape[1]-line_width,slice.shape[1],1):
				slice[x,y]=1
		for y in range(0,slice.shape[1],1):
			for x in range(0,line_width,1):
				slice[x,y]=1
			for x in range(slice.shape[0]-line_width,slice.shape[0],1):
				slice[x,y]=1

#add hole-marking square x from the x-edge and y from the y-edge
hm_x = round(hm_x*px_per_mm_x)
hm_y = round(hm_y*px_per_mm_y)
hm_w_x = round(hm_w*px_per_mm_x/2)
hm_w_y = round(hm_w*px_per_mm_y/2)
if hm_w>0:
	logger.info('Adding hole markings '+str(hm_x)+' mm from the X axis and '+str(hm_y)+' mm from the Y axis')
	for slice in slices:
		for x in range(hm_x-hm_w_x,hm_x+hm_w_x,1):
			for y in range(hm_y-line_width,hm_y+line_width,1):
				slice[x,y]=1
		for y in range(hm_y-hm_w_y,hm_y+hm_w_y,1):
			for x in range(hm_x-line_width,hm_x+line_width,1):
				slice[x,y]=1

		for x in range(slice.shape[0]-hm_x-hm_w_x,slice.shape[0]-hm_x+hm_w_x,1):
			for y in range(hm_y-line_width,hm_y+line_width,1):
				slice[x,y]=1
		for y in range(hm_y-hm_w_y,hm_y+hm_w_y,1):
			for x in range(slice.shape[0]-hm_x-line_width,slice.shape[0]-hm_x+line_width,1):
				slice[x,y]=1

		for x in range(hm_x-hm_w_x,hm_x+hm_w_x,1):
			for y in range(slice.shape[1]-hm_y-line_width,slice.shape[1]-hm_y+line_width,1):
				slice[x,y]=1
		for y in range(slice.shape[1]-hm_y-hm_w_y,slice.shape[1]-hm_y+hm_w_y,1):
			for x in range(hm_x-line_width,hm_x+line_width,1):
				slice[x,y]=1

		for x in range(slice.shape[0]-hm_x-hm_w_x,slice.shape[0]-hm_x+hm_w_x,1):
			for y in range(slice.shape[1]-hm_y-line_width,slice.shape[1]-hm_y+line_width,1):
				slice[x,y]=1
		for y in range(slice.shape[1]-hm_y-hm_w_y,slice.shape[1]-hm_y+hm_w_y,1):
			for x in range(slice.shape[0]-hm_x-line_width,slice.shape[0]-hm_x+line_width,1):
				slice[x,y]=1

#calculate dpi based on layer thinkness and pixels in sliced dimension
dpi=px_per_mm_z/0.0393701
logger.info('Initial image DPI: '+str(dpi))

#resampling_factor
dpi = dpi*resample_i
dpi = int(round(dpi))
logger.info('Image DPI after resampling: '+str(dpi))

images = []

## covert to PNG
if not os.path.exists('./slice_images/'):
    os.makedirs('./slice_images/')
for x in range(0,slice_num,1):
	slice = (slices[x]*255).astype(np.uint8)
	image=Image.fromarray(slice, mode='L')
	image = ImageOps.invert(image)
	image2 = image.convert("LA")
	image2.putalpha(Image.fromarray(slice, mode='L'))
	image2 = image2.resize((image2.size[0]*int(round(scale_x*resample_i)),image.size[1]*int(round(scale_y*resample_i))),resample=Image.LANCZOS)
	image2.save('./slice_images/'+name+'_'+str(x)+'.png',dpi=(dpi,dpi))
logger.info('Recorded image dimensions in pixels: ('+', '.join([str(x) for x in image2.size])+')')

image_x = image2.size[0]/(dpi*0.0393701)
image_y = image2.size[1]/(dpi*0.0393701)

#object physical dimensions
logger.info('Final object dimensions X (in mm): '+str(image_x))
logger.info('Final object dimensions Y (in mm): '+str(image_y))
logger.info('Final object dimensions Z (in mm): '+str(thick*slice_num))

#create printable version
if not(printable == 'False'):
	#A4 210x297
	pdf = FPDF('P', 'mm', printable)
	pdf.set_font('Helvetica', '', 12)
	page_x = pdf.w
	page_y = pdf.h
	logger.info('Page size in X dimension: '+str(page_x))
	logger.info('Page size in Y dimension: '+str(page_y))
	page_x=page_x-margin
	page_y=page_y-margin
	logger.info('Assumed page margin: '+str(margin))
	logger.info('Page size in X dimension after margin: '+str(page_x))
	logger.info('Page size in Y dimension after margin: '+str(page_y))
	max_page = max(page_x,page_y)
	if max(image_x,image_y) < max_page:
		orientation_A = math.floor(page_x/image_x)*math.floor(page_y/image_y)
		orientation_B = math.floor(page_x/image_y)*math.floor(page_y/image_x)
		if orientation_B > orientation_A:
			logger.info('Rotating images for optimal positioning')
			for x in range(0,slice_num,1):
			#rotate all the images 90 degrees
#				images = [x.rotate(90) for x in images]
				with Image.open('./slice_images/'+name+'_'+str(x)+'.png') as image:
					image.rotate(90)
					image.save('./slice_images/'+name+'_'+str(x)+'.png',dpi=(dpi,dpi))
		image_x = image2.size[1]/(dpi*0.0393701)
		image_y = image2.size[0]/(dpi*0.0393701)
		images_per_page = math.floor(page_x/image_x)*math.floor(page_y/image_y)
		logger.info('Number of images per page: '+str(images_per_page))
		buffer_x = image_x*(page_x/image_x-math.floor(page_x/image_x))/(math.floor(page_x/image_x)+1)
		buffer_y = image_y*(page_y/image_y-math.floor(page_y/image_y))/(math.floor(page_y/image_y)+1)
		logger.info('Buffer between images in X dimension: '+str(buffer_x))
		logger.info('Buffer between images in Y dimension: '+str(buffer_y))
		coords = [(pos_x*image_x+(pos_x+1)*buffer_x,pos_y*image_y+(pos_y+1)*buffer_y) for pos_x in range(0,math.floor(page_x/image_x),1) for pos_y in range(0,math.floor(page_y/image_y),1)]
		for page in range(0,math.floor(slice_num/images_per_page)+1,1):
			pdf.add_page(orientation='P')
			for num in range(0,images_per_page,1):
				pos_i = page*images_per_page+num
				if pos_i < slice_num:
					pdf.image('./slice_images/'+name+'_'+str(pos_i)+'.png',x=coords[num][0],y=coords[num][1],w=image_x,h=image_y)
					pdf.text(coords[num][0]+image_x/2,coords[num][1]+image_y+buffer_y/4,name+'_'+str(pos_i))
		pdf.output(name+'_slices.pdf','F')
	else:
		logger.info('images are too large to print')
logger.info('Finished successfully!')
