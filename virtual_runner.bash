rm *sub_map.mrc
rm *masked_map.mrc
rm *contoured_sub_map.mrc
rm *model_maker.log
rm *.png
rm *resampled_map.mrc
rm *slice_map.mrc
rm *plate_map.mrc

source modelmaker/bin/activate

python3 vintage_map_virtual.py -mp LaINDY_gaussian_2.mrc -rm 4 -ct_num 1 -ct_min 1.64 -ax X -n LaINDY -sl 15

deactivate