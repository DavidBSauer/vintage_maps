rm *sub_map.mrc
rm *masked_map.mrc
rm *contoured_sub_map.mrc
rm *model_maker.log
rm *.png
rm *resampled_map.mrc
rm *.pdf
rm -rf slice_images/

source modelmaker/bin/activate

python3 vintage_map_physical.py -mp map.mrc -mk mask.mrc -st 2.0 -rm 1 -ri 1 -ct_num 3 -ct_min 0.085 -ax X -n test -b_w 1 -b_x 100 -b_y 100  -sl 25 -p A4 -hm_x 10 -hm_y 10 -hm_w 4

deactivate