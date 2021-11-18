import chimera
from chimera import runCommand as RC
#preset publication 1
RC('preset apply pub 1')

#slices
## step1, color #000000000000, opas 1.0, thres 0.005
#open ~/slice_map.mrc
RC('volume #0 style surface level 0.005 color black step 1 transparency 0.0')

#plates
## step1, color #00000000e275, opas 0.015, thres 0.99
#open ~/plate_map.mrc
RC('volume #1 style surface level 0.99 color #00000000e275 step 1 transparency 0.99')

#model
#open ~/6wu1.pdb
RC('focus #0')
RC('move cofr mod #2')
