name='modelmaker'

virtualenv $name

virtualenv -p /usr/bin/python3 $name

source $name/bin/activate

pip3 install mrcfile scikit-image Pillow numpy matplotlib fpdf

deactivate

#to remove virtualenv "sudo rm -rf name"
