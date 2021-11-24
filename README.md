# Vintage map maker
From a .mrc map, create a virtual vintage-style representation or templates for a physical version. Tested on MacOS, but should work on Linux also (perhaps with minor tweaks).

## Virtual vintage-style maps
1. Install the necessary packages. 
This assumes you have virtualenv installed.

```
./setup.bash
```

2. Edit the map generation parameters.
You can run the script from command line, but I found it easier to use the script virtual_runner.bash. 

See the parameters using
```
python3 vintage_map_virtual.py -h
```
Edit the virtual_runner.bash as you see fit.

Note: Resampling the map takes time and the data requirements grow as the cube of the factor provided. I typically use 1-2 for testing, and 4-5 for the final run.

3. Run the virtual map maker. 
After editing the virtual_runner.bash script, run it using
```
./virtual_runner.bash
```

4. Chimera rendering
The virtual_runner.bash should automatically launch Chimera. You will need to then open the chimera_script.py, which should then re-draw the maps in a vintage style.

## Physical templates
1. Install the necessary packages. 
This assumes you have virtualenv installed.

```
./setup.bash
```

2. Edit the map generation parameters.
See the parameters using
```
python3 physical_map_virtual.py -h
```
Edit the physical_runner.bash as you see fit.

Note: I would suggest using the virtual map maker to optimize everything before running creating a physical template.  

Note: The critical parameter here is the physical slice thickness. This then determines all other real-world dimensions.

Note: Resampling the map takes time and the data requirements grow as the cube of the factor provided. I typically use 1-2 for testing, and 4-5 for the final run.

3. Run the physical map maker. 
After editing the physical_runner.bash script, run it using
```
./physical_runner.bash
