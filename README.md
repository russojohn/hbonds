# hbonds
Compute hydrogen bonds for SW type potentials

The executable is called 'bonds_map_topologic'. It needs a configurations in the following format

first line: [time stamp] [number of atoms] [box_x] [box_y] [box_z]
other lines: [x] [y] [z] for each atom


In the directory I included a sample "pos.dat" file.

To run the program you then need to specify the SW lambda parameter (which for mW water is 23.15). So you can run the program as

./bonds_map_topologic pos.dat 23.15

It will give you the list of all neighbours for each particle.
