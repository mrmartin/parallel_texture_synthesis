# parallel_texture_synthesis
Parallel Controllable Texture Synthesis
by Martin

Input texture:
![input](https://github.com/mrmartin/parallel_texture_synthesis/raw/master/input_texture.png)

Output:
![output](https://github.com/mrmartin/parallel_texture_synthesis/raw/master/output_texture.png)

Implements the 2005 paper by Lefebvre and Hoppe (Lefebvre, Sylvain, and Hugues Hoppe. "Parallel controllable texture synthesis." ACM Transactions on Graphics (ToG). Vol. 24. No. 3. ACM, 2005.) in Matlab. 

Includes the Gaussian Stack approach of section 3.1, including the bells and whistles of coherent synthesis and k-coherent search.

Includes some basic CPU parallelism

This project also contains new research in terrain synthesis with realistic watershed basins, performed with synthesize_terrain.m. This includes a waterflow preprocessing step, the candidate neighbour preprocessing, and synthesis with additional flow constraints.
