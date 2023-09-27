# MarkField_Generator
A Matlab demo for the paper "Position-Sensing Marker: Mathematic Definition, Generation, and Representation"
************************************************************************************************************
GUIDE
1. "example_10x10_3x3_2.m": generate a 10x10 lattice F with 3x3 G and k=2.
2. "example_50x50_3x3_3.m": generate a 50x50 lattice F with 3x3 G and k=3.
3. "example_80x80_4x4_2.m": generate a 80x80 lattice F with 4x4 G and k=2.
4. "example_torus_4x4_2.m": generate an F on a torus with 4x4 G and k=2 (the Fig.3 in the paper).
************************************************************************************************************
QUICK START
The paper is unavailable yet, thus we offer some comments below for a quick start:
1. F is marker field.
2. G is ID tag shape.
3. k is alphabet size.
4. The index matrix D is an MxN array representing the co-normal isomorphisms between F and G, specifically, the j-th element of the i-th ID tag is the D(i,j)-th element of F (check the "D_3_10.mat" as an example).
************************************************************************************************************
NOTE
1. To generate marker fields on arbitrary solids, we use Blender to mesh the surface (https://www.blender.org).
2. You must build the index matrix D by yourselves. This is easy if F is regular. For more complex F, we use the VF3 algorithm (https://github.com/MiviaLab/vf3lib);
3. Due to some matrix operations in this new generator, we can not utilize the built-in functions of OpenCV to boost the speed, as we did in the previous work "fast-bWFC". Therefore, the C++ version requires a higher programming skill that we do not have. It is worth noting that, in fact, the new generator and fast-bWFC share the same complexity (see Sec.6.2 of the paper), and they should have similar speeds with proper code optimization. Currently, the Matlab version is practical enough, which can also generate large marker fields in minutes, and it is faster than the Matlab version of fast-bWFC.
 
