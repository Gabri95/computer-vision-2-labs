# Assignment 3: 3D Mesh Generation and Texturing


In this assignment we are going to work at 3D Mesh Reconstruction and Texturing from 2D images
and depth maps. For this purpose, we are going to use C++ libraries provided in PointCloudLibrary
(PCL). The first step towards 3D reconstruction is to generate a 3D Point cloud of the scene from
depth images, as we have seen in the previous assignments. Once the point cloud is generated
we can choose a surface reconstruction algorithm to interpolate points by generating a 3D surface
connecting them. Examples of such algorithms are Poisson Reconstruction and Marching Cube.
Some post-processing techniques are then required, including smoothing surfaces and closing holes
in volumes by connecting surfaces (watertighting). Once the 3D Mesh is generated, we will proceed
with Texturing. This final step aims to map color information from 2D images to RGB polygons or
points in the 3D reconstruction.

---



Let's name the project's main directory "ROOT"

Enter in the folder "ROOT/Code"

Run "./builder.sh [-c] [params]" to compile the source code and link the libraries

    -c          -> if set, it re-builds also the whole project (including running cmake for the dependencies etc..)
                                   (N.B.: you might want to change the compiler "CC" at line 8)
    
    [params]    -> if any other parameter is set, the executable "final" is run and these parameters are passed to it
    
    
In order to just run the program:

./final ../3dframes RECONSTRUCTION_MODE TEXTURE_MODE [depth_threshold [param [bandpass_filter [hole_size [plot_point_cloud [plot_mesh]]]]]]

where the parameters in "[]" are optional (N.B.: order matters)


In order to re-run the experiments:

./experiments.sh


