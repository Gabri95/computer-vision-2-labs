# Labs for Computer Vision 2 course, MSc AI @ UvA 2018/2019.


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  
Solutions and implementation from [Davide Belli](https://github.com/davide-belli) and [Gabriele Cesa](https://github.com/Gabri95).

---

## Assignment 2: Structure-from-Motion

In this assignment we are going to implement Structure-from-Motion algorithm with the goal
to reconstruct a 3D-structure from a set of images picturing the same object. At first we are
going to implement Eight-point algorithm to compute the Fundamental matrix describing the 3D
transformation from a view to the next one. To construct this algorithm, we start from a basic version
of it, improving it successively with Normalization of data and applying RANSAC to find the best
transformation. Afterwards, we are going to iteratively match pairs of views in the dataset in order
to create a single chain structure connecting all the views. The match graph created in this way
is represented as a sparse point-view matrix. This matrix will be compared with the ground-truth
dense matrix representing the transformation of the whole point cloud. Finally, we will use these
matrices for affine Structure-from-Motion. To solve the problem with the sparse matrix, we are
going to work with dense blocks in the matrix. By computing once again SVD composition, the
matrices representing Structure and Motion can be derived. Some improvements among which affine
ambiguity removal will be discussed and implemented to conclude this assignment.

---

## Assignment 3: 3D Mesh Generation and Texturing


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




## Copyright

Copyright © 2019 Davide Belli, Gabriele Cesa.

<p align=“justify”>
This project is distributed under the <a href="LICENSE">MIT license</a>.  
Please follow the <a href="http://student.uva.nl/en/content/az/plagiarism-and-fraud/plagiarism-and-fraud.html">UvA regulations governing Fraud and Plagiarism</a> in case you are a student.
</p>

