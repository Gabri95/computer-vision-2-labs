# Labs for Computer Vision 2 course, MSc AI @ UvA 2018/2019.


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  
Solutions and implementation from [Davide Belli](https://github.com/davide-belli), [Gabriele Cesa](https://github.com/Gabri95) and [Lukáš Jelínek](https://github.com/Lukx19).

---

## Assignment 1: Iterative Closest Point

In this assignment, we are going to work with Iterative Closest Point algorithms, particularly focusing
on possible improvements, pitfalls, applications and theoretical concepts. Initially, we implement
a basic version of ICP, which is ran both on dummy points-cloud data and on actual real-world
data. Then, we try different ways to improve efficiency and accuracy of ICP, including some ideas
proposed in [1] and [2]. Next, we employ ICP algorithm to estimate the camera poses in different
frames. Once the camera pose is detected, the points-cloud from those frames can be merged. This
technique can be used for practical applications such as generating 3D models out of a set of pictures
from the same environment/object. In our case, we try to reconstruct the 3D model of the subject
pictured from different angles in the real-world dataset. In addition, we are experimenting with
running ICP on non-consecutive frames and with using the merged points-cloud up to the current
frame as target in camera pose estimation. Possible changes are visualized, commenting the reasons
behind improvements and pitfalls. Based on the observations and knowledge gained so far, we try
to explain how drawbacks in ICP algorithm can be overcome, including additional improvements
discussed in the additional documentation or proposed in relative literatures.

[1] S. Rusinkiewicz and M. Levoy. Efficient variants of the ICP algorithm. Third International Conference on 3D Digital Imaging and Modeling (3DIM), 2001.
[2] Z. Zhang. Iterative Point Matching for Registration of Free-Form Curves and Surfaces. International Journal of Computer Vision, 13:2, 119-152, 1994.

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

Copyright © 2019 Davide Belli, Gabriele Cesa and Lukáš Jelínek.

<p align=“justify”>
This project is distributed under the <a href="LICENSE">MIT license</a>.  
Please follow the <a href="http://student.uva.nl/en/content/az/plagiarism-and-fraud/plagiarism-and-fraud.html">UvA regulations governing Fraud and Plagiarism</a> in case you are a student.
</p>

