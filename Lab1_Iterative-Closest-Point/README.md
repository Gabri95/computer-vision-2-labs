# Assignment 1: Iterative Closest Point

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
