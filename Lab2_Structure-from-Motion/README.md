# Assignment 2: Structure-from-Motion

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

To run the experiments, execute "experiments_*.m" scripts in the respective subdirectories of ./src
You can find already produced results under ./output


Please add the following content to the project before running the code:

VLFeat library:
  ./src/vlfeat

Dataset
  ./Data/House
