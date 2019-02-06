//////////// FILES INCLUDED///////////////////////

- ##########.jpeg file       : RGB images recorded
- ##########.pcd file        : point clouds recorded
- ##########_camera.xml file : camera parameters
- ##########_depth.png  file : depth images recorded
- ##########_normal.pcd file : normals extracted
- ##########_mask.jpeg  file : object masks


/////////// M-FILE ///////////////////////////////

- readPcd.m
This m-file is to read provided pcd files by matlab. This files is provided for windows.
This file requires a minor change for other OS to read point cloud data correctly.
The new line command (for windows) in "line 85" is ['\n'], this should be replaced with ['\r \n'] for other OS.

Please contact TAs if it still doesn't work.

Example usage :
pointCloud = readPcd(pcd_img_path);

Note: a returned pointCloud is a matrix of dimension Nx3, but not 3xN as source.mat and target.mat

///////////////////////////////////////////////////

Note: The provided point clouds does not only consist of person (it also contains background). Please apply a distance threshold to remove background (e.g. remove points further than 2meters to the camera).


//////////////////////////////////////////////////

To run experiments of the first exercise run individual_experiments.m to get values presented in the table. The out is stored in folder output0. You can also run grid search of parameters on living room  dataset by running icp_experimets.m and observing folder output1.

To run merging just run merging_scenes_experiments.m and observe folder output2. This run takes a lot of time so we have included generated results in the folder.