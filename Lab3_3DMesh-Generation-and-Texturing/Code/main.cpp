
#include <stdlib.h>
#include <boost/format.hpp>
#include <iostream>

#include <pcl/common/transforms.h>
#include <pcl/features/integral_image_normal.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/filters/passthrough.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/search/kdtree.h>
#include <pcl/surface/marching_cubes.h>
#include <pcl/surface/marching_cubes_hoppe.h>
#include <pcl/surface/marching_cubes_rbf.h>
#include <pcl/surface/poisson.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/surface/impl/texture_mapping.hpp>
#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_lib_io.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <opencv2/core/eigen.hpp>
#include <opencv2/core/mat.hpp>
#include <opencv2/opencv.hpp>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <pcl/surface/mls.h>
#include <pcl/surface/vtk_smoothing/vtk_mesh_smoothing_windowed_sinc.h>
#include <pcl/surface/vtk_smoothing/vtk_utils.h>
#include <pcl/visualization/common/actor_map.h>
#include <vtkFillHolesFilter.h>

#include "Frame3D/Frame3D.h"

#include <algorithm>  // std::max
#include <sys/stat.h> // mkdir
#include <sstream>





using namespace std;

using Pcl = pcl::PointCloud<pcl::PointXYZRGBNormal>;


boost::shared_ptr<pcl::visualization::PCLVisualizer> rgbVis(pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr cloud);

pcl::PointCloud<pcl::PointXYZ>::Ptr mat2IntegralPointCloud(const cv::Mat &depth_mat, const float focal_length,
                                                           const float max_depth);

pcl::PointCloud<pcl::PointXYZRGB>::Ptr mat2IntegralPointCloudRGB(const cv::Mat &depth_mat, const cv::Mat &rgb_image,
                                                                 const float focal_length, const float max_depth);

pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr computeNormalsWithTexture(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud) {
  // This function computes normals given a point cloud
  // !! Please note that you should remove NaN values from the pointcloud after
  // computing the surface normals.
  pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr cloud_normals(new pcl::PointCloud<pcl::PointXYZRGBNormal>);  // Output
                                                                                                            // datasets
  pcl::IntegralImageNormalEstimation<pcl::PointXYZRGB, pcl::PointXYZRGBNormal> ne;
  ne.setNormalEstimationMethod(ne.AVERAGE_3D_GRADIENT);
  ne.setMaxDepthChangeFactor(0.02f);
  ne.setNormalSmoothingSize(10.0f);
  ne.setInputCloud(cloud);
  ne.compute(*cloud_normals);
  pcl::copyPointCloud(*cloud, *cloud_normals);
  return cloud_normals;
}

pcl::PointCloud<pcl::PointNormal>::Ptr computeNormals(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud) {
  // This function computes normals given a point cloud
  // !! Please note that you should remove NaN values from the pointcloud after computing the surface normals.
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_normals(new pcl::PointCloud<pcl::PointNormal>);  // Output datasets
  pcl::IntegralImageNormalEstimation<pcl::PointXYZ, pcl::PointNormal> ne;
  ne.setNormalEstimationMethod(ne.AVERAGE_3D_GRADIENT);
  ne.setMaxDepthChangeFactor(0.02f);
  ne.setNormalSmoothingSize(10.0f);
  ne.setInputCloud(cloud);
  ne.compute(*cloud_normals);
  pcl::copyPointCloud(*cloud, *cloud_normals);
  return cloud_normals;
}

pcl::PointCloud<pcl::PointXYZRGB>::Ptr transformPointCloud(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud,
                                                           const Eigen::Matrix4f &transform) {
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr transformed_cloud(new pcl::PointCloud<pcl::PointXYZRGB>());
  pcl::transformPointCloud(*cloud, *transformed_cloud, transform);
  return transformed_cloud;
}

template <class T>
typename pcl::PointCloud<T>::Ptr transformPointCloudNormals(typename pcl::PointCloud<T>::Ptr cloud,
                                                            const Eigen::Matrix4f &transform) {
  typename pcl::PointCloud<T>::Ptr transformed_cloud(new typename pcl::PointCloud<T>());
  pcl::transformPointCloudWithNormals(*cloud, *transformed_cloud, transform);
  return transformed_cloud;
}


pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr mergingPointClouds(Frame3D frames[], float depth_threshold) {
  Pcl::Ptr modelCloud(new Pcl());

  for (int i = 0; i < 8; i++) {
    std::cout << boost::format("Merging frame %d") % i << std::endl;

    Frame3D frame = frames[i];
    cv::Mat depthImage = frame.depth_image_;
    double focalLength = frame.focal_length_;
    const Eigen::Matrix4f cameraPose = frame.getEigenTransform();

    // TODO(Student): Merge the i-th frame using predicted camera pose
    // to the global point cloud. ~ 20 lines.

    pcl::PointCloud<pcl::PointXYZ>::Ptr point_cloud = mat2IntegralPointCloud(depthImage, focalLength, depth_threshold);

    pcl::PointCloud<pcl::PointNormal>::Ptr point_cloud_with_normals = computeNormals(point_cloud);

    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr new_pcl(new pcl::PointCloud<pcl::PointXYZRGBNormal>());

    new_pcl->reserve(point_cloud_with_normals->size());
    for (auto &&pt : *point_cloud_with_normals) {
      if (!isnan(pt.normal_x) && !isnan(pt.normal_y) && !isnan(pt.normal_z)) {
        pcl::PointXYZRGBNormal new_pt;
        new_pt.x = pt.x;
        new_pt.y = pt.y;
        new_pt.z = pt.z;
        new_pt.normal_x = pt.normal_x;
        new_pt.normal_y = pt.normal_y;
        new_pt.normal_z = pt.normal_z;
        new_pt.b = new_pt.r = new_pt.g = 128;

        new_pcl->push_back(new_pt);
      }
    }

    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr transformed_point_cloud = transformPointCloudNormals<pcl::PointXYZRGBNormal>(new_pcl, cameraPose);

    *modelCloud += (*transformed_point_cloud);
  }

  return modelCloud;
}

pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr mergingPointCloudsWithTexture(Frame3D frames[], float depth_threshold) {
  pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr modelCloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);

  for (int i = 0; i < 8; i++) {
    std::cout << boost::format("Merging frame %d") % i << std::endl;

    Frame3D frame = frames[i];
    cv::Mat depthImage = frame.depth_image_;
    double focalLength = frame.focal_length_;
    const Eigen::Matrix4f cameraPose = frame.getEigenTransform();

    // TODO(Student): The same as mergingPointClouds but now with texturing. ~ 50 lines.

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_cloud = mat2IntegralPointCloudRGB(depthImage, frame.rgb_image_, focalLength, depth_threshold);


    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr point_cloud_with_normals = computeNormalsWithTexture(point_cloud);

    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr new_pcl(new pcl::PointCloud<pcl::PointXYZRGBNormal>());

    new_pcl->reserve(point_cloud_with_normals->size());
    for (auto &&pt : *point_cloud_with_normals) {
      if (!isnan(pt.normal_x) && !isnan(pt.normal_y) && !isnan(pt.normal_z)) {
        new_pcl->push_back(pt);
      }
    }

    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr transformed_point_cloud =
        transformPointCloudNormals<pcl::PointXYZRGBNormal>(new_pcl, cameraPose);

    *modelCloud += (*transformed_point_cloud);
  }
  return modelCloud;
}

// Different methods of constructing mesh
enum CreateMeshMethod { PoissonSurfaceReconstruction = 0, MarchingCubes = 1 };

// Create mesh from point cloud using one of above methods
pcl::PolygonMesh createMesh(pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr pointCloud, CreateMeshMethod method,
                            int method_specific_param, float pass, float hole_size) {
  std::cout << "Creating meshes" << std::endl;

  // The variable for the constructed mesh
  pcl::PolygonMesh triangles;

  pcl::Poisson<pcl::PointXYZRGBNormal> psr;
  pcl::MarchingCubesHoppe<pcl::PointXYZRGBNormal> mcr;

  pcl::PolygonMesh::Ptr unsmoothed_triangles(new pcl::PolygonMesh());

  int depth, resolution;

  switch (method) {
    case PoissonSurfaceReconstruction:
      // TODO(Student): Call Poisson Surface Reconstruction. ~ 5 lines.

      /* Around 8 */
      depth = method_specific_param;
      /* It's the maximum depth of the octree used to build the indicator function for the surface.
       * It is directly related to the maximum resolution allowed.
      /* Higher values enable a more fine-graned reconstruction, but more sensitive to noise
      /* Low values produce smoother surfaces but lose many details (e.g. nose, eyes, ...)
       */
      psr.setDepth(depth);

      /* the ratio between the diameter of the cube used for reconstruction and the diameter of the samples' bounding
       * cube */
      /*
       * 1 seems too small as part of one arm is not reconstructed properly
       * too large values deform the lower-part of the body and add surfaces outside the body
       * the default one and something aroud 1.25 works well
       */
      //~ psr.setScale(1.25);

      /* minimum number of sample points that should fall within an octree node as the octree construction is adapted to
       * sampling density */
      /*
       * The documentation suggests:
       * 		"For noise-free samples, small values in the range [1.0 - 5.0] can be used.
       * 		For more noisy samples, larger values in the range [15.0 - 20.0] may be needed to provide
       * 		a smoother, noise-reduced, reconstruction."
       * In our case the point cloud doesn't seem very noisy, and we saw law values work pretty well,
       * with little difference with the default one.
       * Increasing more the value, results in gradually smoother results
       */
      //~ psr.setSamplesPerNode(3);

      psr.setOutputPolygons(false);
      psr.setInputCloud(pointCloud);
      psr.reconstruct(*unsmoothed_triangles);

      break;
    case MarchingCubes:
      // TODO(Student): Call Marching Cubes Surface Reconstruction. ~ 5 lines.

      /* Something between 50 - 100*/
      resolution = method_specific_param;

      /* It's the number of cubes the space is divided in (number of ticks in each of the axes)
       * Higher values mean more cubes and, therefore, a more detailed reconstruction. However, using many cubes
       * the results are extremely more sensitive to noise. Moreover, using too many cubes requires too much memory.
       * Indeed, the default value results in the RAM to be completely filled before terminating.
       * Low values produces very bad results, preserving only the general shape of the surface
       */
      mcr.setGridResolution(resolution, resolution, resolution);

      /* They don't seem to improve the results if changed w.r.t. the default ones*/

      /* Is a threshold (in [0, 1]) for choosing the position of a point w.r.t the surface
       * Too high values (> 0.1) result in an empty mesh because no/few cubes with points on both side of the surface
       * are found. With 0.1 the mesh seems to be "reverted" w.r.t. the cloud (it builds 2 meshes with the shapes of the
       * space on the left and on the right of the body) With 0.0 (should be the default one) the results make more
       * sense
       */
      //~ mcr.setIsoLevel(0.1);

      /* enlarge the resulting mesh, the default value matches the size of the point cloud better*/
      //~ mcr.setPercentageExtendGrid(0.7);

      mcr.setInputCloud(pointCloud);
      mcr.reconstruct(*unsmoothed_triangles);

      break;
  }

  if (pass > 0) {
    /* It seems that it doesn't affect significantly the results.
     * Actually, there is no visual relevant difference.
     */

    printf("Smoothing\n");
    pcl::MeshSmoothingWindowedSincVTK vtk;
    vtk.setPassBand(pass);
    vtk.setInputMesh(unsmoothed_triangles);
    vtk.setNumIter(30);
    vtk.setFeatureEdgeSmoothing(true);
    vtk.setFeatureAngle(M_PI / 10);
    vtk.setBoundarySmoothing(true);
    vtk.setEdgeAngle(M_PI / 10);
    vtk.process(triangles);
  } else {
    triangles = *unsmoothed_triangles;
  }

  if (hole_size > 0.) {
    printf("Filtering\n");
    vtkSmartPointer<vtkPolyData> input;
    pcl::VTKUtils::mesh2vtk(triangles, input);

    vtkSmartPointer<vtkFillHolesFilter> fillHolesFilter = vtkSmartPointer<vtkFillHolesFilter>::New();

    fillHolesFilter->SetInputData(input);

    /* It looks like any value greater than or equal to 1 works, and there is no visual difference using higher values*/
    fillHolesFilter->SetHoleSize(hole_size);
    fillHolesFilter->Update();

    vtkSmartPointer<vtkPolyData> polyData = fillHolesFilter->GetOutput();

    pcl::VTKUtils::vtk2mesh(polyData, triangles);
  }

  return triangles;
}


void textureMesh(pcl::PolygonMesh *mesh, Frame3D frames[]) {
	
	pcl::PointCloud<pcl::PointXYZ> mesh_cloud;
	pcl::fromPCLPointCloud2(mesh->cloud, mesh_cloud);

	pcl::TextureMapping<pcl::PointXYZ> texture_map;
	pcl::PointCloud<pcl::PointXYZRGB> color_mesh_cloud;
	pcl::copyPointCloud(mesh_cloud, color_mesh_cloud);
	
	Eigen::Vector2f uv_coords;
	pcl::PointXYZRGB color_point;
	pcl::PointXYZ point;
	cv::Vec3b pixel;


	int colors[color_mesh_cloud.points.size()][4];
	for (int i=0; i<color_mesh_cloud.points.size(); i++){
		colors[i][0] = colors[i][1] = colors[i][2] = colors[i][3] = 0;
	}

	const double RESOLUTION = 0.01;
	
	for (int f = 0; f < 8; f++) {

		printf("Coloring using Frame %d\n", f);
	
		Frame3D frame = frames[f];

		// resize image to depth-map size (i.e. the range of coordinates in the point cloud)
		cv::Mat rgb_image_;
		cv::resize(frame.rgb_image_, rgb_image_, frame.depth_image_.size());
		
		
		pcl::TextureMapping<pcl::PointXYZ>::Camera camera;
		camera.focal_length = frame.focal_length_;
		camera.width = rgb_image_.cols;
		camera.height = rgb_image_.rows;


		// transform mesh points into camera frame
		auto cam_pose = frame.getEigenTransform();
		pcl::PointCloud<pcl::PointXYZ>::Ptr cam_cloud(new pcl::PointCloud<pcl::PointXYZ>());
		pcl::transformPointCloud(mesh_cloud, *cam_cloud, Eigen::Affine3f(cam_pose).inverse().matrix());
	
		pcl::TextureMapping<pcl::PointXYZ>::Octree::Ptr tree(new pcl::TextureMapping<pcl::PointXYZ>::Octree(RESOLUTION));
		tree->setInputCloud(cam_cloud);
		tree->addPointsFromInputCloud();
	
		// loop over all points (vertices) in the mesh
		int verbose = color_mesh_cloud.points.size() / 10;
		for (int vertex_id = 0; vertex_id < color_mesh_cloud.points.size(); vertex_id++) {
		
			if (vertex_id % verbose == 0)
				printf("\t%d %% done\n", 100*vertex_id / int(color_mesh_cloud.points.size()));
		
			point = cam_cloud->points.at(vertex_id);
			
			// check if the point is visible by the camera (it is not occluded)
			if (texture_map.isPointOccluded(point, tree))
				continue;
		
			// check if the point is inside the image
			if (!texture_map.getPointUVCoordinates(point, camera, uv_coords))
				continue;
		
			// retrieve the corresponding pixel
			int u = (int)(uv_coords[0] * camera.width);
			int v = (int)(camera.height - (uv_coords[1] * camera.height));
			
			pixel = rgb_image_.at<cv::Vec3b>(v, u);


			/* In order to average colors from multiple images
			/* we keep track of the sum of the colors and the number
			 * of frames in which the point was visible 
			 */
			colors[vertex_id][3]++;
			colors[vertex_id][0] += pixel[2];
			colors[vertex_id][1] += pixel[1];
			colors[vertex_id][2] += pixel[0];
			
		}
	}
	
	/* 
	 * Compute the average of the pixels from the view in which each point appears
	 */
	for (int i=0; i<color_mesh_cloud.points.size(); i++){
		if (colors[i][3] > 0){
			color_mesh_cloud.points.at(i).r = (uint8_t)(colors[i][0] / colors[i][3]);
			color_mesh_cloud.points.at(i).g = (uint8_t)(colors[i][1] / colors[i][3]);
			color_mesh_cloud.points.at(i).b = (uint8_t)(colors[i][2] / colors[i][3]);
		}
	}
	
	// put the colored points back into the mesh
	pcl::toPCLPointCloud2(color_mesh_cloud, mesh->cloud);
}



int main(int argc, char *argv[]) {
    if (argc < 4) {
      std::cout << "./final [3DFRAMES PATH] [RECONSTRUCTION MODE] [TEXTURE_MODE]" << std::endl;
      return 0;
    }

    srand(time(NULL));
    

    const CreateMeshMethod reconMode = static_cast<CreateMeshMethod>(std::stoi(argv[2]));

    // Loading 3D frames
    Frame3D frames[8];
    for (int i = 0; i < 8; ++i) {
      frames[i].load(boost::str(boost::format("%s/%05d.3df") % argv[1] % i));
    }

    float depth_threshold = 1.0;
    int param;

    switch (reconMode) {
      case PoissonSurfaceReconstruction:
        param = 8;
        break;
      case MarchingCubes:
        param = 100;
        break;
    }

    float hole_size = 100.0;
    float pass = -1;

    bool visualize_points = true;
    bool visualize_mesh = true;

    if (argc > 4) {
      depth_threshold = atof(argv[4]);
    }

    if (argc > 5) {
      param = atoi(argv[5]);
    }

    if (argc > 6) {
      pass = atof(argv[6]);
    }

    if (argc > 7) {
      hole_size = atof(argv[7]);
    }

    if (argc > 8) {
      visualize_points = atoi(argv[8]);
    }
    
    if (argc > 9) {
      visualize_mesh = atoi(argv[9]);
    }

	// Building the Experiment Directory
	std::string root_experiments = "./experiments/";
	mkdir(root_experiments.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
	
	std::stringstream dir_name;
	dir_name << root_experiments << "exp_" << (int)reconMode << "_" << argv[3][0] << "_" << depth_threshold << "_" << param << "_" << pass << "_" << hole_size << "/";
    
    mkdir(dir_name.str().c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
	
	cout << endl << endl << "Experiment: " << dir_name.str() << endl;
	
    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr texturedCloud;
    pcl::PolygonMesh triangles;

    if (argv[3][0] == 't') {
      // SECTION 4: Coloring 3D Model
      // Create one point cloud by merging all frames with texture using
      // the rgb images from the frames

      texturedCloud = mergingPointCloudsWithTexture(frames, depth_threshold);

      // Create a mesh from the textured cloud using a reconstruction method,
      // Poisson Surface or Marching Cubes
	  triangles = createMesh(texturedCloud, reconMode, param, pass, hole_size);

	  // Add texture to the mesh
	  textureMesh(&triangles, frames);

    } else {
      // SECTION 3: 3D Meshing & Watertighting

      // Create one point cloud by merging all frames with texture using
      // the rgb images from the frames
      texturedCloud = mergingPointClouds(frames, depth_threshold);
      
      // Create a mesh from the textured cloud using a reconstruction method,
      // Poisson Surface or Marching Cubes
	  triangles = createMesh(texturedCloud, reconMode, param, pass, hole_size);
    }
	
	
	// Storing results to files
	std::stringstream cloud_file;
	cloud_file << dir_name.str() << "cloud.pcd";
    
	pcl::io::savePCDFile(cloud_file.str().c_str(), *texturedCloud, true);
	
	std::stringstream mesh_file;
	mesh_file << dir_name.str() << "mesh.pcd";
	
	pcl::io::savePolygonFile(mesh_file.str().c_str(), triangles);
        

    if (visualize_points || visualize_mesh) {
		// Sample code for visualization.

		// Show viewer
		std::cout << "Finished texturing" << std::endl;
		boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
		
		if (visualize_points){
			// Add colored point cloud to viewer, because it does not support colored
			// meshes
			pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal> rgb(texturedCloud);
			viewer->addPointCloud<pcl::PointXYZRGBNormal>(texturedCloud, rgb, "cloud");
			viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "cloud");
		}
		
		if (visualize_mesh) {
			// Add mesh
			viewer->addPolygonMesh(triangles, "meshes", 0);
		}
		
		viewer->setBackgroundColor(1, 1, 1);
		viewer->addCoordinateSystem(1.0);
		viewer->initCameraParameters();

		// Keep viewer open
		while (!viewer->wasStopped()) {
		  viewer->spinOnce(100);
		  boost::this_thread::sleep(boost::posix_time::microseconds(100000));
		}
	}
	
    return 0;
  }

  pcl::PointCloud<pcl::PointXYZRGB>::Ptr mat2IntegralPointCloudRGB(const cv::Mat &depth_mat, const cv::Mat &rgb_image,
                                                                   const float focal_length, const float max_depth) {
    // This function converts a depth image to a point cloud
    assert(depth_mat.type() == CV_16U);

    //~ cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE/10 );// Create a window for display.
    //~ cv::imshow( "Display window", rgb_image);
    //~ cv::waitKey(20);

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_cloud(new pcl::PointCloud<pcl::PointXYZRGB>());
    const int half_width = depth_mat.cols / 2;
    const int half_height = depth_mat.rows / 2;
    const float inv_focal_length = 1.0 / focal_length;

    double downscale_x = (double)rgb_image.cols / (double)depth_mat.cols;
    double downscale_y = (double)rgb_image.rows / (double)depth_mat.rows;

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr rgb_cloud(new pcl::PointCloud<pcl::PointXYZRGB>());
    for (int i = 0; i < rgb_image.rows; i++) {
      for (int j = 0; j < rgb_image.cols; j++) {
        cv::Vec3b pixel = rgb_image.at<cv::Vec3b>(i, j);
        rgb_cloud->points.emplace_back(pixel[2], pixel[1], pixel[0]);
        // expects that y axis start at [0,0] in top left position and goes down in depth image
        // x axis  goes to the right
        rgb_cloud->points.back().x = (double)j / downscale_x;
        rgb_cloud->points.back().y = (double)i / downscale_y;
        rgb_cloud->points.back().z = 0;
      }
    }

    pcl::KdTreeFLANN<pcl::PointXYZRGB> kdtree;
    kdtree.setInputCloud(rgb_cloud);

    int K = 1;
    float radius = 2;

    std::vector<int> idxs_radius;
    std::vector<float> distances_radius;

    //~ cv::Mat image(depth_mat.rows, depth_mat.cols, rgb_image.type());

    point_cloud->points.reserve(depth_mat.rows * depth_mat.cols);
    for (int y = 0; y < depth_mat.rows; y++) {
      for (int x = 0; x < depth_mat.cols; x++) {
        float z = depth_mat.at<ushort>(cv::Point(x, y)) * 0.001;
        if (z < max_depth && z > 0) {
          pcl::PointXYZRGB search_pt(0, 0, 0);
          search_pt.x = x;
          search_pt.y = y;

          pcl::PointXYZRGB new_pt;

          if (kdtree.radiusSearch(search_pt, radius, idxs_radius, distances_radius) > 0) {
            pcl::PointXYZRGB pt = rgb_cloud->points[idxs_radius[0]];

            new_pt = pcl::PointXYZRGB(pt);

          } else {
            new_pt = pcl::PointXYZRGB(0, 0, 0);
          }

          new_pt.x = static_cast<float>(x - half_width) * z * inv_focal_length;
          new_pt.y = static_cast<float>(y - half_height) * z * inv_focal_length;
          new_pt.z = z;

          point_cloud->points.push_back(new_pt);

        } else {
          point_cloud->points.emplace_back(0, 0, 0);
          point_cloud->points.back().x = x;
          point_cloud->points.back().y = y;
          point_cloud->points.back().z = NAN;
        }

        //~ cv::Vec3b color;
        //~ color[0] = point_cloud->points.back().r;
        //~ color[1] = point_cloud->points.back().g;
        //~ color[2] = point_cloud->points.back().b;

        //~ image.at<cv::Vec3b>(cv::Point(x,y)) = color;
      }
    }

    point_cloud->width = depth_mat.cols;
    point_cloud->height = depth_mat.rows;

    //~ cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE/10 );// Create a window for display.
    //~ cv::imshow( "Display window", image);
    //~ cv::waitKey(20);

    return point_cloud;
  }

  pcl::PointCloud<pcl::PointXYZ>::Ptr mat2IntegralPointCloud(const cv::Mat &depth_mat, const float focal_length,
                                                             const float max_depth) {
    // This function converts a depth image to a point cloud
    assert(depth_mat.type() == CV_16U);
    pcl::PointCloud<pcl::PointXYZ>::Ptr point_cloud(new pcl::PointCloud<pcl::PointXYZ>());
    const int half_width = depth_mat.cols / 2;
    const int half_height = depth_mat.rows / 2;
    const float inv_focal_length = 1.0 / focal_length;
    point_cloud->points.reserve(depth_mat.rows * depth_mat.cols);
    for (int y = 0; y < depth_mat.rows; y++) {
      for (int x = 0; x < depth_mat.cols; x++) {
        float z = depth_mat.at<ushort>(cv::Point(x, y)) * 0.001;
        if (z < max_depth && z > 0) {
          point_cloud->points.emplace_back(static_cast<float>(x - half_width) * z * inv_focal_length,
                                           static_cast<float>(y - half_height) * z * inv_focal_length, z);
        } else {
          point_cloud->points.emplace_back(x, y, NAN);
        }
      }
    }
    point_cloud->width = depth_mat.cols;
    point_cloud->height = depth_mat.rows;
    return point_cloud;
  }

  boost::shared_ptr<pcl::visualization::PCLVisualizer> rgbVis(pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr cloud) {
    // --------------------------------------------
    // -----Open 3D viewer and add point cloud-----
    // --------------------------------------------
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
    viewer->setBackgroundColor(1, 1, 1);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(cloud);
    viewer->addPointCloud<pcl::PointXYZRGB>(cloud, rgb, "sample cloud");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "sample cloud");
    viewer->addCoordinateSystem(1.0);
    viewer->initCameraParameters();
    return (viewer);
  }

  
