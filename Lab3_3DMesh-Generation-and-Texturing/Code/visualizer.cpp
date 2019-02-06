
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


int main(int argc, char *argv[]) {
    if (argc < 2) {
      std::cout << "./final [file]" << std::endl;
      return 0;
    }

    srand(time(NULL));
    
    pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr texturedCloud;
    
    pcl::PolygonMesh triangles;
    
    cout << pcl::io::loadPolygonFile(argv[1], triangles) << endl;

	        

	// Sample code for visualization.

	// Show viewer
	std::cout << "Finished texturing" << std::endl;
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
	
	// Add colored point cloud to viewer, because it does not support colored
	// meshes
	//~ pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBNormal> rgb(texturedCloud);
	//~ viewer->addPointCloud<pcl::PointXYZRGBNormal>(texturedCloud, rgb, "cloud");
	//~ viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "cloud");

	// Add mesh
	viewer->setBackgroundColor(1, 1, 1);
	viewer->addPolygonMesh(triangles, "meshes", 0);
	viewer->addCoordinateSystem(1.0);
	viewer->initCameraParameters();

	// Keep viewer open
	while (!viewer->wasStopped()) {
	  viewer->spinOnce(100);
	  boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}

    return 0;
  }  
