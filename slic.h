#ifndef SLIC_H
#define SLIC_H

/* slic.h.
 *
 * Written by: Pascal Mettes.
 * and adapted by using the DGtal library (see. http://libdgtal.org and https://github.com/DGtal-team/DGtal)
 *
 * This file contains the class elements of the class Slic. This class is an
 * implementation of the SLIC Superpixel algorithm by Achanta et al. [PAMI'12,
 * vol. 34, num. 11, pp. 2274-2282].
 *
 * This implementation is created for the specific purpose of creating
 * over-segmentations in an OpenCV-based environment.
 */

#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

#define NR_ITERATIONS 10

typedef struct Color_st {
    double r, g, b;
} Color_t;

typedef struct Center_st {
    double x, y, z;
    Color_t color;
} Center_t;


typedef std::vector<std::vector<std::vector<double>>> vec3dd;
typedef std::vector<std::vector<std::vector<int>>>    vec3di;
typedef std::vector<std::vector<std::vector<bool>>>   vec3db;


typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned int> Image3D;


class Slic {
  
private:
  /* The cluster assignments and distance values for each pixel. */
  vec3di clusters;
  vec3dd distances;
        
  /* The LAB and xy values of the centers. */
  std::vector<Center_t> centers;
  
  /* The number of occurrences of each center. */
  std::vector<int> center_counts;
        
  /* The step size per cluster, and the colour (nc) and distance (ns)
   * parameters. */
  int step, nc, ns;
        
  /* Compute the distance between a center and an individual pixel. */
  double compute_dist(int ci, DGtal::Z3i::Point pixel, Color_t color);
  /* Find the pixel with the lowest gradient in a 3x3 surrounding. */
  DGtal::Z3i::Point find_local_minimum(Image3D &image, DGtal::Z3i::Point center);
        
  /* Remove and initialize the 3d vectors. */
  void clear_data();
  void init_data(Image3D &image);
  
public:
  /**
   * Helpers 
   */
  /* Returns a Center structure set to point p and color at that point. */
  Center_t to_center_type(const Image3D& image, const DGtal::Z3i::Point& p);
  /* Returns the color at position (x,y) in the image */
  Color_t get_color_at(const Image3D &image, int x, int y, int z);

  /* Convert dgtal:Color to GrayScale */
  double color_to_grayscale(Color_t color);



  /* Class constructors and destructor. */
  Slic();
  ~Slic();

  /*return width of picture*/
  unsigned int get_width(Image3D &);

  /*return height of picture*/
  unsigned int get_height(Image3D &);

  /*return depth of picture*/
  unsigned int get_depth(Image3D &);

  /* Generate an over-segmentation for an image. */
  void generate_superpixels(Image3D &, int step, int nc);
  /* Enforce connectivity for an image. */
  void create_connectivity(Image3D &image);
        
  /* Draw functions. Resp. displayal of the centers and the contours. */
  void display_center_grid( Image3D &image, DGtal::Color &colour);
  void display_contours( Image3D &image, DGtal::Color &colour);
  void colour_with_cluster_means(Image3D &image);
};

#endif
