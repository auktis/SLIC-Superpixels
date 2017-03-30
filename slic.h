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


//typedef struct CvScalar {
//  double val[4];
//} CvScalar;
//
//
//
//typedef struct ToCvScalarFct{
//  inline
//  CvScalar operator() (unsigned int aVal) const
//  {
//    CvScalar result;
//    DGtal::Color c (aVal, 255);
//    result.val[0] = (double)(c.red());
//    result.val[1] = (double)(c.green());
//    result.val[2] = (double)(c.blue());
//    return result;
//  }
//  
//} ToCvScalarFct;

#define NR_ITERATIONS 10

typedef struct Color_st {
    double r, g, b;
} Color_t;

typedef struct Center_st {
    double x, y;
    Color_t color;
} Center_t;


typedef std::vector<std::vector<double> > vec2dd;
typedef std::vector<std::vector<int> >    vec2di;
typedef std::vector<std::vector<bool> >   vec2db;


typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned int> Image2D;


class Slic {
  
private:
  /* The cluster assignments and distance values for each pixel. */
  vec2di clusters;
  vec2dd distances;
        
  /* The LAB and xy values of the centers. */
  std::vector<Center_t> centers;
  
  /* The number of occurrences of each center. */
  std::vector<int> center_counts;
        
  /* The step size per cluster, and the colour (nc) and distance (ns)
   * parameters. */
  int step, nc, ns;
        
  /* Compute the distance between a center and an individual pixel. */
  double compute_dist(int ci, DGtal::Z2i::Point pixel, Color_t color);
  /* Find the pixel with the lowest gradient in a 3x3 surrounding. */
  DGtal::Z2i::Point find_local_minimum(Image2D &image, DGtal::Z2i::Point center);
        
  /* Remove and initialize the 2d vectors. */
  void clear_data();
  void init_data(Image2D &image);
  

  /**
   * Helpers 
   */
  /* Returns a Center structure set to point p and color at that point. */
  Center_t to_center_type(const Image2D& image, const DGtal::Z2i::Point& p);
  /* Returns the color at position (x,y) in the image */
  Color_t get_color_at(const Image2D &image, int x, int y);

  /* Convert dgtal:Color to GrayScale */
  double color_to_grayscale(Color_t color);


public:
  /* Class constructors and destructor. */
  Slic();
  ~Slic();

  /*return width of picture*/
  unsigned int get_width(Image2D &);

  /*return height of picture*/
  unsigned int get_height(Image2D &);

  /* Generate an over-segmentation for an image. */
  void generate_superpixels(Image2D &, int step, int nc);
  /* Enforce connectivity for an image. */
  void create_connectivity(Image2D &image);
        
  /* Draw functions. Resp. displayal of the centers and the contours. */
  void display_center_grid( Image2D &image, DGtal::Color &colour);
  void display_contours( Image2D &image, DGtal::Color &colour);
  void colour_with_cluster_means(Image2D &image);
};

#endif
