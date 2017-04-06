#include "slic.h"


/*
 * Constructor. Nothing is done here.
 */
Slic::Slic()
{
  
}



/*
 * Destructor. Clear any present data.
 */
Slic::~Slic()
{
  clear_data();
}



/*
 * Clear the data as saved by the algorithm.
 *
 * Input : -
 * Output: -
 */
void Slic::clear_data()
{
  clusters.clear();
  distances.clear();
  centers.clear();
  center_counts.clear();
}


Center_t Slic::to_center_type(const Image3D& image, const DGtal::Z3i::Point& p)
{
  Center_t c;
  c.x = p[0];
  c.y = p[1];
  c.z = p[2];

  Color_t cl = get_color_at(image, p[0], p[1], p[2]);
  c.color = cl;
  
  return c;
}


Color_t Slic::get_color_at(const Image3D &image, int x, int y, int z)
{
  Color_t c;
  DGtal::Color color(image(DGtal::Z3i::Point(x, y, z)), 255);
  c.r = color.red();
  c.g = color.green();
  c.b = color.blue();
  
  return c;
}

/*
 * Converts to GrayScale
 */
double Slic::color_to_grayscale(Color_t color)
{
  return round(0.299*color.r + 0.587*color.g + 0.114*color.b);
}

/*
 * input: Image3D
 * output: width of input image
 */
unsigned int Slic::get_width(Image3D& image)
{
  return 1 + image.domain().upperBound()[0] - image.domain().lowerBound()[0];
}

/*
 * input: Image3D
 * output: height of input image
 */
unsigned int Slic::get_height(Image3D& image)
{
  return 1 + image.domain().upperBound()[1] - image.domain().lowerBound()[1];
}

/*
 * input: Image3D
 * output: depth of input image
 */
unsigned int Slic::get_depth(Image3D& image)
{
  return 1 + image.domain().upperBound()[2] - image.domain().lowerBound()[2];
}



/*
 * Initialize the cluster centers and initial values of the pixel-wise cluster
 * assignment and distance values.
 *
 * Input : The image (IplImage*).
 * Output: -
 */
void Slic::init_data(Image3D& image)
{
  unsigned int imageWidth = get_width(image);
  unsigned int imageHeight = get_height(image);
  unsigned int imageDepth = get_depth(image);

  /* Initialize the cluster and distance matrices. */
  for (size_t i = 0; i < imageWidth; i++) {
    std::vector<std::vector<int> > cluster_row;
    std::vector<std::vector<double> > dist_row;

    for (size_t j = 0; j < imageHeight; j++) {
      std::vector<int> cluster_aisle;
      std::vector<double> dist_aisle;
      
      for (size_t k = 0; k < imageDepth; k++) {
        cluster_aisle.push_back(-1);
        dist_aisle.push_back(std::numeric_limits<double>::max());
      }
      
      cluster_row.push_back(cluster_aisle);
      dist_row.push_back(dist_aisle);
    }

    clusters.push_back(cluster_row);
    distances.push_back(dist_row);
  }


  /* Initialize the centers and counters. */
  for (size_t i = step; i < imageWidth - step / 2; i += step) {
    for (size_t j = step; j < imageHeight - step / 2; j += step) {
      for (size_t k = step; k < imageDepth - step / 2; k += step) {
        Center_t center;
        /* Find the local minimum (gradient-wise). */
        DGtal::Z3i::Point nc = find_local_minimum(image, DGtal::Z3i::Point(i, j, k));

        center = to_center_type(image, nc);

        /* Append to vector of centers. */
        centers.push_back(center);
        center_counts.push_back(0);
      }
    }
  }
}



/*
 * Compute the distance between a cluster center and an individual pixel.
 *
 * Input : The cluster index (int), the pixel (CvPoint), and the Lab values of
 *         the pixel (CvScalar).
 * Output: The distance (double).
 */
double Slic::compute_dist(int ci, DGtal::Z3i::Point pixel, Color_t color)
{
  double dc = sqrt(pow(centers[ci].color.r - color.r, 2) 
                 + pow(centers[ci].color.g - color.g, 2) 
                 + pow(centers[ci].color.b - color.b, 2));
  double ds = sqrt(pow(centers[ci].x - pixel[0], 2) 
                 + pow(centers[ci].y - pixel[1], 2)
                 + pow(centers[ci].z - pixel[2], 2));

  return sqrt(pow(dc / nc, 2) + pow(ds / ns, 2));
}




/*
 * Find a local gradient minimum of a pixel in a 3x3 neighborhood. This
 * method is called upon initialization of the cluster centers.
 *
 * Input : The image (IplImage*) and the pixel center (CvPoint).
 * Output: The local gradient minimum (CvPoint).
 */
DGtal::Z3i::Point Slic::find_local_minimum(Image3D& image, DGtal::Z3i::Point center)
{
  double min_grad = std::numeric_limits<double>::max();
  DGtal::Z3i::Point loc_min = center;

  for (int i = center[0] - 1; i < center[0] + 3; i++) {
    for (int j = center[1] - 1; j < center[1] + 3; j++) {
        for (int k = center[1] - 1; k < center[2] + 3; k++) {
          Color_t c1 = get_color_at(image, i,   j+1, k+1);
          Color_t c2 = get_color_at(image, i+1, j, k+1);
          Color_t c3 = get_color_at(image, i,   j, k);

          /* Convert colour values to grayscale values. */
          double i1 = c1.r;
          double i2 = c2.r;
          double i3 = c3.r;

  //      double i1 = color_to_grayscale(c1);
  //      double i2 = color_to_grayscale(c2);
  //      double i3 = color_to_grayscale(c3);


        /* Compute horizontal and vertical gradients and keep track of the
         minimum. */
        if (sqrt(pow(i1 - i3, 2)) + sqrt(pow(i2 - i3, 2)) < min_grad) {
          min_grad = fabs(i1 - i3) + fabs(i2 - i3);
          loc_min[0] = i;
          loc_min[1] = j;
          loc_min[2] = k;
      }
    }
  }
}

  return loc_min;
}



/*
 * Compute the over-segmentation based on the step-size and relative weighting
 * of the pixel and colour values.
 *
 * Input : The Lab image (IplImage*), the stepsize (int), and the weight (int).
 * Output: -
 */
void Slic::generate_superpixels(Image3D& image, int step, int nc)
{
  this->step = step;
  this->nc = nc;
  this->ns = step;

  unsigned int imageWidth = get_width(image);
  unsigned int imageHeight = get_height(image);
  unsigned int imageDepth = get_depth(image);

  /* Clear previous data (if any), and re-initialize it. */
  clear_data();
  init_data(image);

  /* Run EM for 10 iterations (as prescribed by the algorithm). */
  for (int l = 0; l < NR_ITERATIONS; l++) {
    /* Reset distance values. */
    for (size_t i = 0; i < imageWidth; i++) {
      for (size_t j = 0; j < imageHeight; j++) {
        for (size_t k = 0; k < imageHeight; k++) {
          distances[i][j][k] = std::numeric_limits<double>::max();
        }
      }
    }


    for (size_t j = 0; j < centers.size(); j++) {
      /* Only compare to pixels in a 2 x step by 2 x step region. */
      for (int k = centers[j].x - step; k < centers[j].x + step; k++) {
        for (int l = centers[j].y - step; l < centers[j].y + step; l++) {
          for (int m = centers[j].z - step; m < centers[j].z + step; m++) {
            if (k >= 0 && (size_t)k < imageWidth && l >= 0 && (size_t)l < imageHeight && m >= 0 && (size_t)m < imageDepth) {
              Color_t color = get_color_at(image, k, l, m);
              double d = compute_dist(j, DGtal::Z3i::Point(k, l, m), color);

              /* Update cluster allocation if the cluster minimizes the
                 distance. */
              if (d < distances[k][l][m]) {
                distances[k][l][m] = d;
                clusters[k][l][m] = j;
              }
            }
          }
        }
      }

      /* Clear the center values. */
      for (size_t j = 0; j < centers.size(); j++) {
        centers[j].x = centers[j].y = centers[j].z= 0;
        centers[j].color.r = centers[j].color.g = centers[j].color.b = 0;
        center_counts[j] = 0;
      }

      /* Compute the new cluster centers. */
      for (size_t j = 0; j < imageWidth; j++) {
        for (size_t k = 0; k < imageHeight; k++) {
          for (size_t l = 0; l < imageDepth; l++) {
            int c_id = clusters[j][k][l];


            if (c_id != -1) {
              Color_t color = get_color_at(image, j, k, l);

              centers[c_id].color.r += color.r;
              centers[c_id].color.g += color.g;
              centers[c_id].color.b += color.b;
              centers[c_id].x += j;
              centers[c_id].y += k;
              centers[c_id].z += l;


              center_counts[c_id] += 1;
            }
          }
        }
      }

      /* Normalize the clusters. */
      for (size_t j = 0; j < centers.size(); j++) {
        centers[j].color.r /= center_counts[j];
        centers[j].color.g /= center_counts[j];
        centers[j].color.b /= center_counts[j];
        centers[j].x /= center_counts[j];
        centers[j].y /= center_counts[j];
        centers[j].z /= center_counts[j];
      }
    }
  }
}



/*
 * Enforce connectivity of the superpixels. This part is not actively discussed
 * in the paper, but forms an active part of the implementation of the authors
 * of the paper.
 *
 * Input : The image (Image3D&).
 * Output: -
 */
void Slic::create_connectivity(Image3D& image)
{
  int label = 0, adjlabel = 0;
  unsigned int imageWidth = get_width(image);
  unsigned int imageHeight = get_height(image);
  unsigned int imageDepth = get_depth(image);

  const int lims = (imageWidth * imageHeight * imageDepth) / (centers.size());

  const int dx6[6] = { -1,  0,  1,  0,  0,  0};
  const int dy6[6] = { 0 , -1,  0,  1,  0,  0};
  const int dz6[6] = { 0 ,  0,  0,  0, -1,  0};

  /* Initialize the new cluster matrix. */
  vec3di new_clusters;

  for (size_t i = 0; i < imageWidth; i++) {
    std::vector<std::vector<int>> nc;

    for (size_t j = 0; j < imageHeight; j++) {
      std::vector<int> tmpnc;
      tmpnc.push_back(-1);
        for (size_t K = 0; K < imageDepth; K++) {
          nc.push_back(tmpnc);
        }

    }
    new_clusters.push_back(nc);
  }

  for (size_t i = 0; i < imageWidth; i++) {
    for (size_t j = 0; j < imageHeight; j++) {
      for (size_t l = 0; l < imageDepth; l++) {
        if (new_clusters[i][j][l] == -1) {
          std::vector<DGtal::Z3i::Point> elements;
          elements.push_back(DGtal::Z3i::Point(i, j, l));

          /* Find an adjacent label, for possible use later. */
          for (int k = 0; k < 4; k++) {
            int x = elements[0][0] + dx6[k], y = elements[0][1] + dy6[k], z = elements[0][2] + dz6[k];

            if (x >= 0 && (size_t)x < imageWidth && y >= 0 && (size_t)y < imageHeight && z >= 0 && (size_t)z < imageDepth ) {
              if (new_clusters[x][y][z] >= 0) {
                adjlabel = new_clusters[x][y][z];
              }
            }
          }

          int count = 1;

          for (int c = 0; c < count; c++) {
            for (int k = 0; k < 4; k++) {
              int x = elements[c][0] + dx6[k], y = elements[c][1] + dy6[k], z = elements[c][2] + dz6[k];

              if (x >= 0 && (size_t)x < imageWidth && y >= 0 && (size_t)y < imageHeight && z >= 0 && (size_t)z < imageDepth) {
                if (new_clusters[x][y][z] == -1 && clusters[i][j][l] == clusters[x][y][z]) {
                  elements.push_back(DGtal::Z3i::Point(x, y, z));
                  new_clusters[x][y][z] = label;
                  count += 1;
                }
              }
            }
          }

          /* Use the earlier found adjacent label if a segment size is
             smaller than a limit. */
          if (count <= lims >> 2) {
            for (int c = 0; c < count; c++) {
              new_clusters[elements[c][0]][elements[c][1]][elements[c][2]] = adjlabel;
            }

            label -= 1;
          }

          label += 1;
        }
      }
    }
  }
}



/*
 * Display the cluster centers.
 *
 * Input : The image to display upon (IplImage*) and the colour (CvScalar).
 * Output: -
 */
void Slic::display_center_grid(Image3D& image, DGtal::Color& colour)
{
  for (size_t i = 0; i < centers.size(); i++) {
    image.setValue(DGtal::Z3i::Point(centers[i].x, centers[i].y, centers[i].z), colour.getRGB());
  }
}


/*
 * Display a single pixel wide contour around the clusters.
 *
 * Input : The target image (IplImage*) and contour colour (CvScalar).
 * Output: -
 */
void Slic::display_contours(Image3D& image, DGtal::Color& colour)
{
  const int dx14[14] = {-1, -1,  0,  1,  1,  1,  0, -1,  0,  0,  0,  0,  0,  0};
  const int dy14[14] = { 0, -1, -1, -1,  0,  1,  1,  1, -1,  0,  1, -1,  0,  1};
  const int dz14[14] = { 0,  0,  0,  0,  0,  0,  0,  0, -1, -1, -1,  1,  1,  1};
  unsigned int imageWidth = get_width(image);
  unsigned int imageHeight = get_height(image);
  unsigned int imageDepth = get_depth(image);

  /* Initialize the contour vector and the matrix detailing whether a pixel
   * is already taken to be a contour. */
  std::vector<DGtal::Z3i::Point> contours;
  vec3db istaken;

  for (size_t i = 0; i < imageWidth; i++) {
    std::vector<std::vector<bool>> nb_aisle;

    for (size_t j = 0; j < imageHeight; j++) {
      std::vector<bool> nb;
      
      for (size_t k = 0; k < imageDepth; k++) {
        nb.push_back(false);
      }
      nb_aisle.push_back(nb);
    }
    istaken.push_back(nb_aisle);
  }

  /* Go through all the pixels. */
  for (size_t i = 0; i < imageWidth; i++) {
    for (size_t j = 0; j < imageHeight; j++) {
      for (size_t k = 0; k < imageDepth; k++) {
        int nr_p = 0;

        /* Compare the pixel to its 14 neighbors. */
        for (int k = 0; k < 14; k++) {
          int x = i + dx14[k]; 
          int y = j + dy14[k];
          int z = k + dz14[k];

          if (x >= 0 && (size_t)x < imageWidth && 
              y >= 0 && (size_t)y < imageHeight &&
              z >= 0 && (size_t)z < imageDepth) 
          {
            if (istaken[x][y][z] == false && clusters[i][j][k] != clusters[x][y][z]) {
              nr_p += 1;
            }
          }
        }

        /* Add the pixel to the contour list if desired. */
        if (nr_p >= 2) {
          contours.push_back(DGtal::Z3i::Point(i, j, k));
          istaken[i][j][k] = true;
        }
      }
    }
  }

  /* Draw the contour pixels. */
  for (size_t i = 0; i < contours.size(); i++) {
    image.setValue(DGtal::Z3i::Point(contours[i][0], contours[i][1], contours[i][2]), colour.getRGB());
  }
}





/*
 * Give the pixels of each cluster the same colour values. The specified colour
 * is the mean RGB colour per cluster.
 *
 * Input : The target image (IplImage*).
 * Output: -
 */
void Slic::colour_with_cluster_means(Image3D& image)
{
  std::vector<Color_t> colors(centers.size());
  unsigned int imageWidth = get_width(image);
  unsigned int imageHeight = get_height(image);
  unsigned int imageDepth = get_depth(image);

  /* Gather the colour values per cluster. */
  for (size_t i = 0; i < imageWidth; i++) {
    for (size_t j = 0; j < imageHeight; j++) {
      for (size_t k = 0; k < imageDepth; k++) {
        int index = clusters[i][j][k];
        Color_t color = get_color_at(image, i, j, k);

        colors[index].r += color.r;
        colors[index].g += color.g;
        colors[index].b += color.b;
      }
    }
  }

  /* Divide by the number of pixels per cluster to get the mean colour. */
  for (size_t i = 0; i < colors.size(); i++) {
    colors[i].r /= center_counts[i];
    colors[i].g /= center_counts[i];
    colors[i].b /= center_counts[i];
  }

  /* Fill in. */
  for (size_t i = 0; i < imageWidth; i++) {
    for (size_t j = 0; j < imageHeight; j++) {
      for (size_t k = 0; k < imageDepth; k++) {
        DGtal::Color col;
        Color_t ncolor = colors[clusters[i][j][k]];
        col.setRGBf(ncolor.r / 255.0, ncolor.g / 255.0, ncolor.b / 255.0 );
        image.setValue(DGtal::Z3i::Point(i, j, k), col.getRGB());
      }
    }
  }
}






