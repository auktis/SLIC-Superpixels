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

Center_t Slic::toCenterType(const Image2D& image, const DGtal::Z2i::Point& p)
{
  Center_t c;
  c.x = p[0];
  c.y = p[1];
  Color_t cl = getColorAt(image, p[0], p[1]);
  c.color = cl;
  
  return c;
}

Color_t Slic::getColorAt(const Image2D &image, int x, int y)
{
  Color_t c;
  DGtal::Color color(image(DGtal::Z2i::Point(x, y)), 255);
  c.r = color.red();
  c.g = color.green();
  c.b = color.blue();
  
  return c;
}




/*
 * Initialize the cluster centers and initial values of the pixel-wise cluster
 * assignment and distance values.
 *
 * Input : The image (IplImage*).
 * Output: -
 */
void Slic::init_data(Image2D& image)
{
  unsigned int imageWidth = get_width(image);
  unsigned int imageHeight = get_height(image);

  /* Initialize the cluster and distance matrices. */
  for (size_t i = 0; i < imageWidth; i++) {
    std::vector<int> cluster_row;
    std::vector<double> dist_row;

    for (size_t j = 0; j < imageHeight; j++) {
      cluster_row.push_back(-1);
      dist_row.push_back(std::numeric_limits<double>::max());
    }

    clusters.push_back(cluster_row);
    distances.push_back(dist_row);
  }


  /* Initialize the centers and counters. */
  for (size_t i = step; i < imageWidth - step / 2; i += step) {
    for (size_t j = step; j < imageHeight - step / 2; j += step) {
      Center_t center;
      /* Find the local minimum (gradient-wise). */
      DGtal::Z2i::Point nc = find_local_minimum(image, DGtal::Z2i::Point(i, j));

      center = toCenterType(image, nc);

      /* Append to vector of centers. */
      centers.push_back(center);
      center_counts.push_back(0);
    }
  }
}




/*
 * Convert dgtal:Color to GrayScale
 */
double Slic::colorToGrayscale(DGtal::Color &color)
{
  return round(0.299*color.red() + 0.587*color.green() + 0.114*color.blue());
}




/*
 * Compute the distance between a cluster center and an individual pixel.
 *
 * Input : The cluster index (int), the pixel (CvPoint), and the Lab values of
 *         the pixel (CvScalar).
 * Output: The distance (double).
 */
double Slic::compute_dist(int ci, DGtal::Z2i::Point pixel, Color_t color)
{
  double dc = sqrt(pow(centers[ci].color.r - color.r, 2) 
                 + pow(centers[ci].color.g - color.g, 2) 
                 + pow(centers[ci].color.b - color.b, 2));
  double ds = sqrt(pow(centers[ci].x - pixel[0], 2) + pow(centers[ci].y - pixel[1], 2));

  return sqrt(pow(dc / nc, 2) + pow(ds / ns, 2));
}




/*
 * Find a local gradient minimum of a pixel in a 3x3 neighbourhood. This
 * method is called upon initialization of the cluster centers.
 *
 * Input : The image (IplImage*) and the pixel center (CvPoint).
 * Output: The local gradient minimum (CvPoint).
 */
DGtal::Z2i::Point Slic::find_local_minimum(Image2D& image, DGtal::Z2i::Point center)
{
  double min_grad = std::numeric_limits<double>::max();
  DGtal::Z2i::Point loc_min = center;

  for (int i = center[0] - 1; i < center[0] + 2; i++) {
    for (int j = center[1] - 1; j < center[1] + 2; j++) {
      Color_t c1 = getColorAt(image, i,   j+1);
      Color_t c2 = getColorAt(image, i+1, j  );
      Color_t c3 = getColorAt(image, i,   j  );
      /* Convert colour values to grayscale values. */
      double i1 = c1.r;
      double i2 = c2.r;
      double i3 = c3.r;
      /*double i1 = c1.val[0] * 0.11 + c1.val[1] * 0.59 + c1.val[2] * 0.3;
       double i2 = c2.val[0] * 0.11 + c2.val[1] * 0.59 + c2.val[2] * 0.3;
       double i3 = c3.val[0] * 0.11 + c3.val[1] * 0.59 + c3.val[2] * 0.3;*/

      /* Compute horizontal and vertical gradients and keep track of the
       minimum. */
      if (sqrt(pow(i1 - i3, 2)) + sqrt(pow(i2 - i3, 2)) < min_grad) {
        min_grad = fabs(i1 - i3) + fabs(i2 - i3);
        loc_min[0] = i;
        loc_min[1] = j;
      }
    }
  }


  return loc_min;
}



/*
 * input: Image2D
 * output: width of input image
 */
unsigned int Slic::get_width(Image2D& image)
{
  return 1 + image.domain().upperBound()[0] - image.domain().lowerBound()[0];
}



/*
 * input: Image2D
 * output: height of input image
 */
unsigned int Slic::get_height(Image2D& image)
{
  return 1 + image.domain().upperBound()[1] - image.domain().lowerBound()[1];
}


/*
 * Compute the over-segmentation based on the step-size and relative weighting
 * of the pixel and colour values.
 *
 * Input : The Lab image (IplImage*), the stepsize (int), and the weight (int).
 * Output: -
 */
void Slic::generate_superpixels(Image2D& image, int step, int nc)
{
  this->step = step;
  this->nc = nc;
  this->ns = step;

  unsigned int imageWidth = get_width(image);
  unsigned int imageHeight = get_height(image);

  /* Clear previous data (if any), and re-initialize it. */
  clear_data();
  init_data(image);

  /* Run EM for 10 iterations (as prescribed by the algorithm). */
  for (int i = 0; i < NR_ITERATIONS; i++) {
    /* Reset distance values. */
    for (size_t j = 0; j < imageWidth; j++) {
      for (size_t k = 0; k < imageHeight; k++) {
        distances[j][k] = std::numeric_limits<double>::max();
      }
    }


    for (size_t j = 0; j < centers.size(); j++) {
      /* Only compare to pixels in a 2 x step by 2 x step region. */
      for (int k = centers[j].x - step; k < centers[j].x + step; k++) {
        for (int l = centers[j].y - step; l < centers[j].y + step; l++) {
          if (k >= 0 && (size_t)k < imageWidth && l >= 0 && (size_t)l < imageHeight) {
            Color_t color = getColorAt(image, k, l);
            double d = compute_dist(j, DGtal::Z2i::Point(k, l), color);

            /* Update cluster allocation if the cluster minimizes the
               distance. */
            if (d < distances[k][l]) {
              distances[k][l] = d;
              clusters[k][l] = j;
            }
          }
        }
      }
    }

    /* Clear the center values. */
    for (size_t j = 0; j < centers.size(); j++) {
      centers[j].x = centers[j].y = 0;
      centers[j].color.r = centers[j].color.g = centers[j].color.b = 0;
      center_counts[j] = 0;
    }

    /* Compute the new cluster centers. */
    for (size_t j = 0; j < imageWidth; j++) {
      for (size_t k = 0; k < imageHeight; k++) {
        int c_id = clusters[j][k];

        if (c_id != -1) {
          Color_t color = getColorAt(image, j, k);

          centers[c_id].color.r += color.r;
          centers[c_id].color.g += color.g;
          centers[c_id].color.b += color.b;
          centers[c_id].x += j;
          centers[c_id].y += k;
          
          center_counts[c_id] += 1;
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
    }
  }
}




/*
 * Enforce connectivity of the superpixels. This part is not actively discussed
 * in the paper, but forms an active part of the implementation of the authors
 * of the paper.
 *
 * Input : The image (Image2D&).
 * Output: -
 */
void Slic::create_connectivity(Image2D& image)
{
  int label = 0, adjlabel = 0;
  unsigned int imageWidth = get_width(image);
  unsigned int imageHeight = get_height(image);

  const int lims = (imageWidth * imageHeight) / (centers.size());

  const int dx4[4] = { -1,  0,  1,  0};
  const int dy4[4] = { 0 , -1,  0,  1};

  /* Initialize the new cluster matrix. */
  vec2di new_clusters;

  for (size_t i = 0; i < imageWidth; i++) {
    std::vector<int> nc;

    for (size_t j = 0; j < imageHeight; j++) {
      nc.push_back(-1);
    }

    new_clusters.push_back(nc);
  }

  for (size_t i = 0; i < imageWidth; i++) {
    for (size_t j = 0; j < imageHeight; j++) {
      if (new_clusters[i][j] == -1) {
        std::vector<DGtal::Z2i::Point> elements;
        elements.push_back(DGtal::Z2i::Point(i, j));

        /* Find an adjacent label, for possible use later. */
        for (int k = 0; k < 4; k++) {
          int x = elements[0][0] + dx4[k], y = elements[0][1] + dy4[k];

          if (x >= 0 && (size_t)x < imageWidth && y >= 0 && (size_t)y < imageHeight) {
            if (new_clusters[x][y] >= 0) {
              adjlabel = new_clusters[x][y];
            }
          }
        }

        int count = 1;

        for (int c = 0; c < count; c++) {
          for (int k = 0; k < 4; k++) {
            int x = elements[c][0] + dx4[k], y = elements[c][1] + dy4[k];

            if (x >= 0 && (size_t)x < imageWidth && y >= 0 && (size_t)y < imageHeight) {
              if (new_clusters[x][y] == -1 && clusters[i][j] == clusters[x][y]) {
                elements.push_back(DGtal::Z2i::Point(x, y));
                new_clusters[x][y] = label;
                count += 1;
              }
            }
          }
        }

        /* Use the earlier found adjacent label if a segment size is
           smaller than a limit. */
        if (count <= lims >> 2) {
          for (int c = 0; c < count; c++) {
            new_clusters[elements[c][0]][elements[c][1]] = adjlabel;
          }

          label -= 1;
        }

        label += 1;
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
void Slic::display_center_grid(Image2D& image, DGtal::Color& colour)
{
  for (size_t i = 0; i < centers.size(); i++) {
    image.setValue(DGtal::Z2i::Point(centers[i].x, centers[i].y), colour.getRGB());
  }
}


/*
 * Display a single pixel wide contour around the clusters.
 *
 * Input : The target image (IplImage*) and contour colour (CvScalar).
 * Output: -
 */
void Slic::display_contours(Image2D& image, DGtal::Color& colour)
{
  const int dx8[8] = { -1, -1,  0,  1, 1, 1, 0, -1};
  const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
  unsigned int imageWidth = get_width(image);
  unsigned int imageHeight = get_height(image);

  /* Initialize the contour vector and the matrix detailing whether a pixel
   * is already taken to be a contour. */
  std::vector<DGtal::Z2i::Point> contours;
  vec2db istaken;

  for (size_t i = 0; i < imageWidth; i++) {
    std::vector<bool> nb;

    for (size_t j = 0; j < imageHeight; j++) {
      nb.push_back(false);
    }

    istaken.push_back(nb);
  }

  /* Go through all the pixels. */
  for (size_t i = 0; i < imageWidth; i++) {
    for (size_t j = 0; j < imageHeight; j++) {
      int nr_p = 0;

      /* Compare the pixel to its 8 neighbours. */
      for (int k = 0; k < 8; k++) {
        int x = i + dx8[k], y = j + dy8[k];

        if (x >= 0 && (size_t)x < imageWidth && y >= 0 && (size_t)y < imageHeight) {
          if (istaken[x][y] == false && clusters[i][j] != clusters[x][y]) {
            nr_p += 1;
          }
        }
      }

      /* Add the pixel to the contour list if desired. */
      if (nr_p >= 2) {
        contours.push_back(DGtal::Z2i::Point(i, j));
        istaken[i][j] = true;
      }
    }
  }

  /* Draw the contour pixels. */
  for (size_t i = 0; i < contours.size(); i++) {
    image.setValue(DGtal::Z2i::Point(contours[i][0], contours[i][1]), colour.getRGB());

  }
}





/*
 * Give the pixels of each cluster the same colour values. The specified colour
 * is the mean RGB colour per cluster.
 *
 * Input : The target image (IplImage*).
 * Output: -
 */
void Slic::colour_with_cluster_means(Image2D& image)
{
  std::vector<Color_t> colors(centers.size());
  unsigned int imageWidth = get_width(image);
  unsigned int imageHeight = get_height(image);

  /* Gather the colour values per cluster. */
  for (size_t i = 0; i < imageWidth; i++) {
    for (size_t j = 0; j < imageHeight; j++) {
      int index = clusters[i][j];
      Color_t color = getColorAt(image, i, j);

      colors[index].r += color.r;
      colors[index].g += color.g;
      colors[index].b += color.b;
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
      DGtal::Color col;
      Color_t ncolor = colors[clusters[i][j]];
      col.setRGBf(ncolor.r / 255.0, ncolor.g / 255.0, ncolor.b / 255.0 );
      image.setValue(DGtal::Z2i::Point(i, j), col.getRGB());
    }
  }
}






