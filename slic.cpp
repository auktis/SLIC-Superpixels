#include "slic.h"


/*
 * Constructor. Nothing is done here.
 */
Slic::Slic() {

}



/*
 * Destructor. Clear any present data.
 */
Slic::~Slic() {
  clear_data();
}



/*
 * Clear the data as saved by the algorithm.
 *
 * Input : -
 * Output: -
 */
void Slic::clear_data() {
  clusters.clear();
  distances.clear();
  centers.clear();
  center_counts.clear();
}



/*
 * Initialize the cluster centers and initial values of the pixel-wise cluster
 * assignment and distance values.
 *
 * Input : The image (IplImage*).
 * Output: -
 */
void Slic::init_data(Image2D &image) {
  unsigned int imageWidth = 1+image.domain().upperBound()[0]-image.domain().lowerBound()[0];
  unsigned int imageHeight = 1+image.domain().upperBound()[1]-image.domain().lowerBound()[1];

  /* Initialize the cluster and distance matrices. */
  for (int i = 0; i < imageWidth; i++) { 
    std::vector<int> cr;
    std::vector<double> dr;
    for (int j = 0; j < imageHeight; j++) {
      cr.push_back(-1);
      dr.push_back(std::numeric_limits<double>::max());
    }
    clusters.push_back(cr);
    distances.push_back(dr);
  }
  ToCvScalarFct toCv;
  /* Initialize the centers and counters. */
  for (int i = step; i < imageWidth - step/2; i += step) {
    for (int j = step; j < imageHeight - step/2; j += step) {
      std::vector<double> center;
      /* Find the local minimum (gradient-wise). */
      DGtal::Z2i::Point nc = find_local_minimum(image, DGtal::Z2i::Point(i,j));
      CvScalar colour = toCv(image(nc));
    
      /* Generate the center vector. */
      center.push_back(colour.val[0]);
      center.push_back(colour.val[1]);
      center.push_back(colour.val[2]);
      center.push_back(nc[0]);
      center.push_back(nc[1]);
            
      /* Append to vector of centers. */
      centers.push_back(center);
      center_counts.push_back(0);
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
double Slic::compute_dist(int ci, DGtal::Z2i::Point pixel, CvScalar colour) {
  double dc = sqrt(pow(centers[ci][0] - colour.val[0], 2) + pow(centers[ci][1] - colour.val[1], 2) + pow(centers[ci][2] - colour.val[2], 2));
  double ds = sqrt(pow(centers[ci][3] - pixel[0], 2) + pow(centers[ci][4] - pixel[1], 2));
  
  return sqrt(pow(dc / nc, 2) + pow(ds / ns, 2));
}




/*
 * Find a local gradient minimum of a pixel in a 3x3 neighbourhood. This
 * method is called upon initialization of the cluster centers.
 *
 * Input : The image (IplImage*) and the pixel center (CvPoint).
 * Output: The local gradient minimum (CvPoint).
 */
DGtal::Z2i::Point Slic::find_local_minimum(Image2D &image, DGtal::Z2i::Point center) {

  double min_grad = std::numeric_limits<double>::max();
  DGtal::Z2i::Point loc_min = center;
  ToCvScalarFct toScal;

  for (int i = center[0]-1; i < center[0]+2; i++) {
    for (int j = center[1]-1; j < center[1]+2; j++) {
      CvScalar c1 = toScal(image(DGtal::Z2i::Point(i, j+1)));
      CvScalar c2 = toScal(image(DGtal::Z2i::Point(i+1, j)));
      CvScalar c3 = toScal(image(DGtal::Z2i::Point(i, j)));
      /* Convert colour values to grayscale values. */
      double i1 = c1.val[0];
      double i2 = c2.val[0];
      double i3 = c3.val[0];
      /*double i1 = c1.val[0] * 0.11 + c1.val[1] * 0.59 + c1.val[2] * 0.3;
       double i2 = c2.val[0] * 0.11 + c2.val[1] * 0.59 + c2.val[2] * 0.3;
       double i3 = c3.val[0] * 0.11 + c3.val[1] * 0.59 + c3.val[2] * 0.3;*/
      
      /* Compute horizontal and vertical gradients and keep track of the
       minimum. */
      if (sqrt(pow(i1 - i3, 2)) + sqrt(pow(i2 - i3,2)) < min_grad) {
        min_grad = fabs(i1 - i3) + fabs(i2 - i3);
        loc_min[0] = i;
        loc_min[1] = j;
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
void Slic::generate_superpixels(Image2D &image, int step, int nc) {
  this->step = step;
  this->nc = nc;
  this->ns = step;
  unsigned int imageWidth = 1+image.domain().upperBound()[0]-image.domain().lowerBound()[0];
  unsigned int imageHeight = 1+image.domain().upperBound()[1]-image.domain().lowerBound()[1];
    
  /* Clear previous data (if any), and re-initialize it. */
  clear_data();
  init_data(image);
    
  /* Run EM for 10 iterations (as prescribed by the algorithm). */
  for (int i = 0; i < NR_ITERATIONS; i++) {
    /* Reset distance values. */
    for (int j = 0; j < imageWidth; j++) {
      for (int k = 0;k < imageHeight; k++) {
        distances[j][k] = std::numeric_limits<double>::max();
      }
    }
    ToCvScalarFct toScal;
    for (int j = 0; j < (int) centers.size(); j++) {
      /* Only compare to pixels in a 2 x step by 2 x step region. */
      for (int k = centers[j][3] - step; k < centers[j][3] + step; k++) {
        for (int l = centers[j][4] - step; l < centers[j][4] + step; l++) {
          if (k >= 0 && k < imageWidth && l >= 0 && l < imageHeight) {
            unsigned int colour = image(DGtal::Z2i::Point(k, l));
            double d = compute_dist(j, DGtal::Z2i::Point(k,l), toScal(colour)); //compute_dist(j, cvPoint(k,l), colour);
                        
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
    for (int j = 0; j < (int) centers.size(); j++) {
      centers[j][0] = centers[j][1] = centers[j][2] = centers[j][3] = centers[j][4] = 0;
      center_counts[j] = 0;
    }
        
    /* Compute the new cluster centers. */
    for (int j = 0; j < imageWidth; j++) {
      for (int k = 0; k < imageHeight; k++) {
        int c_id = clusters[j][k];
                
        if (c_id != -1) {
          CvScalar colour = toScal(image(DGtal::Z2i::Point(j, k)));
          
          centers[c_id][0] += colour.val[0];
          centers[c_id][1] += colour.val[1];
          centers[c_id][2] += colour.val[2];
          centers[c_id][3] += j;
          centers[c_id][4] += k;
                  
          center_counts[c_id] += 1;
        }
      }
    }

    /* Normalize the clusters. */
    for (int j = 0; j < (int) centers.size(); j++) {
      centers[j][0] /= center_counts[j];
      centers[j][1] /= center_counts[j];
      centers[j][2] /= center_counts[j];
      centers[j][3] /= center_counts[j];
      centers[j][4] /= center_counts[j];
    }
  }
}




/*
 * Enforce connectivity of the superpixels. This part is not actively discussed
 * in the paper, but forms an active part of the implementation of the authors
 * of the paper.
 *
 * Input : The image (IplImage*).
 * Output: -
 */
void Slic::create_connectivity(Image2D &image) {
  int label = 0, adjlabel = 0;
  unsigned int imageWidth = 1+image.domain().upperBound()[0]-image.domain().lowerBound()[0];
  unsigned int imageHeight = 1+image.domain().upperBound()[1]-image.domain().lowerBound()[1];
    
  const int lims = (imageWidth * imageHeight) / ((int)centers.size());
    
  const int dx4[4] = {-1,  0,  1,  0};
  const int dy4[4] = { 0, -1,  0,  1};
    
  /* Initialize the new cluster matrix. */
  vec2di new_clusters;
  for (int i = 0; i < imageWidth; i++) { 
    std::vector<int> nc;
    for (int j = 0; j < imageHeight; j++) {
      nc.push_back(-1);
    }
    new_clusters.push_back(nc);
  }

  for (int i = 0; i < imageWidth; i++) {
    for (int j = 0; j < imageHeight; j++) {
      if (new_clusters[i][j] == -1) {
        std::vector<DGtal::Z2i::Point> elements;
        elements.push_back(DGtal::Z2i::Point(i, j));
              
        /* Find an adjacent label, for possible use later. */
        for (int k = 0; k < 4; k++) {
          int x = elements[0][0] + dx4[k], y = elements[0][1] + dy4[k];
                
          if (x >= 0 && x < imageWidth && y >= 0 && y < imageHeight) {
            if (new_clusters[x][y] >= 0) {
              adjlabel = new_clusters[x][y];
            }
          }
        }
              
        int count = 1;
        for (int c = 0; c < count; c++) {
          for (int k = 0; k < 4; k++) {
            int x = elements[c][0] + dx4[k], y = elements[c][1] + dy4[k];
                        
            if (x >= 0 && x < imageWidth && y >= 0 && y < imageHeight) {
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
void Slic::display_center_grid(Image2D &image, DGtal::Color &colour) {
  for (int i = 0; i < (int) centers.size(); i++) {
    image.setValue ( DGtal::Z2i::Point(centers[i][3], centers[i][4]), colour.getRGB());
  }
}


/*
 * Display a single pixel wide contour around the clusters.
 *
 * Input : The target image (IplImage*) and contour colour (CvScalar).
 * Output: -
 */
void Slic::display_contours(Image2D &image, DGtal::Color &colour) {
  const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
  const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
  unsigned int imageWidth = 1+image.domain().upperBound()[0]-image.domain().lowerBound()[0];
  unsigned int imageHeight = 1+image.domain().upperBound()[1]-image.domain().lowerBound()[1];
        
  /* Initialize the contour vector and the matrix detailing whether a pixel
   * is already taken to be a contour. */
  std::vector<DGtal::Z2i::Point> contours;
  vec2db istaken;
  for (int i = 0; i < imageWidth; i++) { 
    std::vector<bool> nb;
    for (int j = 0; j < imageHeight; j++) {
      nb.push_back(false);
    }
    istaken.push_back(nb);
  }

  /* Go through all the pixels. */
  for (int i = 0; i < imageWidth; i++) {
    for (int j = 0; j < imageHeight; j++) {
      int nr_p = 0;
            
      /* Compare the pixel to its 8 neighbours. */
      for (int k = 0; k < 8; k++) {
        int x = i + dx8[k], y = j + dy8[k];
                
        if (x >= 0 && x < imageWidth && y >= 0 && y < imageHeight) {
          if (istaken[x][y] == false && clusters[i][j] != clusters[x][y]) {
            nr_p += 1;
          }
        }
      }
            
      /* Add the pixel to the contour list if desired. */
      if (nr_p >= 2) {
        contours.push_back(DGtal::Z2i::Point(i,j));
        istaken[i][j] = true;
      }
    }
  }
    
  /* Draw the contour pixels. */
  for (int i = 0; i < (int)contours.size(); i++) {
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
void Slic::colour_with_cluster_means(Image2D &image) {
  std::vector<CvScalar> colours(centers.size());
  unsigned int imageWidth = 1+image.domain().upperBound()[0]-image.domain().lowerBound()[0];
  unsigned int imageHeight = 1+image.domain().upperBound()[1]-image.domain().lowerBound()[1];
    
  /* Gather the colour values per cluster. */
  for (int i = 0; i < imageWidth; i++) {
    for (int j = 0; j < imageHeight; j++) {
      int index = clusters[i][j];
      ToCvScalarFct toScal;
      CvScalar colour = toScal(image(DGtal::Z2i::Point(i,j)));
      
      colours[index].val[0] += colour.val[0];
      colours[index].val[1] += colour.val[1];
      colours[index].val[2] += colour.val[2];
    }
  }
    
  /* Divide by the number of pixels per cluster to get the mean colour. */
  for (int i = 0; i < (int)colours.size(); i++) {
    colours[i].val[0] /= center_counts[i];
    colours[i].val[1] /= center_counts[i];
    colours[i].val[2] /= center_counts[i];
  }
    
  /* Fill in. */
  for (int i = 0; i < imageWidth; i++) {
    for (int j = 0; j < imageHeight; j++) {
      DGtal::Color col;
      CvScalar ncolour = colours[clusters[i][j]];
      col.setRGBf(ncolour.val[0]/255.0, ncolour.val[1]/255.0, ncolour.val[2]/255.0 );
      image.setValue(DGtal::Z2i::Point(i, j), col.getRGB());
    }
  }
}






