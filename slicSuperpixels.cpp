
/*
* test_slic.cpp.
*
* Written by: Pascal Mettes.
* adapted to the DGtal lib by B. Kerautret 24/04/15
* This file creates an over-segmentation of a provided image based on the SLIC
* superpixel algorithm, as implemented in slic.h and slic.cpp.
*/

#include <iostream>
#include <fstream>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/readers/VolReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include "DGtal/io/writers/VolWriter.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "slic.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


struct IdColor {
  Color operator()(const unsigned int& aValue) const {
    return DGtal::Color(aValue);
  }
};

int main(int argc, char** argv)
{

  po::options_description general_opt("Allowed options are: ");

  general_opt.add_options()
  ("help,h", "display this message")
  ("input,i", po::value<std::string>(), "input image file")
  ("output,o", po::value<std::string>(), "output regions file")
  ("n,n", po::value<int>()->default_value(1000), "the number of super pixels of the resulting image.")
  ("w,w", po::value<int>()->default_value(100), "the weight factor")
  ("noDisplayContour", "to avoid to display region contours.")
  ("displayCenters,c", "display the centers of the regions.")
  ("displayMeanColor,m", "display regions with mean color.");

  bool parseOK = true;
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  } catch (const std::exception& ex) {
    trace.info() << "Error checking program options: " << ex.what() << std::endl;
    parseOK = false;
  }

  po::notify(vm);

  if (vm.count("help") || argc <= 1 || !parseOK) {
    trace.info() << "Apply the slic superpixels segmentatation on input image" << std::endl << "Options: " << std::endl
                 << general_opt << "\n";
    return 0;
  }



  string nameImage = vm["input"].as<string>();
  string nameImageOutput = vm["output"].as<string>();
  Image3D image = VolReader<Image3D>::importVol(nameImage);


  Slic slic;
  unsigned int imageWidth = slic.get_width(image);
  unsigned int imageHeight = slic.get_height(image);
  unsigned int imageDepth = slic.get_depth(image);
  
  cout << imageWidth << " " << imageHeight << " " << imageDepth << endl;
  
  
//  for (size_t i = 0; i < imageWidth; i++) {
//    for (size_t j = 0; j < imageHeight; j++) {
//      Color_t c = slic.get_color_at(image, i, j);
//      std::cout << c.r << " " << c.g << " "<< c.b << " | ";
//    }
//    std::cout << std::endl;
//  }

  int nr = vm["n"].as<int>();
  int nc = vm["w"].as<int>();
  double step = sqrt((imageWidth * imageHeight * imageDepth) / (double) nr);

  /* Perform the SLIC superpixel algorithm. */

  slic.generate_superpixels(image, step, nc);
  slic.create_connectivity(image);
  DGtal::Color c(0, 0, 204);
  DGtal::Color c2(255, 100, 50);

  if (!vm.count("noDisplayContour")) {
    slic.display_contours(image, c);
  }

  if (vm.count("displayMeanColor")) {
    slic.colour_with_cluster_means(image);
  }

  if (vm.count("displayCenters")) {
    slic.display_center_grid(image, c2);
  }

  IdColor id;
  //VolWriter<Image3D, IdColor>::exportVol(nameImageOutput, image, id);

  
  return 1;
}
