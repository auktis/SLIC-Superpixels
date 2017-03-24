
/*
* test_slic.cpp.
*
* Written by: Pascal Mettes.
* adated to the DGtal lib by B. Kerautret 24/04/15
* This file creates an over-segmentation of a provided image based on the SLIC
* superpixel algorithm, as implemented in slic.h and slic.cpp.
*/

#include <iostream>
#include <fstream>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include "DGtal/io/writers/PPMWriter.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "slic.h"


using namespace DGtal;
namespace po = boost::program_options;



struct IdColor {
  Color operator()(const unsigned int& aValue) const {
    return DGtal::Color(aValue);
  }
};



int
main(int argc, char** argv)
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



  std::string nameImage = vm["input"].as<std::string>();
  std::string nameImageOutput = vm["output"].as<std::string>();
  Image2D image = GenericReader<Image2D>::import(nameImage);
  unsigned int imageWidth = image.domain().upperBound()[0] - image.domain().lowerBound()[0] + 1;
  unsigned int imageHeight = image.domain().upperBound()[1] - image.domain().lowerBound()[1] + 1;


  int nr = vm["n"].as<int>();
  int nc = vm["w"].as<int>();
  double step = sqrt((imageWidth * imageHeight) / (double) nr);


  /* Perform the SLIC superpixel algorithm. */
  Slic slic;
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
  PPMWriter<Image2D, IdColor >::exportPPM(nameImageOutput, image, id);


  return 1;
}
