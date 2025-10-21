#include <iostream>

#include "src/SceneParser.hpp"
#include "src/Pathtracer.h"
#include "src/PpmWriter.hpp"
#include "src/PngWriter.h"

int main(int argc, char *argv[]) {
    std::string usage_instructions = "Usage: <inputFilePath> <outputFilePath>";

    if(argc != 3) {
        std::cerr << usage_instructions << std::endl;
        return EXIT_FAILURE;
    }
    // first command-line argument: input_file_path
    std::string inputFilePath(argv[1]);
    // Second argument: output_file_path
    std::string outputFilePath(argv[2]);

    // parse the scene information from the input file
    SceneParser sceneParser;
    Scene scene = sceneParser.parseFile(inputFilePath);
    std::cout << scene << std::endl;

    // send the results to a ppm file for output
    std::string magicNumber = "P3";
    int numColumns = 512;
    int numRows = 512;
    int maxColorVal = 255;

    scene.camera.updateFov(numColumns, numRows);

    Pathtracer pathtracer = Pathtracer(scene);

    // do the raytracing process to get the color results for each pixel
    int*** pixelColors = pathtracer.pathtrace(numColumns, numRows);

    PngWriter::writePng(outputFilePath, numColumns, numRows, maxColorVal, pixelColors);


    // deallocate memory
    for(int i = 0; i < numRows; ++i) {
        for(int j = 0; j < numColumns; ++j) {
            delete[] pixelColors[i][j];
        }
    }
    for(int i = 0; i < numRows; ++i) {
        delete[] pixelColors[i];
    }
    delete[] pixelColors;

    return 0;
}
