//
//  PpmWriter.hpp
//  raytracer_2
//
//  Created by Sofia Iannicelli on 2/15/23.
//

#ifndef PpmWriter_hpp
#define PpmWriter_hpp

#include <cstdio>
#include <fstream>

class PpmWriter {
public:
    static std::ofstream writePpm(const std::string& outFilepath, const std::string& magicNum,
                                  int numCols, int numRows, int maxColorVal, int*** pixelColors) {
        std::ofstream outfile = setUpFile(outFilepath);
        writeHeaders(outfile, magicNum, numCols, numRows, maxColorVal);

        // write the pixel values
        for(int i = 0; i < numRows; ++i) {
            for(int j = 0; j < numCols; ++j) {
                outfile << pixelColors[i][j][0] << " ";
                outfile << pixelColors[i][j][1] << " ";
                outfile << pixelColors[i][j][2] << " ";
            }
        }

        return outfile;
    }

private:
    static std::ofstream setUpFile(const std::string& outFilepath) {
        // opening file to write output to
        std::ofstream outfile;
        outfile.open(outFilepath);
        return outfile;
    }

    static void writeHeaders(std::ofstream& outfile, const std::string& magicNum,
                              int numCols, int numRows, int maxColorVal) {
        // write file headers
        outfile << magicNum << std::endl; // magic number
        outfile << numCols << " " << numRows << std::endl; // set num cols and rows
        outfile << maxColorVal << std::endl; // maximum color value
    }

};

#endif /* PpmWriter_hpp */
