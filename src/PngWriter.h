//
// Created by Sofia Iannicelli on 3/4/25.
//

#ifndef PATHTRACER_PNGWRITER_H
#define PATHTRACER_PNGWRITER_H

#include "ImageUtils.h"
#include <iostream>

class PngWriter {
public:
    static void writePng(const std::string& outFilepath,
                                  int numCols, int numRows, int maxColorVal,
                                  int*** pixelColors) {
        int nChannels = 3;

        // flatten into unsigned char array
        auto* pixelColorsC = new unsigned char[numRows * numCols * nChannels];
        for(int i = 0; i < numRows; ++i) {
            for(int j = 0; j < numCols; ++j) {
                for(int k = 0; k < nChannels; ++k) {
                    int index = i * numCols * nChannels + j * nChannels + k;
                    pixelColorsC[index] = (unsigned char) pixelColors[i][j][k];
                    // memcpy(pixelColorsC + index, (unsigned char *)(pixelColors + index),
                    //       sizeof(unsigned char));
                }
            }
        }

        ImageUtils::write_png(outFilepath.c_str(),
                              numCols,
                              numRows,
                              nChannels,
                              pixelColorsC,
                              sizeof(unsigned char) * nChannels * numCols);
    }
};

#endif //PATHTRACER_PNGWRITER_H
