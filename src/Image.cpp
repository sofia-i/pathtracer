//
// Created by Sofia Iannicelli on 2/18/25.
//

#include "Image.h"
#include "ImageUtils.h"
#include <cassert>

Image::Image(const std::string& filename) : filename(filename), isLoaded(false) { }

Image::Image(const char *filename) : filename(filename), isLoaded(false) { }

std::shared_ptr<Image> Image::getLoadedImage(const std::string& filename) {
    std::shared_ptr<Image> img = std::make_shared<Image>(filename);
    img->image = ImageUtils::load_img(filename.c_str(), &img->width, &img->height, &img->channels, 3);
    img->isLoaded = true;
    return img;
}

Eigen::Vector3d Image::getPixelValue(double u, double v) const {
    if(!isLoaded) {
        throw std::runtime_error("Requested pixel from non-loaded Image.");
    }
    int x = std::min(int(this->width * u), this->width - 1);
    int y = std::min(int(this->height * (1 - v)), this->height - 1);
    return getPixelValue(x, y);
}

Eigen::Vector3d Image::getPixelValue(int x, int y) const {
    if(!isLoaded) {
        throw std::runtime_error("Requested pixel from non-loaded Image.");
    }
    assert(x < this->width);
    assert(y < this->height);

    unsigned char *p = image + (3 * (y * this->width + x));
    unsigned char r = p[0];
    unsigned char g = p[1];
    unsigned char b = p[2];
    double rd = double(r) / 256;
    double gd = double(g) / 256;
    double bd = double(b) / 256;
    return {rd, gd, bd};
}


