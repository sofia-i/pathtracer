//
// Created by Sofia Iannicelli on 2/18/25.
//

#ifndef RAYTRACER_2_IMAGE_H
#define RAYTRACER_2_IMAGE_H

#include <string>
#include <memory>
#include <Eigen/Eigen>

class Image {
public:
    explicit Image(const std::string& filename);
    explicit Image(const char *filename);
    ~Image() {
        if(image != NULL) {
            delete[] image;
        }
    }

    /**
     * Get image object from image file
     * @param filename filepath to the image
     * @return a shared pointer to the new Image object
     */
    static std::shared_ptr<Image> getLoadedImage(const std::string& filename);

    /**
     * Get the r,g,b value of the image from texture coordinates
     * @param u Texture coordinate in range [0, 1]
     * @param v Texture coordinate in range [0, 1]
     * @return r,g,b intensity on scale [0, 1]
     */
    Eigen::Vector3d getPixelValue(double u, double v) const;
    /**
     * Get r,g,b value of image from width and height
     * @param x x pixel coordinate
     * @param y y pixel coordinate
     * @return
     */
    Eigen::Vector3d getPixelValue(int x, int y) const;

private:
    bool isLoaded;  // Whether the image contents have been loaded into memory
    int width;  // width pixel count
    int height;  // heigh pixel count
    int channels;  // number of channels loaded in
    std::string filename;
    unsigned char* image;
};


#endif //RAYTRACER_2_IMAGE_H
