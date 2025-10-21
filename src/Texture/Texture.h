//
// Created by Sofia Iannicelli on 2/19/25.
//

#ifndef RAYTRACER_2_TEXTURE_H
#define RAYTRACER_2_TEXTURE_H

#include "../Image.h"

class Texture {
public:
    virtual Eigen::Vector3d getIntensity(double u, double v, const Eigen::Vector3d& point) const = 0;

    virtual std::string toString() const = 0;

    friend std::ostream& operator<<(std::ostream& os, const Texture& tex) {
        os << tex.toString();
        return os;
    }
};

class ColorTexture : public Texture {
public:
    explicit ColorTexture(Eigen::Vector3d color) : color(color) { }

    Eigen::Vector3d getIntensity(double u, double v, const Eigen::Vector3d& point) const override {
        return color;
    }

    std::string toString() const override {
        std::stringstream ss;
        ss << color[0] << " " << color[1] << " " << color[2];
        return ss.str();
    }

private:
    Eigen::Vector3d color;
};

class ImageTexture : public Texture {
public:
    explicit ImageTexture(const std::string& filename) {
        image = Image::getLoadedImage(filename);
    }

    Eigen::Vector3d getIntensity(double u, double v, const Eigen::Vector3d& point) const override {
        return image->getPixelValue(u, v);
    }

    std::string toString() const override {
        return "image texture";
    }

private:
    std::shared_ptr<Image> image;
};


#endif //RAYTRACER_2_TEXTURE_H
