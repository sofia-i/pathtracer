//
// Created by Sofia Iannicelli on 3/6/25.
//

#ifndef PATHTRACER_AREALIGHT_H
#define PATHTRACER_AREALIGHT_H

#include "Light.h"

class AreaLight : public Light {
public:
    AreaLight(Eigen::Vector3d color,
              Eigen::Vector3d center, Eigen::Vector3d aim, Eigen::Vector3d up,
              double width, double height, double sampleResolution) :
            Light(int(width/sampleResolution) * int(height/sampleResolution)),
            color(color),
            center(center),
            width(width),
            height(height),
            wSampleCount(int(width/sampleResolution)),
            hSampleCount(int(height/sampleResolution)) {
        normal = (aim - center).normalized();
        uAxis = (up.normalized().cross(this->normal)).normalized();
        vAxis = this->normal.cross(uAxis);
        blPos = center - (width / 2.) * uAxis - (height / 2. * vAxis);

        uInc = width/wSampleCount;
        vInc = height/hSampleCount;

        dist = rn.get_dist(-0.5, 0.5);
        // dist2 = std::uniform_real_distribution<double>(0., 1.);
        jitterAmt = 0.9;
    }

    ~AreaLight() override = default;
    AreaLight(const AreaLight& other) = default;  // II. copy constructor
    AreaLight(AreaLight&& other) noexcept = default;  // IV. move constructor

    /*
    Eigen::Vector3d getDirectionToLight(Eigen::Vector3d point, Eigen::Vector3d hitNormal) override {
        return (center - point).normalized();
    }
     */

#if PATHTOLIGHT
    void getPathToLight(const Eigen::Vector3d& point, int index, bool& hit, Eigen::Vector3d& direction,
                        double& distance) override {
        int vIdx = index / hSampleCount;
        int uIdx = index % hSampleCount;
        Eigen::Vector3d lightPt = blPos + (uIdx * uInc * uAxis) + (vIdx * vInc * vAxis);
        // lightPt += jitterAmt * std::min(uInc, vInc) * Eigen::Vector3d::getRandom(dist);
        lightPt += jitterAmt * std::min(uInc, vInc) * rn.get_random_vec(dist);
        Eigen::Vector3d vectorToLight = lightPt - point;

        // check back face
        if(normal.dot(vectorToLight) > 0) {
            hit = false;
            return;
        }

        hit = true;
        distance = vectorToLight.norm();
        direction = vectorToLight.normalized();
    }
#endif

    LightHit getLightHit(const Eigen::Vector3d &point, const Eigen::Vector3d &pointNormal) override {
        double randU = rn.get_random_double_from_dist(dist);
        double randV = rn.get_random_double_from_dist(dist);

        Eigen::Vector3d target = center;
        target += randU * width * uAxis;
        target += randV * height * vAxis;

        Eigen::Vector3d vectorToLight = target - point;

        // check back face
        if(normal.dot(vectorToLight) > 0) {
            return {false};
        }

        return {vectorToLight.normalized(), color, vectorToLight.norm()};
    }

    /*
    void getPathToLight(const Eigen::Vector3d &point, const Eigen::Vector3d& pointNormal,
                        bool &hit, Eigen::Vector3d &direction, double &distance) override {
        double randU = rn.get_random_double_from_dist(dist);
        double randV = rn.get_random_double_from_dist(dist);

        Eigen::Vector3d target = center;
        target += randU * width * uAxis;
        target += randV * height * vAxis;

        Eigen::Vector3d vectorToLight = target - point;

        // check back face
        if(normal.dot(vectorToLight) > 0) {
            hit = false;
            return;
        }

        hit = true;
        distance = vectorToLight.norm();
        direction = vectorToLight.normalized();
    }
     */

    /*
    double getDistToLight(const Eigen::Vector3d &point) override {
        // TODO: change based on index?
        return (center - point).norm();
    }
     */

    /*
    Eigen::Vector3d getLightColor(Ray ray) const override {
        return color;
    }
     */

    std::string toString() const override {
        std::stringstream ss;
        ss << "Area Light" << std::endl;
        ss << "\t" << "Color: " << color[0] << " " << color[1] << " " << color[2] << std::endl;
        ss << "\t" << "Center: " << center[0] << " " << center[1] << " " << center[2] << std::endl;
        return ss.str();
    }

private:
    Eigen::Vector3d color;

    Eigen::Vector3d center;

    Eigen::Vector3d uAxis;
    Eigen::Vector3d vAxis;
    Eigen::Vector3d normal;

    Eigen::Vector3d blPos;  // bottom left position

    int wSampleCount;
    int hSampleCount;

    double width;
    double height;

    double uInc;
    double vInc;

    RandomNumber rn;
    std::uniform_real_distribution<double> dist;
    // std::uniform_real_distribution<double> dist2;
    double jitterAmt;
};

#endif //PATHTRACER_AREALIGHT_H
