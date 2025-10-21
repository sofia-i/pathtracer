//
// Created by Sofia Iannicelli on 3/6/25.
//

#ifndef PATHTRACER_POINTLIGHT_H
#define PATHTRACER_POINTLIGHT_H

#include "Light.h"

class PointLight : public Light {
public:
    PointLight(const Eigen::Vector3d& lightColor, const Eigen::Vector3d& position) : Light(),
                                                                                     color(lightColor), position(position) { }

    ~PointLight() override = default;
    PointLight(const PointLight& other) = default;  // II. copy constructor
    PointLight(PointLight&& other) noexcept = default;  // IV. move constructor

    LightHit getLightHit(const Eigen::Vector3d &point, const Eigen::Vector3d &pointNormal) override {
        Eigen::Vector3d vectorToLight = position - point;
        return {vectorToLight.normalized(), color, vectorToLight.norm()};
    }

#if PATHTOLIGHT
    void getPathToLight(const Eigen::Vector3d& point, int index, bool& hit, Eigen::Vector3d& direction,
                        double& distance) override {
        Eigen::Vector3d vectorToLight = position - point;
        hit = true;
        distance = vectorToLight.norm();
        direction = vectorToLight.normalized();
    }
#endif

    /*
    Eigen::Vector3d getDirectionToLight(Eigen::Vector3d point, Eigen::Vector3d hitNormal) override {
        return (position - point).normalized();
    }


    void getPathToLight(const Eigen::Vector3d &point, const Eigen::Vector3d& pointNormal,
                        bool &hit, Eigen::Vector3d &direction, double &distance) override {
        Eigen::Vector3d vectorToLight = position - point;
        hit = true;
        distance = vectorToLight.norm();
        direction = vectorToLight.normalized();
    }


    double getDistToLight(const Eigen::Vector3d &point) override {
        return (position - point).norm();
    }


    Eigen::Vector3d getLightColor(Ray ray) const override {
        return color;
    }
     */

    std::string toString() const override {
        std::stringstream ss;
        ss << "Point Light" << std::endl;
        ss << "\t" << "color: " << color << std::endl;
        ss << "\t" << "Position: " << position << std::endl;
        return ss.str();
    }

private:
    Eigen::Vector3d position;
    Eigen::Vector3d color;

};

#endif //PATHTRACER_POINTLIGHT_H
