//
// Created by Sofia Iannicelli on 3/6/25.
//

#ifndef PATHTRACER_DIRECTIONALLIGHT_H
#define PATHTRACER_DIRECTIONALLIGHT_H

#include "Light.h"

class DirectionalLight : public Light{
public:
    DirectionalLight(const Eigen::Vector3d& lightColor, const Eigen::Vector3d& directionToLight) :
            Light(),
            color(lightColor)
    {
        this->directionToLight = directionToLight.normalized();
    }

    ~DirectionalLight() override = default;  // I. destructor
    DirectionalLight(const DirectionalLight& other) = default;  // II. copy constructor
    DirectionalLight(DirectionalLight&& other) noexcept = default;  // IV. move constructor

    /*
    Eigen::Vector3d getDirectionToLight(Eigen::Vector3d point, Eigen::Vector3d hitNormal) override {
        return directionToLight;
    }
     */

#if PATHTOLIGHT
    void getPathToLight(const Eigen::Vector3d& point, int index, bool& hit, Eigen::Vector3d& direction,
                        double& distance) override {
        hit = true;
        direction = directionToLight;
        distance = DBL_MAX;
    }
#endif

    LightHit getLightHit(const Eigen::Vector3d &point, const Eigen::Vector3d &pointNormal) override {
        return {directionToLight, color, DBL_MAX};
    }

    /*
    void getPathToLight(const Eigen::Vector3d &point, const Eigen::Vector3d& pointNormal,
                        bool &hit, Eigen::Vector3d &direction, double &distance) override {
        hit = true;
        direction = directionToLight;
        distance = DBL_MAX;
    }


    double getDistToLight(const Eigen::Vector3d& point) override {
        return DBL_MAX;
    }


    Eigen::Vector3d getLightColor(Ray ray) const override {
        return color;
    }
     */

    std::string toString() const override {
        std::stringstream ss;
        ss << "Directional Light" << std::endl;
        ss << "\t" << "color: " << color.x() << " " << color.y() << " " << color.z() << std::endl;
        ss << "\t" << "Direction to light: " << directionToLight.x() << " " << directionToLight.y() << " "
           << directionToLight.z() << std::endl;
        return ss.str();
    }

private:
    Eigen::Vector3d directionToLight;
    Eigen::Vector3d color;

};

#endif //PATHTRACER_DIRECTIONALLIGHT_H
