//
// Created by Sofia Iannicelli on 1/15/25.
//

#ifndef RAYTRACER_2_LIGHT_H
#define RAYTRACER_2_LIGHT_H

#include <Eigen/Eigen>
#include <cfloat>
#include "../MathUtils/Utils.h"
#include "../Texture/Texture.h"
#include "../Ray.hpp"

#define PATHTOLIGHT 0

struct LightHit {
    LightHit(const Eigen::Vector3d& direction, const Eigen::Vector3d& intensity, const double& distance) :
            direction(direction), intensity(intensity), distance(distance), hit(true) {}
    /**
     * to be used for misses
     * @param hit
     */
    LightHit(bool hit) : hit(hit) {}

    Eigen::Vector3d direction;
    Eigen::Vector3d intensity;
    double distance;
    bool hit;

};

class Light {
public:
    Light() : shadowRayCount(1) {}
    explicit Light(int shadowRayCount) : shadowRayCount(shadowRayCount) {}
    // explicit Light(std::shared_ptr<Texture> texture) : texture(texture), shadowRayCount(1) {}
    // Light(std::shared_ptr<Texture> texture, int shadowRayCount) : texture(texture),
    //                                                     shadowRayCount(shadowRayCount) {}

    virtual ~Light() = default; // I. destructor
    Light(const Light& other) = default; // II. copy constructor
    Light(Light&& other) noexcept = default;// IV. move constructor

    // virtual Eigen::Vector3d getDirectionToLight(Eigen::Vector3d point, Eigen::Vector3d hitNormal) = 0;

#if PATHTOLIGHT
    virtual void getPathToLight(const Eigen::Vector3d& point, int index, bool& hit, Eigen::Vector3d& direction,
                                double& distance) = 0;
#endif

    /*
    virtual void getPathToLight(const Eigen::Vector3d& point, const Eigen::Vector3d& pointNormal,
                                bool& hit, Eigen::Vector3d& direction, double& distance) = 0;
    */

    // virtual double getDistToLight(const Eigen::Vector3d& point) = 0;

    virtual LightHit getLightHit(const Eigen::Vector3d& point, const Eigen::Vector3d& pointNormal) = 0;

    // virtual Eigen::Vector3d getLightColor(Ray ray) const = 0;

    virtual std::string toString() const = 0;

    friend std::ostream& operator<<(std::ostream& os, const Light& light) {
        os << light.toString();
        return os;
    }

    const int shadowRayCount;

private:
    // std::shared_ptr<Texture> texture;

};


#endif //RAYTRACER_2_LIGHT_H
