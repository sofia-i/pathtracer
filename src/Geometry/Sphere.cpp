//
//  Sphere.cpp
//  raytracer_2
//
//  Created by Sofia Iannicelli on 3/4/23.
//

#include "Sphere.h"
#include "../Ray.hpp"

double Sphere::findRayGeoIntersectionT(Ray ray) {
    // extract information from the ray and sphere
    Eigen::Vector3d ray_o = ray.getOrigin();
    Eigen::Vector3d ray_d = ray.getDirection();

    Eigen::Vector3d sphere_c = this->getCenter();

    // attempting to solve for intersection - compute B and C in the quadratic formula
    double b = 2 * (ray_d.x() * ray_o.x() - ray_d.x() * sphere_c.x() +
                    ray_d.y() * ray_o.y() - ray_d.y() * sphere_c.y() +
                    ray_d.z() * ray_o.z() - ray_d.z() * sphere_c.z());
    double c = pow(ray_o.x(), 2) - 2 * ray_o.x() * sphere_c.x() + pow(sphere_c.x(), 2) +
               pow(ray_o.y(), 2) - 2 * ray_o.y() * sphere_c.y() + pow(sphere_c.y(), 2) +
               pow(ray_o.z(), 2) - 2 * ray_o.z() * sphere_c.z() + pow(sphere_c.z(), 2) -
               pow(radius, 2);

    // compute the discriminant of the quadratic formula
    double discriminant = pow(b, 2) - 4 * c;
    // if the discriminant is negative, there is no intersection point
    if(discriminant < 0) {
        return -1.0;
    }
    else {
        // calculate smaller intersection parameter
        double t0 = (-b - sqrt(discriminant)) / 2;
        // if positive, update t and return
        if(t0 > 0) {
            return t0;
        }
        else {
            // calculate larger t-value
            double t1 = (-b + sqrt(discriminant)) / 2;
            return t1;
        }
    }
}

RayHit Sphere::findRayHit(Ray ray) {
    double t = findRayGeoIntersectionT(ray);

    // if not intersected, return false
    if(t <= 0) return RayHit::Miss();

    // Calculate intersection information
    Eigen::Vector3d hitPoint = ray.getOrigin() + t * ray.getDirection();
    // calculate hit normal
    Eigen::Vector3d hitNormal = (hitPoint - this->getCenter()).normalized();
    // determine back face (if ray starts inside sphere)
    bool backFace = (ray.getOrigin() - center).norm() < radius;

    double u, v;
    getUV(hitPoint, u, v);

    return {true, t, u, v, hitPoint, hitNormal, backFace, brdf};
}

Extent Sphere::findExtent() {
    return {
            center.x() - radius,
            center.x() + radius,
            center.y() - radius,
            center.y() + radius,
            center.z() - radius,
            center.z() + radius
    };
}

void Sphere::getUV(const Eigen::Vector3d &point, double &u, double &v) {
    Eigen::Vector3d unitP = (point - center) / radius;

    double theta = std::acos(-unitP.y());
    double phi = std::atan2(-unitP.z(), unitP.x()) + M_PI;

    u = phi / (2 * M_PI);
    v = theta / M_PI;
}
