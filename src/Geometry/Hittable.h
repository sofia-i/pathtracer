//
// Created by Sofia Iannicelli on 2/4/25.
//

#ifndef RAYTRACER_2_HITTABLE_H
#define RAYTRACER_2_HITTABLE_H

#include "../Texture/Material.h"
#include "../Texture/BRDF.h"
#include "../Ray.hpp"

struct RayHit {
    RayHit() = default;
    RayHit(bool isHit, double t, double u, double v,
           const Eigen::Vector3d& intersectPt, const Eigen::Vector3d& normal,
           bool backFace, const std::shared_ptr<BRDF>& brdf) :
            isHit(isHit), t(t), u(u), v(v), point(intersectPt), normal(normal), backFace(backFace), brdf(brdf) {}

    static RayHit Miss() { return RayHit(false); }

    bool isHit;
    double t;
    double u;
    double v;
    Eigen::Vector3d point;
    Eigen::Vector3d normal;
    bool backFace;
    std::shared_ptr<BRDF> brdf;

private:
    explicit RayHit(bool isHit) : isHit(isHit) {}
};

class Hittable {
public:
    virtual RayHit findRayHit(Ray ray) = 0;
};

#endif //RAYTRACER_2_HITTABLE_H
