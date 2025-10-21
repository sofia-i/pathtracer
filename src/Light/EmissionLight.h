//
// Created by Sofia Iannicelli on 3/6/25.
//

#ifndef PATHTRACER_EMISSIONLIGHT_H
#define PATHTRACER_EMISSIONLIGHT_H

#include "Light.h"
#include "../Texture/Texture.h"
#include "../Geometry/Hittable.h"

class EmissionLight : public Light {
public:
    EmissionLight(const std::shared_ptr<Hittable>& geo, const std::shared_ptr<Texture>& emissionTexture) :
        geo(geo), emissionTexture(emissionTexture) {}

    LightHit getLightHit(const Eigen::Vector3d &point, const Eigen::Vector3d &pointNormal) override {
        double x = pointNormal.x();
        double y = pointNormal.y();
        if(x > 0 && y < 0.1) {
            // std::cerr << "stop here" << std::endl;
        }
        // FIXME: looking for the normal connecting with the geo isn't right
        RayHit lightHit = geo->findRayHit(Ray(point, pointNormal));
        if(lightHit.isHit) {
            Eigen::Vector3d toLight = lightHit.point - point;
            return {toLight.normalized(),
                    emissionTexture->getIntensity(lightHit.u, lightHit.v, lightHit.point),
                    toLight.norm()};
        }
        else {
            return {false};
        }
    }

    std::string toString() const override {
        std::stringstream ss;
        ss << "Emmissive Light" << std::endl;
        ss << "\t" << "texture: " << emissionTexture << std::endl;
        ss << "\t" << "geo: " << geo << std::endl;
        return ss.str();
    }

private:
    std::shared_ptr<Hittable> geo;
    std::shared_ptr<Texture> emissionTexture;
};

#endif //PATHTRACER_EMISSIONLIGHT_H
