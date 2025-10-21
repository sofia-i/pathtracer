//
//  Scene.hpp
//  raytracer_2
//
//  Created by Sofia Iannicelli on 2/15/23.
//

#ifndef Scene_hpp
#define Scene_hpp

#include <cstdio>
#include <vector>
#include "Camera.h"
#include "Geometry/Geometry.h"
#include "BoundingVolume/BoundingBox.h"
#include "BoundingVolume/BoundingVolumeHierarchy.h"
#include "Light/Light.h"
#include "MathUtils/Utils.h"

class Scene {
public:
    Scene(Camera &camera, Eigen::Vector3d ambient_light, std::shared_ptr<Texture> background)  :
            camera(camera), ambient_light(ambient_light), background(background),
            ambient_ior(IndexOfRefraction::AIR) { }

    void process() {
        bvh = std::make_shared<MedianSplit>();
        bvh->constructHierarchy(geo);
    }

    Camera getCamera() { return camera; }
    Eigen::Vector3d getAmbientLight() const { return ambient_light; }
    double getAmbientIor() const { return ambient_ior; }

    Eigen::Vector3d getBackgroundIntensity(Ray ray) const {
        // point to origin, unit vector
        Eigen::Vector3d d = ray.getDirection();

        Eigen::Vector2d uv;
        uv[0] = 0.5 + atan2(d[2], d[0]) / (2 * M_PI);
        uv[1] = 0.5 + asin(d[1]) / M_PI;
        return background->getIntensity(uv[0], uv[1], Eigen::Vector3d::Zero());
    }

    std::vector<std::shared_ptr<Geometry>> geo;
    std::vector<std::shared_ptr<Light>> lights;
    std::shared_ptr<BoundingVolumeHierarchy> bvh;

    friend std::ostream& operator<<(std::ostream& os, Scene const &scene) {
        os << "Scene: " << std::endl;
        os << "\t" << scene.camera << std::endl;
        os << "\tAmbient Light: " << scene.ambient_light[0] << " " <<
            scene.ambient_light[1] << " " << scene.ambient_light[2] << std::endl;
        os << "\tBackground: " << *(scene.background) << std::endl;
        os << std::endl;
        for(auto&& light: scene.lights) {
            os << "\t" << *light << std::endl;
        }
        for(auto&& g : scene.geo) {
            os << "\tGeo: " << g->toString() << std::endl;
        }
        return os;
    }

    Camera camera;

private:
    Eigen::Vector3d ambient_light;
    std::shared_ptr<Texture> background;
    // Eigen::Vector3d backgroundColor;

    double ambient_ior;  // index of refraction

};

#endif /* Scene_hpp */
