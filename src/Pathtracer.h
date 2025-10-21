//
// Created by Sofia Iannicelli on 3/3/25.
//

#ifndef RAYTRACER_2_PATHTRACER_H
#define RAYTRACER_2_PATHTRACER_H

#include "Scene.hpp"

struct WorldSpaceCoord {
    double maxU;
    double maxV;

    double uInc;
    double vInc;

    Eigen::Vector3d uAxis;
    Eigen::Vector3d vAxis;

    WorldSpaceCoord(double maxU, double maxV, double uInc, double vInc,
                    Eigen::Vector3d uAxis, Eigen::Vector3d vAxis) :
            maxU(maxU), maxV(maxV), uInc(uInc), vInc(vInc),
            uAxis(uAxis), vAxis(vAxis) {}
};

class Pathtracer {
public:
    explicit Pathtracer(Scene scene) : scene(scene), maxDepth(5) {
        // iors.reserve(3);
    }
    Pathtracer(Scene scene, int maxDepth) : scene(scene), maxDepth(maxDepth) {}

    int*** pathtrace(int numCols, int numRows);
    WorldSpaceCoord calculateWorldSpaceCoords(int numCols, int numRows);

private:
    Eigen::Vector3d getPathLi(const Eigen::Vector3d& target);
    Eigen::Vector3d getPathLi(Ray ray, int bounceCount, std::vector<double> iors);

    RayHit getClosestIntersection(const Ray& ray);
    double getInShadow(const Eigen::Vector3d& point, const Eigen::Vector3d& direction,
                       const std::shared_ptr<Light>& light, const double& distToLight);

    Scene scene;

    int raysPerPixelPerSide = 1;
    int pathsPerPixel = 100;
    int maxDepth;

    bool LOG_TIME = true;

    // std::vector<double> iors;

};


#endif //RAYTRACER_2_PATHTRACER_H
