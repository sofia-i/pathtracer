//
//  Raytracer.hpp
//  raytracer_2
//
//  Created by Sofia Iannicelli on 2/15/23.
//

#ifndef Raytracer_hpp
#define Raytracer_hpp

#include <cstdio>
#include <stack>
#include "Scene.hpp"
#include "Camera.h"
#include "Ray.hpp"
#include "Geometry/Hittable.h"

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

class Raytracer {
public:
    explicit Raytracer(Scene scene) : scene(scene), raysPerPixelPerSide(2) { }
    Raytracer(Scene scene, int raysPerPixelPerSide) : scene(scene), raysPerPixelPerSide(raysPerPixelPerSide) { }

    ~Raytracer();
    
    int*** raytrace(int numCols, int numRows);
    
    Scene* getScene() {
        return &scene;
    }
    
    Camera getCamera() {
        return scene.getCamera();
    }

private:
    const bool LOG_TIME = true;
    const int LOG_INTERVAL = 500;
    const bool USE_BOUNDING_VOLUME = true;

    const int MAX_NUM_RAYS = 5;
    const double EPSILON = 2e-8;
    const int raysPerPixelPerSide;

    // TODO: should scene be a member?
    Scene scene;

    /**
     *
     * @param target world-space coordinate to raytrace through
     * @return color tuple
     */
    Eigen::Vector3i getRayResult(Eigen::Vector3d target);
    Eigen::Vector3i getRayResult(Ray ray, int rayCount, std::stack<double>& iors);

    RayHit getClosestIntersection(const Ray& ray);
    double getInShadow(const Eigen::Vector3d& intersectPt, const std::shared_ptr<Light>& light);

    WorldSpaceCoord calculateWorldSpaceCoords(int numCols, int numRows);
    Eigen::Vector3i illuminationEq(const std::shared_ptr<Material>& mat, const double u, const double v,
                             const Eigen::Vector3d& normal, const Eigen::Vector3d& view,
                             const Eigen::Vector3d& intersectPt);

    inline Eigen::Vector3d getAmbient(const std::shared_ptr<Material>& mat, const Eigen::Vector3d& point,
                                   const double u, const double v);
    inline Eigen::Vector3d getDiffuse(const std::shared_ptr<Material>& mat, const Eigen::Vector3d& point,
                                   const double u, const double v,
                                   const std::shared_ptr<Light>& light,
                                   const Eigen::Vector3d& normal, const Eigen::Vector3d& toLight);
    inline Eigen::Vector3d getSpecular(const std::shared_ptr<Material>& mat, const std::shared_ptr<Light>& light,
                                    const Eigen::Vector3d& normal, const Eigen::Vector3d& toLight,
                                    const Eigen::Vector3d& view);
    inline Eigen::Vector3i getTransmission(const std::shared_ptr<Material>& mat, const Eigen::Vector3d& normal,
                                     const Eigen::Vector3d& rayD, const Eigen::Vector3d& intersectPt,
                                     double iorRatio, int rayCount, std::stack<double>& iors);
    inline Ray getTransmissionRay(const std::shared_ptr<Material>& mat, const Eigen::Vector3d& normal,
                                  const Eigen::Vector3d& rayD, const Eigen::Vector3d& intersectPt, double iorRatio) const;
    inline Eigen::Vector3i getReflection(const std::shared_ptr<Material>& mat, const Eigen::Vector3d& normal,
                                   const Eigen::Vector3d& toView, const Eigen::Vector3d& intersectPt, int rayCount,
                                   std::stack<double>& iors);

    void getIorAcrossIntersection(const std::shared_ptr<Material>& mat, bool isBackFace,
                                  double& iorIn, double& iorOut, double& iorRatio,
                                  std::stack<double>& iors);
    double getPortionReflected(const Eigen::Vector3d& normal, const Eigen::Vector3d& rayD, const double matRefl,
                               const double& iorIn, const double& iorOut, const double& iorRatio);

    void showProgress(int index, int total);

};

#endif /* Raytracer_hpp */
