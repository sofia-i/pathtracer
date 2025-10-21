//
// Created by Sofia Iannicelli on 1/22/25.
//

#ifndef RAYTRACER_2_CYLINDER_H
#define RAYTRACER_2_CYLINDER_H

#include "Geometry.h"

class Cylinder : public Geometry {
public:
    Cylinder(Eigen::Vector3d capCenter1, Eigen::Vector3d capCenter2, double radius,
             const std::shared_ptr<BRDF>& brdf, const std::string& description);

    ~Cylinder() override = default; // I. destructor
    Cylinder(const Cylinder& other) = default; // II. copy constructor
    Cylinder& operator=(const Cylinder& other) = default; // III. copy assignment
    Cylinder(Cylinder&& other) noexcept = default;// IV. move constructor
    Cylinder& operator=(Cylinder&& other) noexcept = default; // V. move assignment

    Eigen::Vector3d getCapCenter1() const { return capCenter0; }
    Eigen::Vector3d getCapCenter2() const { return capCenter1; }
    Eigen::Vector3d getDir() const { return cylinderD; }
    double getRadius() const { return radius; }

    RayHit findRayHit(Ray ray) override;

    Extent findExtent() override;

    void getUV(const Eigen::Vector3d& point, double& u, double& v);

    std::string toString() const override {
        std::string str;
        std::stringstream ss(str);

        ss << getDescription() << std::endl;
        ss << "\tCap Centers: " << capCenter0 << "; " << capCenter1 << std::endl;
        ss << "\tRadius: " << getRadius() << std::endl;
        ss << Geometry::toString() << std::endl;

        return ss.str();
    }

private:
    Eigen::Vector3d capCenter0;
    Eigen::Vector3d capCenter1;
    double radius;

    Eigen::Vector3d cylinderD;
    Eigen::Vector3d capTangent;

    double findRayGeoIntersectionT(Ray ray);
    bool cylinderPtInBounds(const Eigen::Vector3d& pt);
    bool capPtInBounds(const Eigen::Vector3d& capCenter, const Eigen::Vector3d& pt) const;
    double calculateDistToOrigin(const Eigen::Vector3d& pt);
};


#endif //RAYTRACER_2_CYLINDER_H
