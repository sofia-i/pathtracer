//
//  Sphere.h
//  raytracer_2
//
//  Created by Sofia Iannicelli on 2/15/23.
//

#ifndef Sphere_h
#define Sphere_h

#include "Geometry.h"

class Sphere : public Geometry {
public:
    Sphere(Eigen::Vector3d center, double radius, const std::shared_ptr<BRDF>& brdf, const std::string& description) :
        Geometry(brdf, std::move(description)), center(center), radius(radius) {}

    ~Sphere() override = default; // I. destructor
    Sphere(const Sphere& other) = default; // II. copy constructor
    Sphere& operator=(const Sphere& other) = default; // III. copy assignment
    Sphere(Sphere&& other) noexcept = default;// IV. move constructor
    Sphere& operator=(Sphere&& other) noexcept = default; // V. move assignment

    Eigen::Vector3d getCenter() const { return center; }
    double getRadius() const { return radius; }

    RayHit findRayHit(Ray ray) override;

    Extent findExtent() override;

    void getUV(const Eigen::Vector3d& point, double& u, double& v);
    
    std::string toString() const override {
        std::string str;
        std::stringstream ss(str);
        
        ss << getDescription() << std::endl;
        ss << "\tCenter: " << center[0] << " " << center[1] << " " << center[2] << std::endl;
        ss << "\tRadius: " << getRadius() << std::endl;
        ss << Geometry::toString() << std::endl;
         
        return ss.str();
    }

private:
    double findRayGeoIntersectionT(Ray ray);

private:
    Eigen::Vector3d center;
    double radius;
};


#endif /* Sphere_h */
