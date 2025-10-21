//
//  Triangle.h
//  raytracer_2
//
//  Created by Sofia Iannicelli on 3/4/23.
//

#ifndef Triangle_hpp
#define Triangle_hpp

#include <cstdio>
#include <vector>
#include "Geometry.h"

class Triangle : public Geometry {
public:
    Triangle(const std::vector<Eigen::Vector3d>& vertices, const std::vector<Eigen::Vector2d>& uvs,
             std::shared_ptr<BRDF> brdf, const std::string& description);
    Triangle(const std::vector<Eigen::Vector3d>& vertices, std::shared_ptr<BRDF> brdf, const std::string& description);

    ~Triangle() override = default; // I. destructor
    Triangle(const Triangle& other) = default; // II. copy constructor
    Triangle& operator=(const Triangle& other) = default; // III. copy assignment
    Triangle(Triangle&& other) noexcept = default;// IV. move constructor
    Triangle& operator=(Triangle&& other) noexcept = default; // V. move assignment

    RayHit findRayHit(Ray ray) override;

    Extent findExtent() override;

    void getUV(const Eigen::Vector3d& point, double& u, double& v);

    std::vector<Eigen::Vector3d> getVertices() { return vertices; }
    
    std::string toString() const override {
        std::string str;
        std::stringstream ss(str);
        
        ss << getDescription() << std::endl;
        // print out the vertices
        ss << Geometry::toString() << std::endl;
        
        return ss.str();
    }

private:
    double findRayGeoIntersectionT(Ray ray);

private:
    void calculatePlaneNormal();
    double calculateDistToOrigin();

    Eigen::Vector3d planeNormal;
    double distToOrigin;
    bool doubleSided = true;
    std::vector<Eigen::Vector3d> vertices;
    std::vector<Eigen::Vector2d> uvs;
};


#endif /* Triangle_h */
