//
// Created by Sofia Iannicelli on 1/17/25.
//

#ifndef RAYTRACER_2_MESH_H
#define RAYTRACER_2_MESH_H

#include "Geometry.h"
#include <utility>
#include <vector>

struct Face {
    int vIdx[4];
};

class Mesh : public Geometry {
public:
    Mesh(std::vector<Eigen::Vector3d> vertices, std::vector<Face> faces, const std::shared_ptr<BRDF>& brdf,
         std::string description);

    ~Mesh() override = default; // I. destructor
    Mesh(const Mesh& other) = default; // II. copy constructor
    Mesh& operator=(const Mesh& other) = default; // III. copy assignment
    Mesh(Mesh&& other) noexcept = default;// IV. move constructor
    Mesh& operator=(Mesh&& other) noexcept = default; // V. move assignment

    RayHit findRayHit(Ray ray) override;

private:
    double findRayGeoIntersectionT(Ray ray);

private:
    std::vector<Eigen::Vector3d> vertices;
    std::vector<Face> faces;

};


#endif //RAYTRACER_2_MESH_H
