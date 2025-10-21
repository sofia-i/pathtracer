//
// Created by Sofia Iannicelli on 1/17/25.
//

#include "Mesh.h"

Mesh::Mesh(std::vector<Eigen::Vector3d> vertices, std::vector<Face> faces, const std::shared_ptr<BRDF>& brdf,
           std::string description) :
                Geometry(brdf, std::move(description)),
                vertices(std::move(vertices)), faces(std::move(faces)) { }

double Mesh::findRayGeoIntersectionT(Ray ray) {
    // TODO
    return 0;
}

RayHit Mesh::findRayHit(Ray ray) {
    // TODO
    double t = findRayGeoIntersectionT(ray);
}


