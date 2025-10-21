//
//  Triangle.cpp
//  raytracer_2
//
//  Created by Sofia Iannicelli on 3/4/23.
//

#include "Triangle.hpp"
#include "../Ray.hpp"

Triangle::Triangle(const std::vector<Eigen::Vector3d>& vertices, const std::vector<Eigen::Vector2d>& uvs,
         std::shared_ptr<BRDF> brdf, const std::string& description) :
            Geometry(brdf, std::move(description)), vertices(vertices), uvs(uvs) {
    calculatePlaneNormal();
    distToOrigin = calculateDistToOrigin();
}

Triangle::Triangle(const std::vector<Eigen::Vector3d>& vertices, std::shared_ptr<BRDF> brdf, const std::string& description) :
        Triangle(vertices,
                 {{1, 0}, {0, 1}, {0, 0}},
                 brdf,
                 description) {}

double Triangle::calculateDistToOrigin() {
    return planeNormal.dot(-vertices[0]);
}

void Triangle::calculatePlaneNormal() {
    // assuming the vertices are specified in CCW order
    Eigen::Vector3d vector1 = (vertices.at(0) - vertices.at(1)).normalized();
    Eigen::Vector3d vector2 = (vertices.at(2) - vertices.at(1)).normalized();
    
    Eigen::Vector3d normal = vector2.cross(vector1);
    this->planeNormal = normal.normalized();
}

double Triangle::findRayGeoIntersectionT(Ray ray) {
    /*
     * check if the ray intersects the plane containing the triangle
     */

    double d = calculateDistToOrigin();

    // extract out the ray info
    Eigen::Vector3d ray_o = ray.getOrigin();
    Eigen::Vector3d ray_d = ray.getDirection();

    double denominator = planeNormal.dot(ray_d);
    if(denominator == 0) { return -1; }
    double t = -(planeNormal.dot(ray_o) + d) / denominator;

    if(t < 0) {
        return t;
    }

    Eigen::Vector3d intersectionPt = ray.getPointOnRay(t);
    // check if the intersection point is inside the triangle
    for(int i = 0; i < 3; ++i) {
        Eigen::Vector3d testVertex = vertices.at(i);
        Eigen::Vector3d nextVertex = vertices.at((i + 1) % 3);

        Eigen::Vector3d testNormal = (nextVertex - testVertex).cross(intersectionPt - testVertex);

        if(planeNormal.dot(testNormal) < 0) {
            return -1;
        }
    }

    return t;
}

RayHit Triangle::findRayHit(Ray ray) {
    double t = findRayGeoIntersectionT(ray);
    // miss if t is negative
    // FIXME: negative here breaks for triangle on bounding box
    if(t < 0) return RayHit::Miss();

    // Intersected, so calculate intersect information
    Eigen::Vector3d hitNormal;
    bool backFace = false;
    if(ray.getDirection().dot(planeNormal) < 0.0) {
        // The ray and the plane normal are facing in opposite directions (front face)
        hitNormal = planeNormal;
    }
    else {
        // back face
        hitNormal = -planeNormal;
        if(!doubleSided) {
            backFace = true;
        }
    }

    Eigen::Vector3d hitPoint = ray.getPointOnRay(t);

    double u, v;
    getUV(hitPoint, u, v);

    return {true, t, u, v, hitPoint, hitNormal, backFace, brdf};
}

Extent Triangle::findExtent() {
    Extent extent(vertices.at(0));
    for(auto vert : vertices) {
        extent.update(vert);
    }
    return extent;
}

void Triangle::getUV(const Eigen::Vector3d &point, double &u, double &v) {
    // use barycentric coordinates
    // https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
    Eigen::Vector3d v0 = vertices[1] - vertices[0];
    Eigen::Vector3d v1 = vertices[2] - vertices[0];
    Eigen::Vector3d v2 = point - vertices[0];

    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);
    double denom = d00 * d11 - d01 * d01;

    double w1 = (d11 * d20 - d01 * d21) / denom;
    double w2 = (d00 * d21 - d01 * d20) / denom;
    double w0 = 1.0 - w2 - w1;

    u = w0 * uvs[0][0] + w1 * uvs[1][0] + w2 * uvs[2][0];
    v = w0 * uvs[0][1] + w1 * uvs[1][1] + w2 * uvs[2][1];
}

