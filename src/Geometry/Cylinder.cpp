//
// Created by Sofia Iannicelli on 1/22/25.
//

#include "Cylinder.h"

Cylinder::Cylinder(Eigen::Vector3d capCenter1, Eigen::Vector3d capCenter2, double radius,
                   const std::shared_ptr<BRDF>& brdf, const std::string& description) :
        Geometry(brdf, description),
        capCenter0(capCenter1), capCenter1(capCenter2), radius(radius) {
    cylinderD = (capCenter2 - capCenter1).normalized();

    Eigen::Vector3d other = Eigen::Vector3d(0, 0, 1);
    if(cylinderD.cross(other).norm() == 0) {
        other = {0, 1, 0};
    }

    capTangent = (cylinderD.cross(other)).normalized();
}

bool Cylinder::cylinderPtInBounds(const Eigen::Vector3d& pt) {
    return cylinderD.dot(pt - capCenter0) > 0 && cylinderD.dot(pt - capCenter1) < 0;
}

double Cylinder::calculateDistToOrigin(const Eigen::Vector3d& pt) {
    return cylinderD.dot(-pt); // FIXME
}

bool Cylinder::capPtInBounds(const Eigen::Vector3d& capCenter, const Eigen::Vector3d& pt) const {
    Eigen::Vector3d diff = pt - capCenter;
    return diff.dot(diff) < radius * radius;
}

double Cylinder::findRayGeoIntersectionT(Ray ray) {
    const Eigen::Vector3d& rayO = ray.getOrigin();
    const Eigen::Vector3d& rayD = ray.getDirection();
    Eigen::Vector3d cylinderPt = capCenter0;
    double t = -1;

    // Intersect with infinite cylinder
    Eigen::Vector3d aHelper = rayD - rayD.dot(cylinderD) * cylinderD;
    double a = aHelper.dot(aHelper);
    if(a != 0) {
        Eigen::Vector3d deltaP = rayO - cylinderPt;
        double b = 2 * (rayD - rayD.dot(cylinderD) * cylinderD).dot
                (deltaP - deltaP.dot(cylinderD) * cylinderD);
        Eigen::Vector3d cHelper = deltaP - deltaP.dot(cylinderD) * cylinderD;
        double c = cHelper.dot(cHelper) - radius * radius;

        // compute the discriminant of the quadratic formula
        double discriminant = b * b - 4 * a * c;
        // if the discriminant is positive, there are one or two intersection points
        if(discriminant >= 0) {
            double sqrtDisc = sqrt(discriminant);
            // calculate smaller intersection parameter
            double t0 = (-b - sqrtDisc) / (2 * a);
            // check in bounds
            bool t0InBounds = cylinderPtInBounds(ray.getPointOnRay(t0));
            // if positive and in bounds, update t
            if(t0 >= 0 && t0InBounds) {
                t = t0;
            }
            else {
                // calculate larger t-value
                double t1 = (-b + sqrtDisc) / (2 * a);
                // check in bounds
                bool t1InBounds = cylinderPtInBounds(ray.getPointOnRay(t1));
                // if positive and in bounds, update t
                if(t1 >= 0 && t1InBounds) {
                    t = t1;
                }
            }
        }
    }

    // Intersect with each cap plane
    double d;
    double denominator;
    denominator = cylinderD.dot(ray.getDirection());
    // Cap plane 1
    // cap plane normal is the cylinder direction
    d = calculateDistToOrigin(capCenter0);
    if(denominator != 0) {
        double t3 = -(cylinderD.dot(ray.getOrigin()) + d) / denominator;
        // double t3 = -(dot(cylinderD, ray.getOrigin() - capCenter0)) / denominator;
        // check if intersection is inside cap
        bool inBounds = capPtInBounds(capCenter0, ray.getPointOnRay(t3));
        if(t3 >= 0 && inBounds && (t == -1 || t3 < t)) {
            t = t3;
        }
    }

    // Cap plane 2
    d = calculateDistToOrigin(capCenter1);
    if(denominator != 0) {
        double t4 = -(cylinderD.dot(ray.getOrigin()) + d) / denominator;
        // check if intersection is inside cap
        bool inBounds = capPtInBounds(capCenter1, ray.getPointOnRay(t4));
        if(t4 >= 0 && inBounds && (t == -1 || t4 < t)) {
            t = t4;
        }
    }

    return t;
}

RayHit Cylinder::findRayHit(Ray ray) {
    double t = findRayGeoIntersectionT(ray);

    // TODO: ?
    if(t <= 0) return RayHit::Miss();

    Eigen::Vector3d hitPoint = ray.getPointOnRay(t);
    Eigen::Vector3d hitNormal;

    if((hitPoint - capCenter0).norm() <= radius) {
        // cap 0
        hitNormal = -cylinderD;
    }
    else if((hitPoint - capCenter1).norm() <= radius) {
        // cap 1
        hitNormal = cylinderD;
    }
    else {
        // intersected cylinder
        double tAlongCenterline = (hitPoint - capCenter0).dot(cylinderD);
        Eigen::Vector3d centerlinePt = capCenter0 + tAlongCenterline * cylinderD;
        hitNormal = (hitPoint - centerlinePt).normalized();
    }

    bool backFace = hitNormal.dot(ray.getDirection()) > 0;
    if(backFace) {
        hitNormal = -hitNormal;
    }

    double u, v;
    getUV(hitPoint, u, v);

    return {true, t, u, v, hitPoint, hitNormal, backFace, brdf} ;
}

void Cylinder::getUV(const Eigen::Vector3d& point, double& u, double& v){
    Eigen::Vector3d objectSpacePoint = point - capCenter0;

    // get projection onto cap normal
    Eigen::Vector3d normalProj = objectSpacePoint.dot(cylinderD) * cylinderD;
    // get v from height
    v = normalProj.norm() / (capCenter1 - capCenter0).norm();

    // cap projection
    Eigen::Vector3d capProj = objectSpacePoint - normalProj;

    // find angle with cap tangent
    double angle = std::acos(capProj.dot(capTangent) / capProj.norm());
    u = angle / (2 * M_PI);
}

Extent Cylinder::findExtent() {
    // Extent extent = Extent(capCenter0);
    Extent extent{capCenter0};
    // https://www.gamedev.net/forums/topic/338522-bounding-box-for-a-cylinder/#:~:text=Bounds%20in%20direction%20X%20(same%20for%20Y%20and%20Z)%20can%20be%20found%20as%3A%0ALet%20A.X%3CB.X%20(otherwise%20swap%20points)%0AGood%20approximate%20lowest%20bound%20is%20A.X%2Dr%20and%20highest%20is%20B.X%2Br%20(precise%20for%20capsule).%20At%20worst%2C%20in%20one%20direction%20it%20can%20be%20larger%20than%20needed
    for(int d = 0; d < 3; ++d) {
        double a = capCenter0[d];
        double b = capCenter1[d];
        if(a > b) std::swap(a, b);
        extent.corners[0][d] = a - radius;
        extent.corners[1][d] = b + radius;
    }

    return extent;
}
