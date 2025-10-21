//
//  Ray.hpp
//  raytracer_2
//
//  Created by Sofia Iannicelli on 2/16/23.
//

#ifndef Ray_hpp
#define Ray_hpp

#include <cstdio>
#include <Eigen/Eigen>

class Ray {
private:
    Eigen::Vector3d origin;
    Eigen::Vector3d direction;
    Eigen::Vector3d invDir;
    
public:
    Ray(Eigen::Vector3d origin, Eigen::Vector3d direction) {
        this->origin = origin;
        this->direction = direction.normalized();
        this->invDir = this->direction.cwiseInverse();
    }

    Eigen::Vector3d getDirection() { return direction; }
    Eigen::Vector3d getOrigin() { return origin; }
    Eigen::Vector3d getInvDir() { return invDir; }
    Eigen::Vector3d getPointOnRay(double t) const {
        return origin + t * direction;
    }
    
};

#endif /* Ray_hpp */
