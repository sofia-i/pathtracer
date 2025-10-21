//
// Created by Sofia Iannicelli on 2/4/25.
//

#ifndef RAYTRACER_2_EXTENT_H
#define RAYTRACER_2_EXTENT_H

struct Extent {
    Extent() : corners() {}
    Extent(double minX, double maxX, double minY, double maxY, double minZ, double maxZ) {
        corners[0] = {minX, minY, minZ};
        corners[1] = {maxX, maxY, maxZ};
    }
    Extent(double x, double y, double z) : Extent(x, x, y, y, z, z) {}
    explicit Extent(Eigen::Vector3d v) : Extent(v.x(), v.y(), v.z()) { }

    inline double lengthX() const { return corners[1].x() - corners[0].x(); }
    inline double lengthY() const { return corners[1].y() - corners[0].y(); }
    inline double lengthZ() const { return corners[1].z() - corners[0].z(); }

    void update(Eigen::Vector3d point) {
        if(point.x() < corners[0].x()) corners[0][0] = point.x();
        if(point.x() > corners[1].x()) corners[1][0] = point.x();
        if(point.y() < corners[0].y()) corners[0][1] = point.y();
        if(point.y() > corners[1].y()) corners[1][1] = point.y();
        if(point.z() < corners[0].z()) corners[0][2] = point.z();
        if(point.z() > corners[1].z()) corners[1][2] = point.z();
    }

    void update(Extent other) {
        corners[0][0] = std::min(corners[0].x(), other.corners[0].x());
        corners[0][1] = std::min(corners[0].y(), other.corners[0].y());
        corners[0][2] = std::min(corners[0].z(), other.corners[0].z());
        corners[1][0] = std::max(corners[1].x(), other.corners[1].x());
        corners[1][1] = std::max(corners[1].y(), other.corners[1].y());
        corners[1][2] = std::max(corners[1].z(), other.corners[1].z());
    }

    bool containsPoint(const Eigen::Vector3d& point) const {
        for(int axis = 0; axis < 3; ++axis) {
            if(point[axis] < corners[0][axis]) return false;
            if(point[axis] > corners[1][axis]) return false;
        }
        return true;
    }

public:
    Eigen::Vector3d corners[2];
};

#endif //RAYTRACER_2_EXTENT_H
