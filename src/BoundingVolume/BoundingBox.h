//
// Created by Sofia Iannicelli on 1/31/25.
//

#ifndef RAYTRACER_2_BOUNDINGBOX_H
#define RAYTRACER_2_BOUNDINGBOX_H

#include <vector>
#include <memory>
#include "../Geometry/Hittable.h"
#include "../Geometry/Cylinder.h"
#include "../Geometry/Triangle.hpp"
#include "../Geometry/Sphere.h"

typedef unsigned int uint;

struct RayAABBHit {
public:
    RayAABBHit(bool isHit, double tNear, double tFar) : isHit(isHit), tNear(tNear), tFar(tFar) {
        tEnter = std::max(0., tNear);
    }
    static RayAABBHit Miss() { return RayAABBHit(false); }

    bool isHit;
    double tNear;
    double tFar;
    double tEnter;

private:
    explicit RayAABBHit(bool isHit) : isHit(isHit), tNear(-1), tFar(-1), tEnter(-1) {}
};

class AxisAlignedBoundingBox {
public:
    AxisAlignedBoundingBox() = default;
    explicit AxisAlignedBoundingBox(Extent bounds) : bounds(bounds) {}

    RayAABBHit findRayBoxHit(Ray ray);

    Extent getBounds() const { return bounds; }
protected:
    Extent bounds;
};

/**
 * Geometry axis-aligned bounding box
 */
class GeoBoundingBox : public Hittable, public AxisAlignedBoundingBox {
public:
    explicit GeoBoundingBox(Extent extent, const std::shared_ptr<Geometry>& geoPtr);
    explicit GeoBoundingBox(const std::shared_ptr<Geometry>& geoPtr);

public:
    RayHit findRayHit(Ray ray) override;

private:
    bool rayHitsBox(Ray ray);

private:
    std::shared_ptr<Geometry> geoPtr;
};

/**
 * Bounding box for use with median-split
 */
class MSBoundingBox : public AxisAlignedBoundingBox {
public:
    MSBoundingBox() = default;
    explicit MSBoundingBox(Extent bounds);
    MSBoundingBox(Extent bounds, std::vector<uint> gbbIndexes);

    std::vector<uint> getGbbIndexes() { return gbbIndexes; }

private:
    std::vector<uint> gbbIndexes;  // geo bounding box indexes
};

#endif //RAYTRACER_2_BOUNDINGBOX_H
