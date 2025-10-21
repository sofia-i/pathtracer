//
// Created by Sofia Iannicelli on 2/4/25.
//

#ifndef RAYTRACER_2_BOUNDINGVOLUMEHIERARCHY_H
#define RAYTRACER_2_BOUNDINGVOLUMEHIERARCHY_H

#include <vector>
#include "BoundingBox.h"
#include "../Geometry/Geometry.h"

class BoundingVolumeHierarchy : public Hittable {
public:
    virtual void constructHierarchy(const std::vector<std::shared_ptr<Geometry>>& geo) = 0;
};

class MedianSplit : public BoundingVolumeHierarchy {
public:
    void constructHierarchy(const std::vector<std::shared_ptr<Geometry>>& geo) override;

    RayHit findRayHit(Ray ray) override;

private:
    void split(uint splitCount, Extent extent, const std::vector<uint>& gbbIndexes);

    static int findLongestAxis(const Extent& extent) ;
    static bool overlaps(const Extent& extent1, const Extent& extent2);

    RayHit findRayBoxesIntersection(Ray ray, double tCovered = 0.0, int count = 1);

private:
    MSBoundingBox root;
    std::vector<MSBoundingBox> bBoxes;
    std::vector<GeoBoundingBox> geoBoxes;

    // TODO: dynamically set these values
    uint maxSplitCount = 50;
    uint goalGeoCount = 100;

    const double EPSILON = 2.e-8;
};


#endif //RAYTRACER_2_BOUNDINGVOLUMEHIERARCHY_H
