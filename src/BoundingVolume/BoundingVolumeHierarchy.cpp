//
// Created by Sofia Iannicelli on 2/4/25.
//

#include "BoundingVolumeHierarchy.h"
#include <numeric>
#include <cassert>
#include <cmath>

void MedianSplit::constructHierarchy(const std::vector<std::shared_ptr<Geometry>>& geo) {
    if(geo.empty()) return;
    // Surround each geo with a bounding box
    geoBoxes.reserve(geo.size());
    for(const auto & g : geo) {
        geoBoxes.push_back(GeoBoundingBox(g));
    }

    // surround entire scene with bounding box
    Extent extent = geoBoxes[0].getBounds();
    for(const auto & geoBox : geoBoxes) {
        extent.update(geoBox.getBounds());
    }
    root = MSBoundingBox(extent);

    uint iterCount = 0;

    std::vector<uint> gbbIdxs;
    gbbIdxs.resize(geo.size());
    std::iota(gbbIdxs.begin(), gbbIdxs.end(), 0);

    split(++iterCount, extent, gbbIdxs);
}

void MedianSplit::split(uint splitCount, Extent extent, const std::vector<uint>& gbbIndexes) {
    // check if max number of splits has been reached
    if(splitCount > maxSplitCount || gbbIndexes.size() <= goalGeoCount) {
        bBoxes.push_back(MSBoundingBox(extent, gbbIndexes));
        return;
    }

    // find largest magnitude extent to split across
    int largestAxis = findLongestAxis(extent);

    double median = extent.corners[0][largestAxis] +
            0.5 * (extent.corners[1][largestAxis] - extent.corners[0][largestAxis]);
    Extent leftExt = extent;
    Extent rightExt = extent;

    leftExt.corners[1][largestAxis] = median;
    rightExt.corners[0][largestAxis] = median;
    std::vector<uint> lGbbIndexes;
    std::vector<uint> rGbbIndexes;
    for(uint idx : gbbIndexes) {
        if(overlaps(leftExt, geoBoxes[idx].getBounds())) {
            lGbbIndexes.push_back(idx);
        }
        if(overlaps(rightExt, geoBoxes[idx].getBounds())) {
            rGbbIndexes.push_back(idx);
        }
    }

    split(++splitCount, leftExt, lGbbIndexes);
    split(++splitCount, rightExt, rGbbIndexes);
}

bool MedianSplit::overlaps(const Extent& extent1, const Extent& extent2) {
    // projections onto each axis must all overlap
    // https://stackoverflow.com/questions/20925818/algorithm-to-check-if-two-boxes-overlap
    for(int d = 0; d < 3; ++d) {
        if(extent1.corners[1][d] < extent2.corners[0][d]) return false;
        if(extent2.corners[1][d] < extent1.corners[0][d]) return false;
    }
    return true;
}

RayHit MedianSplit::findRayHit(Ray ray) {
    if(!root.findRayBoxHit(ray).isHit) {
        // doesn't intersect any geo
        return RayHit::Miss();
    }

    return findRayBoxesIntersection(ray);
}

RayHit MedianSplit::findRayBoxesIntersection(Ray ray, double tCovered, int count) {
    // FIXME: optimization with excluding already checked?
    // Find closest intersected box
    RayAABBHit closestBoxHit = RayAABBHit::Miss();
    uint intersectBoxIdx = -1;
    for(int i = 0; i < bBoxes.size(); ++i) {
        RayAABBHit bBoxHit = bBoxes[i].findRayBoxHit(ray);
        if(bBoxHit.isHit && (!closestBoxHit.isHit || bBoxHit.tEnter < closestBoxHit.tEnter)) {
            closestBoxHit = bBoxHit;
            intersectBoxIdx = i;
        }
    }

    // return if no box intersected
    if(!closestBoxHit.isHit) return RayHit::Miss();

    // Box found: Intersect with each geo in the box (between enter and exit points)
    double enterT = std::max(0., closestBoxHit.tNear);
    double exitT = closestBoxHit.tFar;
    assert(exitT >= 0);


    Ray geoRay = ray;

    RayHit closestRayHit = RayHit::Miss();
    for(uint idx : bBoxes[intersectBoxIdx].getGbbIndexes()) {
        RayHit hit = geoBoxes[idx].findRayHit(geoRay);
        if(hit.isHit && (hit.t - exitT < EPSILON) && (!closestRayHit.isHit || hit.t < closestRayHit.t)) {
            closestRayHit = hit;
        }
    }

    // hit real geo!
    if(closestRayHit.isHit) {
        closestRayHit.t += tCovered;
        return closestRayHit;
    }

    // No hit - continue to the next box
    // find exit point
    Eigen::Vector3d exitPoint = ray.getPointOnRay(exitT) + (EPSILON * ray.getDirection());
    Ray nextRay = Ray(exitPoint, ray.getDirection());
    return findRayBoxesIntersection(nextRay, tCovered + exitT, ++count);
}

int MedianSplit::findLongestAxis(const Extent &extent) {
    // compare the first two axes
    int longestAxis = extent.lengthX() > extent.lengthY() ? 0 : 1;
    double length = extent.lengthX() > extent.lengthY() ? extent.lengthX() : extent.lengthY();

    // compare winner against z axis
    return length > extent.lengthZ() ? longestAxis : 2;
}

