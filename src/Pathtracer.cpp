//
// Created by Sofia Iannicelli on 3/3/25.
//

#include "Pathtracer.h"
#include <iostream>
#include <chrono>
// #include <omp.h>

const double EPSILON = 2e-8;

int ***Pathtracer::pathtrace(int numCols, int numRows) {
    auto start_time = std::chrono::steady_clock::now();
    // initialize pixelColors multi-dimensional array
    int*** pixelColors = new int**[numRows];
    for(int i = 0; i < numRows; ++i) {
        pixelColors[i] = new int*[numCols];
        for(int j = 0; j < numCols; ++j) {
            pixelColors[i][j] = new int[3];
        }
    }

    WorldSpaceCoord worldCoords  = calculateWorldSpaceCoords(numCols, numRows);

    // random helper
    double invRPPPS = 1. / raysPerPixelPerSide;
    RandomNumber randNum;
    std::uniform_real_distribution<double> dist = randNum.get_dist(-0.5, 0.5);

    int max = 0;

    // handle each pixel
    // omp_set_num_threads(1);
#pragma omp parallel for
    for(int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            if(j == 443 && i == 73) {
                // std::cerr << "stop here" << std::endl;
            }
            Eigen::Vector3d pixelColor(0, 0, 0);
            for (int m = 0; m < raysPerPixelPerSide; ++m)
            {
                for (int n = 0; n < raysPerPixelPerSide; ++n)
                {
                    double jitterU = randNum.get_random_double_from_dist(dist);
                    double jitterV = randNum.get_random_double_from_dist(dist);
                    double uCoord = ((j + 0.5 * invRPPPS) + ((double(m) + jitterU) * invRPPPS)) * worldCoords.uInc;
                    double vCoord = ((i + 0.5 * invRPPPS) + ((double(n) + jitterV) * invRPPPS)) * worldCoords.vInc;

                    pixelColor += getPathLi((uCoord - worldCoords.maxU) * worldCoords.uAxis +
                                               (worldCoords.maxV - vCoord) * worldCoords.vAxis);
                }
            }
            pixelColor /= double(raysPerPixelPerSide * raysPerPixelPerSide);
            pixelColor = pixelColor.cwiseMin(1.).cwiseMax(0.);
            pixelColors[i][j][0] = int(pixelColor[0] * 255);
            pixelColors[i][j][1] = int(pixelColor[1] * 255);
            pixelColors[i][j][2] = int(pixelColor[2] * 255);
            max = std::max(std::max(std::max(pixelColors[i][j][0], max), pixelColors[i][j][1]), pixelColors[i][j][2]);
        }
        std::cout << "row" << i << std::endl;
        // std::cout << "max" << max << std::endl;
    }

    if(LOG_TIME) {
        auto end_time = std::chrono::steady_clock::now();
        auto elapsed = (end_time - start_time);
        std::cerr << "Time spent raytracing: ";
        std::cerr << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << "[ms]";
        std::cerr << " (" << std::chrono::duration_cast<std::chrono::seconds>(elapsed).count() << "[s]" << ")";
        std::cerr << " (" << std::chrono::duration_cast<std::chrono::minutes>(elapsed).count() << "[m]" << ")";
        std::cerr << std::endl;
    }

    return pixelColors;
}

Eigen::Vector3d Pathtracer::getPathLi(const Eigen::Vector3d& target) {
    // start at eye
    Eigen::Vector3d origin = scene.getCamera().getLookFrom();
    Eigen::Vector3d direction = target - origin;
    Ray ray{origin, direction};

    std::vector<double> iors;
    iors.push_back(scene.getAmbientIor());

    int bounceCount = 1;

    // do multiple paths and average
    Eigen::Vector3d liSum = Eigen::Vector3d::Zero();
    for(int i = 0; i < pathsPerPixel; ++i) {
        liSum += getPathLi(ray, bounceCount, iors);
    }
    return (1. / pathsPerPixel) * liSum;
}

Eigen::Vector3d Pathtracer::getPathLi(Ray ray, int bounceCount, std::vector<double> iors) {
    // if the maximum number of rays have been reached, this one will not contribute and stop recursion
    if(bounceCount > maxDepth) {
        return {0, 0, 0};
    }

    // find ray-scene intersection
    RayHit hit = getClosestIntersection(ray);

    // If no geo was intersected, return background
    if(!hit.isHit) {
        return scene.getBackgroundIntensity(ray);
    }

    // something was hit
    Eigen::Vector3d toView = (ray.getOrigin() - hit.point).normalized();

    // get direct contributions at hit
    Eigen::Vector3d totalIntensity = Eigen::Vector3d::Zero();
    for(auto light : scene.lights) {
        LightHit lightHit = light->getLightHit(hit.point, hit.normal);

        Eigen::Vector3d reflectDirection = (((2 * (hit.normal.dot(toView))) * hit.normal) -
                                            toView).normalized();
        Eigen::Vector3d halfVector = (lightHit.direction + toView).normalized();

        // if(!lightHit.hit) continue;
        if(lightHit.hit) {
            double shadowAmt;
            if(lightHit.distance == 0) {
                // on emissive light
                // FIXME?
                shadowAmt = 0;
            }
            else {
                shadowAmt = getInShadow(hit.point, lightHit.direction, light, lightHit.distance);
            }
            totalIntensity += (1. - shadowAmt) *
                              (hit.brdf->evalIndirect(toView, hit.normal, halfVector, lightHit.direction,
                                                      hit.point, hit.u, hit.v).array() * lightHit.intensity.array()).matrix();
        }

        // send out new ray according to brdf
        BRDFInfo brdfInfo{hit.normal, reflectDirection, toView, halfVector, lightHit.direction, hit.backFace};
        double prob;
        bool valid;
        // TODO: is the iors messed up for the next light?
        Eigen::Vector3d nextDirection = hit.brdf->sample(brdfInfo, iors, prob, valid);
        if(valid) {
            Ray nextRay = Ray(hit.point + EPSILON * nextDirection, nextDirection);

            totalIntensity += prob * getPathLi(nextRay, ++bounceCount, iors);
        }
    }

    // FIXME: one new ray per hit or one new ray per hit per light?!?
    /*
    // send out new ray according to brdf
    double prob;
    bool valid;
    BRDFInfo brdfInfo{hit.normal, };
    Eigen::Vector3d nextDirection = hit.brdf->sample(hit.normal, prob, valid);
    if(valid) {
        Ray nextRay = Ray(hit.point + EPSILON * nextDirection, nextDirection);

        totalIntensity += prob * getPathLi(nextRay, ++bounceCount, iors);
    }
     */

    return totalIntensity;
}

double Pathtracer::getInShadow(const Eigen::Vector3d& point, const Eigen::Vector3d& direction,
                   const std::shared_ptr<Light>& light, const double& distToLight) {
    // FIXME: for now, assume that the direction is going to the light
    double inShadow = 0;

    // check for geo in the way of the path to the light
    Eigen::Vector3d shadowRayOrigin = point + EPSILON * direction;
    Ray shadowRay = Ray(shadowRayOrigin, direction);

    // FIXME: don't have to get closest first...
    // FIXME: optimization with excluding already checked?
    double tCovered = 0;
    while(inShadow < 1.) {
        RayHit geoHitInfo = getClosestIntersection(shadowRay);

        // if nothing is hit or if you're past the light, stop looking
        if(!geoHitInfo.isHit || (tCovered + geoHitInfo.t) >= distToLight) break;

        /* handle refractive
        if(geoHitInfo.brdf->getIsRefractive()) {
            inShadowPart += (1 - geoHitInfo.material->getRefractionK());
            // Start new ray after the point already hit
            tCovered += geoHitInfo.t;

            shadowRayOrigin = geoHitInfo.point + EPSILON * shadowRayDirection;
            shadowRay = Ray(shadowRayOrigin, shadowRayDirection);
        }
         */
        inShadow += 1;
        break;
    }
    // inShadow /= double(light->shadowRayCount);

    return inShadow;
}

RayHit Pathtracer::getClosestIntersection(const Ray& ray) {
    return scene.bvh->findRayHit(ray);
}

WorldSpaceCoord Pathtracer::calculateWorldSpaceCoords(int numCols, int numRows) {
    Camera camera = scene.getCamera();

    Eigen::Vector3d viewRay = (camera.getLookAt() - camera.getLookFrom());
    double distToCenter = viewRay.norm();

    // calculate the maximum u and v values based on the FOV
    double maxU = distToCenter * tan(camera.getFovXRad()/2);
    double maxV = distToCenter * tan(camera.getFovYRad()/2);

    // calculate the u and v axes
    Eigen::Vector3d uAxis = viewRay.cross(camera.getCameraLookUp()).normalized();
    Eigen::Vector3d vAxis = uAxis.cross(viewRay).normalized();

    // calculate the u and v increments for a change in pixel
    double uInc = 2 * maxU / numCols;
    double vInc = 2 * maxV / numRows;

    WorldSpaceCoord coord(maxU, maxV, uInc, vInc, uAxis, vAxis);
    return coord;

}
