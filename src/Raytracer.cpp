//
//  Raytracer.cpp
//  raytracer_2
//
//  Created by Sofia Iannicelli on 2/15/23.
//

#include "Raytracer.hpp"
#include "Ray.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <chrono>
#include "MathUtils/Utils.h"

void Raytracer::showProgress(int index, int total) {
    int barWidth = 70;

    std::cout << "[";
    float progress = float(index) / float(total);
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "%" << std::endl;
}

int*** Raytracer::raytrace(int numCols, int numRows) {
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

    // handle each pixel
    for(int i = 0; i < numRows; ++i) {
        for(int j = 0; j < numCols; ++j) {
            Eigen::Vector3d pixelColor(0, 0, 0);
            for(int m = 0; m < raysPerPixelPerSide; ++m) {
                for(int n = 0; n < raysPerPixelPerSide; ++n) {
                    double jitterU = randNum.get_random_double_from_dist(dist);
                    double jitterV = randNum.get_random_double_from_dist(dist);
                    double uCoord = ((j + 0.5 * invRPPPS) + ((double(m) + jitterU) * invRPPPS)) * worldCoords.uInc;
                    double vCoord = ((i + 0.5 * invRPPPS) + ((double(n) + jitterV) * invRPPPS)) * worldCoords.vInc;
                    pixelColor += getRayResult((uCoord - worldCoords.maxU) * worldCoords.uAxis +
                                               (worldCoords.maxV - vCoord) * worldCoords.vAxis);
                }
            }
            pixelColor = pixelColor / double(raysPerPixelPerSide * raysPerPixelPerSide);
            pixelColors[i][j][0] = pixelColor[0];
            pixelColors[i][j][1] = pixelColor[1];
            pixelColors[i][j][2] = pixelColor[2];
            if(LOG_TIME && (i * numRows + j) % LOG_INTERVAL == 0) showProgress(i, numRows - 1);
        }
    }

    auto end_time = std::chrono::steady_clock::now();

    if(LOG_TIME) {
        auto elapsed = (end_time - start_time);
        std::cerr << "Time spent raytracing: ";
        std::cerr << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << "[ms]";
        std::cerr << " (" << std::chrono::duration_cast<std::chrono::seconds>(elapsed).count() << "[s]" << ")";
        std::cerr << " (" << std::chrono::duration_cast<std::chrono::minutes>(elapsed).count() << "[m]" << ")";
        std::cerr << std::endl;
    }
    
    return pixelColors;
}

double Raytracer::getInShadow(const Eigen::Vector3d& intersectPt, const std::shared_ptr<Light>& light) {
    Eigen::Vector3d shadowRayDirection;
    Eigen::Vector3d shadowRayOrigin;
    double inShadow = 0;

    // Iterate over each ray to light (1 for point & direction, multiple for area)
    for(int i = 0; i < light->shadowRayCount; ++i) {
        // calculate shadow ray
        double distToLight;
        bool hitLight;
        // find direction to and distance to light
        light->getPathToLight(intersectPt, i, hitLight, shadowRayDirection, distToLight);

        if(!hitLight) {
            // TODO: ?
            inShadow += 1;
            continue;
        }

        // check for geo in the way of the path to the light
        shadowRayOrigin = intersectPt + EPSILON * shadowRayDirection;
        Ray shadowRay = Ray(shadowRayOrigin, shadowRayDirection);
        double inShadowPart = 0;
        // FIXME: don't have to get closest first...
        // FIXME: optimization with excluding already checked?
        double tCovered = 0;
        while(inShadowPart < 1.) {
            RayHit geoHitInfo = getClosestIntersection(shadowRay);

            // if nothing is hit or if you're past the light, stop looking
            if(!geoHitInfo.isHit || (tCovered + geoHitInfo.t) >= distToLight) break;

            if(geoHitInfo.material->getIsRefractive()) {
                inShadowPart += (1 - geoHitInfo.material->getRefractionK());
                // Start new ray after the point already hit
                tCovered += geoHitInfo.t;

                shadowRayOrigin = geoHitInfo.point + EPSILON * shadowRayDirection;
                shadowRay = Ray(shadowRayOrigin, shadowRayDirection);
            }
            else {
                inShadowPart += 1;
                break;
            }
        }

        inShadow += std::min(inShadowPart, 1.);
    }
    inShadow /= double(light->shadowRayCount);

    return inShadow;
}

inline Eigen::Vector3d Raytracer::getAmbient(const std::shared_ptr<Material>& mat, const Eigen::Vector3d& point,
                                          const double u, const double v) {
    return mat->getAmbientK() * scene.getAmbientLight() * mat->getDiffuseColor(u, v, point);
}

inline Eigen::Vector3d Raytracer::getDiffuse(const std::shared_ptr<Material>& mat,
                                          const Eigen::Vector3d& point,
                                          const double u, const double v,
                                          const std::shared_ptr<Light>& light,
                                          const Eigen::Vector3d& normal, const Eigen::Vector3d& toLight) {
    return mat->getDiffuseK() * light->getLightColor() * mat->getDiffuseColor(u, v, point) * std::max(0.0, dot(normal, toLight));
}

inline Eigen::Vector3d Raytracer::getSpecular(const std::shared_ptr<Material>& mat,
                                           const std::shared_ptr<Light>& light,
                                           const Eigen::Vector3d& normal, const Eigen::Vector3d& toLight,
                                           const Eigen::Vector3d& view) {
    Eigen::Vector3d reflection = ((2 * (dot(normal, toLight))) * normal) - toLight;
    double specularK = mat->getSpecularK();
    double glsK = mat->getGlsK();
    const Eigen::Vector3d& specularColor = mat->getSpecularColor();
    return specularK * light->getLightColor() * specularColor *
        std::pow(std::max(0.0, dot(view, reflection)), glsK);
}

/**
 * Compute illumination equation
 * includes intensity from ambient, specular, diffuse
 * not from reflection, transmission
 * @param mat
 * @param normal
 * @param view
 * @param intersectPt
 * @return
 */
Eigen::Vector3i Raytracer::illuminationEq(const std::shared_ptr<Material>& mat, const double u, const double v,
                                    const Eigen::Vector3d& normal,
                                    const Eigen::Vector3d& view, const Eigen::Vector3d& intersectPt) {
    // all the incoming vectors should be normalized
    assert(normal.isNormalized() && view.isNormalized() && "normal and view rays should be normalized");
    Eigen::Vector3i colorResult;
    Eigen::Vector3d colorSum(0.0, 0.0, 0.0);

    // compute ambient contribution
    colorSum += getAmbient(mat, intersectPt, u, v);

    for(auto&& light : scene.lights) {
        // compute shadow information
        double shadowAmt = getInShadow(intersectPt, light);

        // calculate intensity if not in complete shadow
        if(shadowAmt < 1) {
            Eigen::Vector3d toLight = light->getDirectionToLight(intersectPt);
            // compute diffuse contribution
            colorSum += (1 - shadowAmt) * getDiffuse(mat, intersectPt, u, v, light, normal, toLight);
            // compute specular contribution
            colorSum += (1 - shadowAmt) * getSpecular(mat, light, normal, toLight, view);
        }
    }

    colorResult = toIntVec3(255 * colorSum);
    return colorResult;
}

Eigen::Vector3i Raytracer::getRayResult(Eigen::Vector3d target) {
    // computes result of raytracing through a specified world-space coordinate
    Eigen::Vector3d rayOrigin = scene.getCamera().getLookFrom();
    Eigen::Vector3d rayDirection = target - rayOrigin;
    Ray ray = Ray(rayOrigin, rayDirection);
    std::stack<double> iors;
    iors.push(scene.getAmbientIor());
    return getRayResult(ray, 1, iors);
}

RayHit Raytracer::getClosestIntersection(const Ray& ray) {
    if(USE_BOUNDING_VOLUME) {
        return scene.bvh->findRayHit(ray);
    }
    else {
        Eigen::Vector3d normal;

        int closestGeoIdx = -1;
        RayHit closestHitInfo = RayHit::Miss();

        // iterate over all geo to test each
        for(int i = 0; i < scene.geo.size(); ++i) {
            RayHit hitInfo = scene.geo[i]->findRayHit(ray);
            // if the ray intersects the geo, check to see if the geo is the closest one hit (so far)
            if(hitInfo.isHit && (closestGeoIdx == -1 || hitInfo.t < closestHitInfo.t)) {
                // update the closest intersected geo
                closestGeoIdx = i;
                closestHitInfo = hitInfo;
            }
        }

        return closestHitInfo;
    }
}

void Raytracer::getIorAcrossIntersection(const std::shared_ptr<Material>& mat, bool isBackFace,
                                         double& iorIn, double& iorOut, double& iorRatio,
                                         std::stack<double>& iors) {
    // FIXME: if reflexive overlapping it won't work
    if(isBackFace && iors.size() > 1) {
        // coming out of material
        iorIn = iors.top();
        iors.pop();
        iorOut = iors.top();
        iorRatio = iorIn / iorOut;
    }
    else {
        iorIn = iors.top();
        iorOut = mat->getIOR();
        iorRatio = iorIn / iorOut;
        iors.push(iorOut);
    }
}

double Raytracer::getPortionReflected(const Eigen::Vector3d& normal, const Eigen::Vector3d& rayD, const double matRefl,
                                      const double& iorIn, const double& iorOut, const double& iorRatio) {
    // Schlick’s approximation https://blog.demofox.org/2017/01/09/raytracing-reflection-refraction-fresnel-total-internal-reflection-and-beers-law/
    double r0 = (iorIn - iorOut) / (iorIn + iorOut);
    r0 *= r0;
    double cosIn = -dot(normal, rayD);
    // account for backface
    if(cosIn < 0) {
        cosIn = -cosIn;
    }
    if(iorIn > iorOut) {
        double n = iorIn / iorOut;
        double sinOutSq = n * n * (1.0 - cosIn * cosIn);
        if(sinOutSq > 1.0) {
            return 1;
        }
        cosIn = sqrt(1.0 - sinOutSq);  // FIXME?
    }
    double cosTerm = (1 - cosIn);
    double reflFresnel = r0 + (1. - r0)*(cosTerm * cosTerm * cosTerm * cosTerm * cosTerm);
    if(reflFresnel > 1.0) {
        std::cerr << "gt one" << std::endl;
    }
    return matRefl + (1. - matRefl) * reflFresnel;
}

inline Ray Raytracer::getTransmissionRay(const std::shared_ptr<Material>& mat, const Eigen::Vector3d& normal,
                                         const Eigen::Vector3d& rayD, const Eigen::Vector3d& intersectPt,
                                         const double iorRatio) const {
    double cosIn = dot(normal, rayD);
    Eigen::Vector3d normalRef = normal;
    if(cosIn < 0) {
        // outside
        cosIn = -cosIn;
    }
    else {
        // inside
        normalRef = -normal;
    }
    // Find parallel and orthogonal portions of refraction direction
    Eigen::Vector3d refractDirP = iorRatio * (rayD + cosIn * normalRef);
    Eigen::Vector3d refractDirS = - std::sqrt(1. - std::pow(refractDirP.norm(), 2)) * normalRef;

    Eigen::Vector3d refractDirection = refractDirP + refractDirS;
    // jitter refraction
    refractDirection += mat->getTransJitter() * Eigen::Vector3d::getRandom(-0.5, 0.5);
    refractDirection = refractDirection.normalized();
    Eigen::Vector3d refractOrigin = intersectPt + (EPSILON * refractDirection);
    return Ray(refractOrigin, refractDirection);
}

inline Eigen::Vector3i Raytracer::getTransmission(const std::shared_ptr<Material>& mat, const Eigen::Vector3d& normal,
                                            const Eigen::Vector3d& rayD, const Eigen::Vector3d& intersectPt,
                                            double iorRatio, int rayCount, std::stack<double>& iors) {
    Ray transmissionRay = getTransmissionRay(mat, normal, rayD, intersectPt, iorRatio);
    return mat->getRefractionK() * getRayResult(transmissionRay, ++rayCount, iors);
}

inline Eigen::Vector3i Raytracer::getReflection(const std::shared_ptr<Material>& mat, const Eigen::Vector3d& normal,
                                          const Eigen::Vector3d& toView, const Eigen::Vector3d& intersectPt,
                                          int rayCount, std::stack<double>& iors) {
    if(mat->getRefl() == 0) {
        return {0, 0, 0};
    }
    Eigen::Vector3d reflectRayDirection = (((2 * (normal.dot(toView))) * normal) -
                                                     toView).normalized();
    // jitter reflection direction
    reflectRayDirection += mat->getReflJitter() * Eigen::Vector3d::getRandom(-0.5, 0.5);
    reflectRayDirection = reflectRayDirection.normalized();
    Eigen::Vector3d reflectRayOrigin = intersectPt + (EPSILON * reflectRayDirection);
    Ray reflectionRay = Ray(reflectRayOrigin, reflectRayDirection);

    return mat->getRefl() * getRayResult(reflectionRay, ++rayCount, iors);
}

Eigen::Vector3i Raytracer::getRayResult(Ray ray, int rayCount, std::stack<double>& iors) {
    // if the maximum number of rays have been reached, this one will not contribute and stop recursion
    if(rayCount > MAX_NUM_RAYS) {
        return {0, 0, 0};
    }
    
    // Find the closest geo intersected by the ray
    // Intersection hit = getClosestIntersection(ray);
    RayHit hit = getClosestIntersection(ray);

    // If no geo was intersected, return the background color
    if(!hit.isHit) {
        return toIntVec3(255 * scene.getBackgroundColor());
    }

    Eigen::Vector3i colorResult(0, 0, 0);

    Eigen::Vector3d toView = getUnitVector(ray.getOrigin() - hit.point);
    Eigen::Vector3i primaryResult = illuminationEq(hit.material, hit.u, hit.v, hit.normal, toView, hit.point);
    colorResult += primaryResult;

    Eigen::Vector3i refractionResult(0, 0, 0);
    // Refraction & reflection
    if(hit.material->getIsRefractive()) {
        // find iorRatio (eta)
        double iorRatio;
        double iorIn;
        double iorOut;
        getIorAcrossIntersection(hit.material, hit.backFace, iorIn, iorOut, iorRatio, iors);

        double kTran = 0;
        double kRefl = getPortionReflected(hit.normal, ray.getDirection(),
                                           hit.material->getRefl(), iorIn, iorOut, iorRatio);
        if(kRefl < 1) {
            // compute refraction
            kTran = 1.0 - kRefl;
            refractionResult = getTransmission(hit.material, hit.normal, ray.getDirection(), hit.point, iorRatio,
                                               rayCount, iors);
        }
        // compute results from reflection
        Eigen::Vector3i reflectionResult = getReflection(hit.material, hit.normal, toView, hit.point, rayCount, iors);

        // combine reflection and refraction based on fresnel equations
        colorResult += kRefl * reflectionResult + kTran * refractionResult;
    }
    // Reflection (no refraction)
    else {
        Eigen::Vector3i reflectionResult = getReflection(hit.material, hit.normal, toView, hit.point, rayCount, iors);
        colorResult += reflectionResult;
    }

    // make sure not to have overflow
    colorResult = clip(colorResult, 0, 255);
    
    return colorResult;
}

WorldSpaceCoord Raytracer::calculateWorldSpaceCoords(int numCols, int numRows) {
    Camera camera = scene.getCamera();
    
    Eigen::Vector3d viewRay = (camera.getLookAt() - camera.getLookFrom());
    double distToCenter = viewRay.norm();

    // calculate the maximum u and v values based on the FOV
    double maxU = distToCenter * tan(camera.getFovXRad()/2);
    double maxV = distToCenter * tan(camera.getFovYRad()/2);
    
    // calculate the u and v axes
    Eigen::Vector3d uAxis = getUnitVector(cross(viewRay, camera.getCameraLookUp()));
    Eigen::Vector3d vAxis = getUnitVector(cross(uAxis, viewRay));
    
    // calculate the u and v increments for a change in pixel
    double uInc = 2 * maxU / numCols;
    double vInc = 2 * maxV / numRows;

    WorldSpaceCoord coord(maxU, maxV, uInc, vInc, uAxis, vAxis);
    return coord;
    
}
