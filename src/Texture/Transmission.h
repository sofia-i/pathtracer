//
// Created by Sofia Iannicelli on 3/5/25.
//

#ifndef PATHTRACER_TRANSMISSION_H
#define PATHTRACER_TRANSMISSION_H

#include "BRDF.h"

class Transmission : public BRDF {
public:
    Transmission(double ior) : Transmission(ior, 0.0) {}
    Transmission(double ior, double jitterAmt) : BRDF(), ior(ior), jitterAmt(jitterAmt),
                                                 dist(std::uniform_real_distribution<double>(-0.5, 0.5)) {}

    Eigen::Vector3d sample(BRDFInfo info, std::vector<double>& iors, double &prob, bool &valid) override {
        Eigen::Vector3d rayD = -info.view;
        double cosIn = info.normal.dot(rayD);
        Eigen::Vector3d normalRef = info.normal;
        if(cosIn < 0) {
            // outside
            cosIn = -cosIn;
        }
        else {
            // inside
            normalRef = -info.normal;
        }
        // get ior ratio
        double iorIn, iorOut, iorRatio;
        getIorAcrossIntersection(info.isBackFace, iorIn, iorOut, iorRatio, iors);
        // Find parallel and orthogonal portions of refraction direction
        Eigen::Vector3d refractDirP = iorRatio * (rayD + cosIn * normalRef);
        // NOTE: added max
        Eigen::Vector3d refractDirS = - std::sqrt(std::max(0., 1. - std::pow(refractDirP.norm(), 2))) * normalRef;

        Eigen::Vector3d refractDirection = refractDirP + refractDirS;
        // jitter refraction
        refractDirection += jitterAmt * Eigen::Vector3d{dist(mt), dist(mt), dist(mt)};
        refractDirection.normalize();

        prob = 1.;

        assert(!std::isnan(refractDirection.x()) &&
                !std::isnan(refractDirection.y()) &&
                !std::isnan(refractDirection.z()) &&
                "Ray direction NaN");
        valid = true;

        return refractDirection;
    }

    Eigen::Vector3d sample(const Eigen::Vector3d& view, const Eigen::Vector3d& normal, const bool& isBackFace,
                           const double& iorIn, const double& iorOut, const double& iorRatio,
                           double& prob, bool& valid) {
        Eigen::Vector3d rayD = -view;
        double cosIn = normal.dot(rayD);
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
        // NOTE: added max
        Eigen::Vector3d refractDirS = - std::sqrt(std::max(0., 1. - std::pow(refractDirP.norm(), 2))) * normalRef;

        Eigen::Vector3d refractDirection = refractDirP + refractDirS;
        // jitter refraction
        refractDirection += jitterAmt * Eigen::Vector3d{dist(mt), dist(mt), dist(mt)};
        refractDirection.normalize();

        prob = 1.;

        assert(!std::isnan(refractDirection.x()) &&
               !std::isnan(refractDirection.y()) &&
               !std::isnan(refractDirection.z()) &&
               "Ray direction NaN");
        valid = true;

        return refractDirection;
    }

    Eigen::Vector3d eval(const Eigen::Vector3d& incoming, const Eigen::Vector3d& outgoing,
                         const Eigen::Vector3d& normal, const Eigen::Vector3d& point,
                         double u, double v) override {
        // TODO
    }

    Eigen::Vector3d evalIndirect(const Eigen::Vector3d& view, const Eigen::Vector3d& normal,
                                 const Eigen::Vector3d& half, const Eigen::Vector3d& toLight,
                                 const Eigen::Vector3d& point, double u, double v) override {
        return {0., 0., 0.};
    }

    std::string toString() override {
        std::stringstream ss;
        ss << "Transmission" << std::endl;
        return ss.str();
    }

public:
    /**
     *
     * @param isBackFace
     * @param[out] iorIn
     * @param[out] iorOut
     * @param[out] iorRatio
     * @param iors
     */
    void getIorAcrossIntersection(bool isBackFace, double& iorIn, double& iorOut, double& iorRatio,
                                  std::vector<double>& iors) const {
        // FIXME: if reflexive overlapping it won't work
        if(isBackFace && iors.size() > 1) {
            // coming out of material
            iorIn = iors.back();
            iors.pop_back();
            iorOut = iors.back();
            iorRatio = iorIn / iorOut;
        }
        else {
            iorIn = iors.back();
            iorOut = ior;
            iorRatio = iorIn / iorOut;
            iors.push_back(iorOut);
        }
    }

private:
    std::uniform_real_distribution<double> dist;
    double ior;

    double jitterAmt;
};

#endif //PATHTRACER_TRANSMISSION_H
