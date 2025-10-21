//
// Created by Sofia Iannicelli on 3/4/25.
// based on https://boksajak.github.io/files/CrashCourseBRDF.pdf
//

#ifndef PATHTRACER_BRDF_H
#define PATHTRACER_BRDF_H

#include <Eigen/Eigen>
#include <random>
#include <cmath>
#include "Texture.h"
#include "../MathUtils/Sampler.h"

struct BRDFInfo {
    const Eigen::Vector3d normal;
    const Eigen::Vector3d refl;
    const Eigen::Vector3d view;
    const Eigen::Vector3d half;
    const Eigen::Vector3d toLight;
    const bool isBackFace;
};

class BRDF {
public:
    BRDF() {
        mt = std::mt19937(rd());
        // dist = std::uniform_real_distribution<double>(-1.0, 1.0);
    }

    /**
     * Sample the next ray direction
     * @param[in] info Ray hit information
     * @param[in,out] iors Stack of ior values along path
     * @param[out] prob probability of output direction
     * @param[out] valid Valid sample was found
     * @return outgoing direction
     */
    virtual Eigen::Vector3d sample(BRDFInfo info, std::vector<double>& iors, double& prob, bool& valid) {
        valid = false;
    }

    /**
     *
     * @param incoming
     * @param outgoing
     * @return
     */
    virtual Eigen::Vector3d eval(const Eigen::Vector3d& incoming, const Eigen::Vector3d& outgoing,
                                 const Eigen::Vector3d& normal, const Eigen::Vector3d& point,
                                 double u, double v) = 0;

    /**
     * Get the intensity (color) based on viewing information (point, normal, view ray, normal, etc.)
     * @param view
     * @param normal
     * @param half
     * @param toLight
     * @param point
     * @param u
     * @param v
     * @return Intensity vector
     */
    virtual Eigen::Vector3d evalIndirect(const Eigen::Vector3d& view, const Eigen::Vector3d& normal,
                                         const Eigen::Vector3d& half, const Eigen::Vector3d& toLight,
                                         const Eigen::Vector3d& point, double u, double v) = 0;

    /**
     * evaluates probability density function based on outgoing
     * @param outgoing
     * @return probability density
     */
    // virtual double pdf(Eigen::Vector3d outgoing) = 0;

    virtual std::string toString() = 0;

protected:
    std::random_device rd;
    std::mt19937 mt;


};


#endif //PATHTRACER_BRDF_H
