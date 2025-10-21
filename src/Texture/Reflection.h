//
// Created by Sofia Iannicelli on 3/5/25.
//

#ifndef PATHTRACER_REFLECTION_H
#define PATHTRACER_REFLECTION_H

#include "BRDF.h"

class Reflection : public BRDF {
public:
    Reflection() : Reflection(0.0) {}
    explicit Reflection(double jitterAmt) : BRDF(), jitterAmt(jitterAmt),
                                            dist(std::uniform_real_distribution<double>(-0.5, 0.5)) {}

    Eigen::Vector3d sample(BRDFInfo info, std::vector<double>& iors, double &prob, bool &valid) override {
        return sample(info.normal, prob, valid);
    }

    Eigen::Vector3d sample(const Eigen::Vector3d &local, double &prob, bool &valid) {
        prob = 1.;
        valid = true;
        Eigen::Vector3d sample = local;
        sample += jitterAmt * Eigen::Vector3d{dist(mt), dist(mt), dist(mt)};
        sample.normalize();
        return sample;
    }

    Eigen::Vector3d eval(const Eigen::Vector3d& incoming, const Eigen::Vector3d& outgoing,
                         const Eigen::Vector3d& normal, const Eigen::Vector3d& point,
                         double u, double v) override {
        // TODO
    }

    Eigen::Vector3d evalIndirect(const Eigen::Vector3d& view, const Eigen::Vector3d& normal,
                                 const Eigen::Vector3d& half, const Eigen::Vector3d& toLight,
                                 const Eigen::Vector3d& point, double u, double v) override {
        return {0.0, 0.0, 0.0};
    }

    std::string toString() override {
        std::stringstream ss;
        ss << "Reflection" << std::endl;
        return ss.str();
    }

private:
    std::uniform_real_distribution<double> dist;
    double jitterAmt;

};

#endif //PATHTRACER_REFLECTION_H
