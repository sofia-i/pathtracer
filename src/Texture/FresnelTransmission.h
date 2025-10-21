//
// Created by Sofia Iannicelli on 3/5/25.
//

#ifndef PATHTRACER_FRESNELTRANSMISSION_H
#define PATHTRACER_FRESNELTRANSMISSION_H

#include "Transmission.h"
#include "Reflection.h"

class FresnelTransmission : public BRDF {
public:
    explicit FresnelTransmission(double reflAmt, double ior, double reflJitter, double transJitter) : BRDF(),
            reflection(reflJitter), transmission(ior, transJitter), materialReflectanceAmt(reflAmt),
            dist(std::uniform_real_distribution<double>(0., 1.)) {}

    Eigen::Vector3d sample(BRDFInfo info, std::vector<double> &iors, double &prob, bool &valid) override
    {
        double iorIn, iorOut, iorRatio;
        transmission.getIorAcrossIntersection(info.isBackFace, iorIn, iorOut, iorRatio, iors);

        double reflAmt = getPortionReflected(info.normal, -info.view, iorIn, iorOut);

        double u = dist(mt);
        if(u < reflAmt) {
            // reflection
            Eigen::Vector3d sample = reflection.sample(info, iors, prob, valid);
            prob *= reflAmt;  // FIXME??
            return sample;
        }
        else {
            // transmission
            // return transmission.sample(info, iors, prob, valid);
            Eigen::Vector3d sample = transmission.sample(info.view, info.normal, info.isBackFace,
                                                         iorIn, iorOut, iorRatio, prob, valid);
            prob *= (1. - reflAmt);  // FIXME??
            return sample;
        }
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
        ss << "FresnelTransmission" << std::endl;
        ss << "\t" << "Reflection " << reflection.toString() << std::endl;
        ss << "\t" << "Transmission " << transmission.toString() << std::endl;
        ss << "\t" << "Refl " << materialReflectanceAmt << std::endl;
        return ss.str();
    }

private:
    double getPortionReflected(const Eigen::Vector3d& normal, const Eigen::Vector3d& rayD,
                               const double& iorIn, const double& iorOut) {
        // Schlick’s approximation https://blog.demofox.org/2017/01/09/raytracing-reflection-refraction-fresnel-total-internal-reflection-and-beers-law/
        double r0 = (iorIn - iorOut) / (iorIn + iorOut);
        r0 *= r0;
        double cosIn = -normal.dot(rayD);
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
        /*
        if(reflFresnel > 1.0) {
            std::cerr << "gt one" << std::endl;
        }
         */
        return materialReflectanceAmt + (1. - materialReflectanceAmt) * reflFresnel;
    }

private:
    Reflection reflection;
    Transmission transmission;

    std::uniform_real_distribution<double> dist;
    double materialReflectanceAmt;
};

#endif //PATHTRACER_FRESNELTRANSMISSION_H
