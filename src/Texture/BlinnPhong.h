//
// Created by Sofia Iannicelli on 3/5/25.
//

#ifndef PATHTRACER_BLINNPHONG_H
#define PATHTRACER_BLINNPHONG_H

#include "BRDF.h"
#include "Lambertian.h"
#include "Phong.h"

class BlinnPhong : public BRDF {
public:
    BlinnPhong(std::shared_ptr<Texture> diffuse, std::shared_ptr<Texture> specular,
               std::shared_ptr<Texture> ambient, double shininess, double diffuseK,
               double specularK, double ambientK) :
            BRDF(), lambertian(diffuse), phong(specular, shininess), ambientTexture(ambient),
            diffuseK(diffuseK), specularK(specularK), ambientK(ambientK)
    {
        dist = std::uniform_real_distribution<double>{0., 1.};
    }

    Eigen::Vector3d sample(BRDFInfo info, std::vector<double> &iors, double &prob, bool &valid) override {
        // TODO: validate
        double u = dist(mt);
        if(u < diffuseK) {
            Eigen::Vector3d sample = lambertian.sample(info.normal, prob, valid);
            prob *= diffuseK;
            return sample;
        }
        else if (u < (diffuseK + specularK)) {
            Eigen::Vector3d sample = phong.sample(info.refl, prob, valid);
            prob *= specularK;
            return sample;
        }
        else {
            valid = false;
            return {};
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
        return (diffuseK * lambertian.evalIndirect(view, normal, half, toLight, point, u, v) +
                specularK * phong.evalIndirect(view, normal, half, toLight, point, u, v) +
                ambientK * ambientTexture->getIntensity(u, v, point));
    }

    std::string toString() override {
        std::stringstream ss;
        ss << "BlinnPhong" << std::endl;
        ss << "\t" << "Lambertian" << lambertian.toString() << std::endl;
        ss << "\t" << "Phong" << phong.toString() << std::endl;
        ss << "\t" << "diffuseK: " << diffuseK << " specularK: " << specularK << " ambientK: " << ambientK << std::endl;
        return ss.str();
    }

private:
    std::uniform_real_distribution<double> dist;
    Lambertian lambertian;
    Phong phong;
    std::shared_ptr<Texture> ambientTexture;

    double diffuseK;
    double specularK;
    double ambientK;
};


#endif //PATHTRACER_BLINNPHONG_H
