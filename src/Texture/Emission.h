//
// Created by Sofia Iannicelli on 3/6/25.
//

#ifndef PATHTRACER_EMISSION_H
#define PATHTRACER_EMISSION_H

#include "BRDF.h"
#include "Texture.h"

class Emission : public BRDF {

public:
    explicit Emission(const std::shared_ptr<Texture>& texture) : BRDF(), texture(texture) {}

    Eigen::Vector3d sample(BRDFInfo info, std::vector<double> &iors, double &prob, bool &valid) override {
        valid = false;
    }

    Eigen::Vector3d eval(const Eigen::Vector3d& incoming, const Eigen::Vector3d& outgoing,
                         const Eigen::Vector3d& normal, const Eigen::Vector3d& point,
                         double u, double v) override {
        // TODO
        return {};
    }

    Eigen::Vector3d
    evalIndirect(const Eigen::Vector3d& view, const Eigen::Vector3d& normal,
                 const Eigen::Vector3d& half, const Eigen::Vector3d& toLight,
                 const Eigen::Vector3d& point, double u, double v) override {
        return texture->getIntensity(u, v, point);
    }

    std::shared_ptr<Texture> getTexture() { return texture; }

    std::string toString() override {
        std::stringstream ss;
        ss << "Emission" << std::endl;
        ss << "\t" << "Texture " << texture << std::endl;
        return ss.str();
    }

private:
    std::shared_ptr<Texture> texture;

};

#endif //PATHTRACER_EMISSION_H
