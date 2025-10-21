//
// Created by Sofia Iannicelli on 1/21/25.
//

#include "Material.h"

Material::Material(double diffuseK, double specularK, double ambientK, double glsK,
                   std::shared_ptr<Texture> diffuse, Eigen::Vector3d specularColor,
                   double refl, double reflJitter, double transJitter) :
            diffuseK(diffuseK), specularK(specularK), ambientK(ambientK), glsK(glsK),
            diffuse(diffuse), specularColor(specularColor), refl(refl),
            reflJitter(reflJitter), transJitter(transJitter),
            refractive(false) {}

Material::Material(double diffuseK, double specularK, double ambientK, double glsK,
                    std::shared_ptr<Texture> diffuse, Eigen::Vector3d specularColor,
                    double refl, double reflJitter, double transJitter,
                    double ior, double refractionK) :
            diffuseK(diffuseK), specularK(specularK), ambientK(ambientK), glsK(glsK),
            diffuse(diffuse), specularColor(specularColor), refl(refl),
            reflJitter(reflJitter), transJitter(transJitter),
            ior(ior), refractionK(refractionK) {
    refractive = refractionK > 0.0;
}

Eigen::Vector3d Material::getDiffuseColor(double u, double v, const Eigen::Vector3d& p) const {
    return diffuse->getIntensity(u, v, p);
}
