//
// Created by Sofia Iannicelli on 1/21/25.
//

#ifndef RAYTRACER_2_MATERIAL_H
#define RAYTRACER_2_MATERIAL_H

#include "Texture.h"
#include <sstream>

class Material {
public:
    Material(double diffuseK, double specularK, double ambientK, double glsK,
               std::shared_ptr<Texture> diffuse, Eigen::Vector3d specularColor,
               double refl, double reflJitter, double transJitter);

    Material(double diffuseK, double specularK, double ambientK, double glsK,
               std::shared_ptr<Texture> diffuse, Eigen::Vector3d specularColor,
               double refl, double reflJitter, double transJitter,
               double ior, double refractionK);

    double getDiffuseK() const { return diffuseK; }
    double getSpecularK() const { return specularK; }
    double getAmbientK() const { return ambientK; }
    virtual Eigen::Vector3d getDiffuseColor(double u, double v, const Eigen::Vector3d& p) const;
    Eigen::Vector3d getSpecularColor() const { return specularColor; }
    double getGlsK() const { return glsK; }
    double getRefl() const { return refl; }
    bool getIsRefractive() const { return refractive; }
    double getIOR() const { return ior; }
    double getRefractionK() const { return refractionK; }
    double getReflJitter() const { return reflJitter; }
    double getTransJitter() const { return transJitter; }

    friend std::ostream& operator<<(std::ostream& os, const Material& mat) {
        os << mat.toString();
        return os;
    }

    virtual std::string toString() const {
        std::string s;
        std::stringstream ss(s);

        ss << "\tKd: " << diffuseK << std::endl;
        ss << "\tKs: " << specularK << std::endl;
        ss << "\tKa: " << ambientK << std::endl;
        ss << "\tDiffuse Color: " << *diffuse << std::endl;
        ss << "\tSpecular Color: " << specularColor << std::endl;
        ss << "\tKgls: " << glsK << std::endl;
        ss << "\tReflectivity: " << refl << std::endl;
        if(refractive) {
            ss << "\tTransmission: " << refractionK << std::endl;
            ss << "\tIOR: " << ior << std::endl;
        }

        return ss.str();
    }

protected:
    double diffuseK;  // diffuse coefficient
    double specularK;  // specular coefficient
    double ambientK;  // ambient coefficient
    double glsK;  // gloss coefficient
    // vec3<double> diffuseColor;
    std::shared_ptr<Texture> diffuse;
    Eigen::Vector3d specularColor;
    double refl;

    bool refractive;
    double ior;
    double refractionK;

    double reflJitter;
    double transJitter;

};


#endif //RAYTRACER_2_MATERIAL_H
