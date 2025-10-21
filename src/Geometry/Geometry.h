//
// Created by Sofia Iannicelli on 1/31/25.
//

#ifndef RAYTRACER_2_GEOMETRY_H
#define RAYTRACER_2_GEOMETRY_H

#include <memory>
#include <cstdio>
#include <sstream>
#include <string>
#include <utility>
#include "../Ray.hpp"
#include "Hittable.h"
#include "../Texture/Material.h"
#include "../Texture/BRDF.h"
#include "Extent.h"

class Geometry : public Hittable {
public:
    Geometry(const std::shared_ptr<BRDF>& brdf, std::string description) :
            brdf(brdf), description(std::move(description)) { }

    virtual ~Geometry() = default; // I. destructor
    Geometry(const Geometry& other) = default; // II. copy constructor
    Geometry& operator=(const Geometry& other) = default; // III. copy assignment
    Geometry(Geometry&& other) noexcept = default;// IV. move constructor
    Geometry& operator=(Geometry&& other) noexcept = default; // V. move assignment

    virtual Extent findExtent() = 0;

    std::shared_ptr<BRDF> getBRDF() const { return brdf; }
    std::string getDescription() const { return description; }

    friend std::ostream& operator<<(std::ostream& os, const Geometry& geo) {
        os << geo.toString();
        return os;
    }

    virtual std::string toString() const {
        std::string s;
        std::stringstream ss(s);

        // ss << getDescription() << std::endl;
        ss << "\t" << brdf->toString() << std::endl;

        return ss.str();
    }

protected:
    std::shared_ptr<BRDF> brdf;
    std::string description;

};


#endif //RAYTRACER_2_GEOMETRY_H
