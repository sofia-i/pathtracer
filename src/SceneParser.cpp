//
//  SceneParser.cpp
//  raytracer_2
//
//  Created by Sofia Iannicelli on 2/15/23.
//

#include "SceneParser.hpp"

#include <fstream>
#include <vector>
#include <cassert>
#include <array>
#include "Geometry/Sphere.h"
#include "Geometry/Triangle.hpp"
#include "Geometry/Cylinder.h"
#include "Light/Skylight.h"
#include "Light/DirectionalLight.h"
#include "Light/PointLight.h"
#include "Light/AreaLight.h"
#include "Light/EmissionLight.h"
#include "Texture/Lambertian.h"
#include "Texture/Phong.h"
#include "Texture/Reflection.h"
#include "Texture/Transmission.h"
#include "Texture/FresnelTransmission.h"
#include "Texture/BlinnPhong.h"
#include "Texture/Emission.h"

std::unordered_map<std::string, MaterialElement> MaterialParser::strToElem = {
        {"Kd", DIFFUSE_K},
        {"Ks", SPECULAR_K},
        {"Ka", AMBIENT_K},
        {"Od", DIFFUSE_COLOR},
        {"DiffuseTexture", DIFFUSE_TEXTURE},
        {"Os", SPECULAR_COLOR},
        {"Kgls", GLS_K},
        {"Refl", REFLECTION_K},
        {"rJitter", REFLECTION_JITTER},
        {"Kt", TRANSMISSION_K},
        {"kt", TRANSMISSION_K},
        {"tJitter", TRANSMISSION_JITTER},
        {"ior", IOR}
};
std::unordered_map<MaterialElement, std::string> MaterialParser::elemToString = {
        {DIFFUSE_K, "Kd"},
        {SPECULAR_K, "Ks"},
        {AMBIENT_K, "Ka"},
        {DIFFUSE_COLOR, "Od"},
        {SPECULAR_COLOR, "Os"},
        {GLS_K, "Kgls"},
        {REFLECTION_K, "Refl"},
        {REFLECTION_JITTER, "rJitter"},
        {TRANSMISSION_K, "Kt"},
        {TRANSMISSION_K, "kt"},
        {TRANSMISSION_JITTER, "tJitter"},
        {IOR, "ior"}
};

bool MaterialParser::findElem(std::string str, MaterialElement& elem) {
    auto result = strToElem.find(str);
    if(result == strToElem.end()) {
        return false;
    }

    elem = result->second;
    return true;
}

SceneParser::SceneParser() {
    strToElement["CameraLookAt"] = CAMERA_LOOK_AT;
    strToElement["CameraLookFrom"] = CAMERA_LOOK_FROM;
    strToElement["CameraLookUp"] = CAMERA_LOOK_UP;
    strToElement["FieldOfView"] = FIELD_OF_VIEW;
    strToElement["AmbientLight"] = AMBIENT_LIGHT;
    strToElement["DirectionalLight"] = DIRECTIONAL_LIGHT;
    strToElement["PointLight"] = POINT_LIGHT;
    strToElement["AreaLight"] = AREA_LIGHT;
    strToElement["BackgroundColor"] = BACKGROUND_COLOR;
    strToElement["BackgroundTexture"] = BACKGROUND_TEXTURE;
    strToElement["Sphere"] = SPHERE;
    strToElement["Triangle"] = TRIANGLE;
    strToElement["Lambertian"] = LAMBERTIAN;
    strToElement["Phong"] = PHONG;
    strToElement["Reflection"] = REFLECTION;
    strToElement["Transmission"] = TRANSMISSION;
    strToElement["FresnelTransmission"] = FRESNEL_TRANSMISSION;
    strToElement["BlinnPhong"] = BLINN_PHONG;
    strToElement["Emission"] = EMISSION;
    strToElement["Material"] = MATERIAL;
    strToElement["RefractiveMaterial"] = REFRACTIVE_MATERIAL;
    strToElement["Cylinder"] = CYLINDER;

    elemToStr[CAMERA_LOOK_AT] = "CameraLookAt";
    elemToStr[CAMERA_LOOK_FROM] = "CameraLookFrom";
    elemToStr[CAMERA_LOOK_UP] = "CameraLookUp";
    elemToStr[FIELD_OF_VIEW] = "FieldOfView";
    elemToStr[AMBIENT_LIGHT] = "AmbientLight";
    elemToStr[DIRECTIONAL_LIGHT] = "DirectionalLight";
    elemToStr[POINT_LIGHT] = "PointLight";
    elemToStr[AREA_LIGHT] = "AreaLight";
    elemToStr[BACKGROUND_COLOR] = "BackgroundColor";
    elemToStr[BACKGROUND_TEXTURE] = "BackgroundTexture";
    elemToStr[SPHERE] = "Sphere";
    elemToStr[TRIANGLE] = "Triangle";
    elemToStr[LAMBERTIAN] = "Lambertian";
    elemToStr[PHONG] = "Phong";
    elemToStr[REFLECTION] = "Reflection";
    elemToStr[TRANSMISSION] = "Transmission";
    elemToStr[FRESNEL_TRANSMISSION] = "FresnelTransmission";
    elemToStr[BLINN_PHONG] = "BlinnPhong";
    elemToStr[EMISSION] = "Emission";
    elemToStr[MATERIAL] = "Material";
    elemToStr[REFRACTIVE_MATERIAL] = "RefractiveMaterial";
    elemToStr[CYLINDER] = "Cylinder";
}

Eigen::Vector2d SceneParser::readInVec2(std::ifstream& infile) {
    double x, y;
    infile >> x;
    infile >> y;
    return {x, y};
}

Eigen::Vector2d SceneParser::readInVec2(std::stringstream& instream) {
    double x, y;
    instream >> x;
    instream >> y;
    return {x, y};
}

Eigen::Vector3d SceneParser::readInVector(std::ifstream& infile){
    double x; double y; double z;
    infile >> x;
    infile >> y;
    infile >> z;
    return {x, y, z};
}

Eigen::Vector3d SceneParser::readInVector(std::stringstream &instream) {
    double x, y, z;
    instream >> x;
    instream >> y;
    instream >> z;
    return {x, y, z};
}

Scene SceneParser::parseFile(const std::string& input_file_path) {
    std::ifstream infile;
    infile.open(input_file_path);
    if(!infile.is_open()) {
        std::cerr << "failed to open input file" << std::endl;
    }

    SceneElement requiredElem[] = {CAMERA_LOOK_AT, CAMERA_LOOK_FROM, CAMERA_LOOK_UP,
                                 FIELD_OF_VIEW, AMBIENT_LIGHT, BACKGROUND};
    std::unordered_map<SceneElement, bool> requiredElemFound;
    for(SceneElement elem : requiredElem) {
        requiredElemFound[elem] = false;
    }

    std::string description;
    SceneElement elem;

    Eigen::Vector3d camera_look_at;
    Eigen::Vector3d camera_look_from;
    Eigen::Vector3d camera_look_up;
    double fov;
    Eigen::Vector3d ambient_light;
    std::shared_ptr<Texture> background;
    // Eigen::Vector3d background_color;
    std::vector<std::shared_ptr<Geometry>> geo;
    std::vector<std::shared_ptr<Light>> lights;
    std::vector<std::shared_ptr<BRDF>> brdfs;

    std::vector<uint> emissiveIndexes;
    std::vector<std::shared_ptr<Texture>> emissionTextures;

    std::string line;
    std::stringstream ss;
    while(getline(infile, line)) {
        if(line.empty()) {
            continue;
        }
        // clear string stream
        ss.clear();
        ss.str(std::string());
        // Load the line into string stream
        ss << line;
        ss >> description;

        // check if the line is a comment (starts with #)
        if(description == "#") {
            continue;
        }

        // make sure the description is a recognized scene element
        if(!strToElement.count(description)) {
            throw std::invalid_argument("didn't recognize " + description);
        }

        // Get the scene element to parse
        elem = strToElement[description];
        if(requiredElemFound.count(elem)) {
            requiredElemFound[elem] = true;
        }
        // Parse the element
        switch(elem) {
            case CAMERA_LOOK_AT: {
                camera_look_at  = readInVector(ss);
                break;
            }
            case CAMERA_LOOK_FROM: {
                camera_look_from = readInVector(ss);
                break;
            }
            case CAMERA_LOOK_UP: {
                camera_look_up = readInVector(ss);
                break;
            }
            case FIELD_OF_VIEW: {
                ss >> fov;
                break;
            }
            case AMBIENT_LIGHT: {
                ambient_light = readInVector(ss);
                break;
            }
            case DIRECTIONAL_LIGHT: {
                lights.push_back(std::move(readInDirectionalLight(infile)));
                break;
            }
            case POINT_LIGHT: {
                lights.push_back(readInPointLight(infile));
                break;
            }
            case AREA_LIGHT: {
                lights.push_back(readInAreaLight(infile));
                break;
            }
            case BACKGROUND_COLOR: {
                // background_color = readInVector(ss);
                background = std::make_shared<ColorTexture>(readInVector(ss));
                requiredElemFound[BACKGROUND] = true;
                break;
            }
            case BACKGROUND_TEXTURE: {
                std::string filename;
                ss >> filename;
                background = std::make_shared<ImageTexture>(filename);
                lights.push_back(std::make_shared<Skylight>(background));
                requiredElemFound[BACKGROUND] = true;
                break;
            }
            case LAMBERTIAN: {
                brdfs.push_back(std::move(readInLambertian(infile)));
                break;
            }
            case PHONG: {
                brdfs.push_back(std::move(readInPhong(infile)));
                break;
            }
            case REFLECTION: {
                brdfs.push_back(std::move(readInReflection(infile)));
                break;
            }
            case TRANSMISSION: {
                brdfs.push_back(std::move(readInTransmission(infile)));
                break;
            }
            case FRESNEL_TRANSMISSION: {
                brdfs.push_back(std::move(readInFresnelTransmission(infile)));
                break;
            }
            case BLINN_PHONG: {
                brdfs.push_back(std::move(readInBlinnPhong(infile)));
                break;
            }
            case EMISSION: {
                std::shared_ptr<Emission> emission = readInEmission(infile);
                brdfs.push_back(emission);
                emissiveIndexes.push_back(brdfs.size() - 1);
                emissionTextures.push_back(emission->getTexture());
                break;
            }
            case MATERIAL: case REFRACTIVE_MATERIAL: {
                // materials.push_back(std::move(readInMaterial(infile)));
                break;
            }
            case SPHERE: {
                std::string obj_description;
                geo.push_back(std::move(readInSphere(obj_description, infile, brdfs,
                                                     emissiveIndexes, emissionTextures, lights)));
                break;
            }
            case TRIANGLE: {
                std::string obj_description;
                geo.push_back(std::move(readInTriangle(obj_description, infile, brdfs,
                                                       emissiveIndexes, emissionTextures, lights)));
                break;
            }
            case CYLINDER: {
                std::string obj_description;
                geo.push_back(std::move(readInCylinder(obj_description, infile, brdfs,
                                                       emissiveIndexes, emissionTextures, lights)));
                break;
            }
            default:
                break;
        }
    }

    infile.close();

    // Make sure all of the required elements were input
    for(SceneElement rElem : requiredElem) {
        if(!requiredElemFound[rElem]) {
            throw std::runtime_error("Invalid input: " + elemToStr[rElem] + " missing.");
        }
    }

    Camera camera = Camera(camera_look_at, camera_look_from, camera_look_up, fov);

    // create scene
    Scene scene = Scene(camera, ambient_light, background);
    for(auto & g : geo) {
        scene.geo.push_back(std::move(g));
    }

    for(auto & light : lights) {
        scene.lights.push_back(std::move(light));
    }

    scene.process();
    return scene;
}

std::shared_ptr<Geometry> SceneParser::readInSphere(const std::string& obj_description, std::ifstream& infile,
                                                    const std::vector<std::shared_ptr<BRDF>>& brdfs,
                                                    const std::vector<uint>& emissiveIndexes,
                                                    const std::vector<std::shared_ptr<Texture>>& emissionTextures,
                                                    std::vector<std::shared_ptr<Light>>& lights) {
    std::string description;
    infile >> description;
    Eigen::Vector3d center = readInVector(infile);

    infile >> description;
    double radius;
    infile >> radius;

    infile >> description;
    int brdfIdx;
    infile >> brdfIdx;

    if(brdfIdx < 0 || brdfIdx >= brdfs.size()) {
        std::cerr << "Material index " << brdfIdx << " invalid. Defaulting to 0." << std::endl;
        brdfIdx = 0;
    }

    // create sphere
    std::shared_ptr<Sphere> sphere = std::make_shared<Sphere>(center, radius, brdfs[brdfIdx], obj_description);

    auto it = std::find(emissiveIndexes.begin(), emissiveIndexes.end(), brdfIdx);
    if(it != emissiveIndexes.end()) {
        // brdf is emissive
        auto index = std::distance(emissiveIndexes.begin(), it);
        lights.push_back(std::make_shared<EmissionLight>(sphere, emissionTextures[index]));
    }

    return sphere;
}

std::shared_ptr<Geometry> SceneParser::readInTriangle(const std::string& obj_description, std::ifstream& infile,
                                                      const std::vector<std::shared_ptr<BRDF>>& brdfs,
                                                      const std::vector<uint>& emissiveIndexes,
                                                      const std::vector<std::shared_ptr<Texture>>& emissionTextures,
                                                      std::vector<std::shared_ptr<Light>>& lights) {
    std::string description;
    
    // take in the vertices
    std::vector<Eigen::Vector3d> vertices;
    for(int i = 0; i < 3; ++i) {
        Eigen::Vector3d vertex = readInVector(infile);
        vertices.push_back(vertex);
    }

    bool hasSpecifiedUV = false;

    // take in UVs if present
    std::vector<Eigen::Vector2d> uvs;

    infile >> description;
    if(description == "uv") {
        hasSpecifiedUV = true;

        uvs.push_back(readInVec2(infile));
        uvs.push_back(readInVec2(infile));
        uvs.push_back(readInVec2(infile));

        infile >> description;
    }

    // BRDF index
    int brdfIdx;
    infile >> brdfIdx;

    if(brdfIdx < 0 || brdfIdx >= brdfs.size()) {
        std::cerr << "BRDF index " << brdfIdx << " invalid. Defaulting to 0." << std::endl;
        brdfIdx = 0;
    }


    std::shared_ptr<Triangle> triangle;
    // create triangle
    if(hasSpecifiedUV) {
        triangle = std::make_shared<Triangle>(vertices, uvs, brdfs[brdfIdx], obj_description);
    }
    triangle = std::make_shared<Triangle>(vertices, brdfs[brdfIdx], obj_description);

    // check emission
    auto it = std::find(emissiveIndexes.begin(), emissiveIndexes.end(), brdfIdx);
    if(it != emissiveIndexes.end()) {
        // brdf is emissive
        auto index = std::distance(emissiveIndexes.begin(), it);
        lights.push_back(std::make_shared<EmissionLight>(triangle, emissionTextures[index]));
    }

    return triangle;
}

std::shared_ptr<Geometry> SceneParser::readInCylinder(const std::string& obj_description, std::ifstream& infile,
                                                      const std::vector<std::shared_ptr<BRDF>>& brdfs,
                                                      const std::vector<uint>& emissiveIndexes,
                                                      const std::vector<std::shared_ptr<Texture>>& emissionTextures,
                                                      std::vector<std::shared_ptr<Light>>& lights) {
    std::string description;

    Eigen::Vector3d capCenter1 = readInVector(infile);
    Eigen::Vector3d capCenter2 = readInVector(infile);

    infile >> description;
    double radius;
    infile >> radius;

    infile >> description;
    int brdfIdx;
    infile >> brdfIdx;

    // check material index
    if(brdfIdx < 0 || brdfIdx > (brdfs.size() - 1)) {
        std::cerr << "Invalid material index " << brdfIdx << ". Defaulting to 0" << std::endl;
        brdfIdx = 0;
    }

    // create cylinder
    std::shared_ptr<Cylinder> cylinder = std::make_shared<Cylinder>(
            capCenter1, capCenter2, radius, brdfs[brdfIdx], obj_description);

    // check emission
    auto it = std::find(emissiveIndexes.begin(), emissiveIndexes.end(), brdfIdx);
    if(it != emissiveIndexes.end()) {
        // brdf is emissive
        auto index = std::distance(emissiveIndexes.begin(), it);
        lights.push_back(std::make_shared<EmissionLight>(cylinder, emissionTextures[index]));
    }

    return cylinder;
}

std::shared_ptr<Light> SceneParser::readInDirectionalLight(std::ifstream& infile) {
    std::string description;
    infile >> description;
    Eigen::Vector3d light_color = readInVector(infile);
    infile >> description;
    Eigen::Vector3d to_light = readInVector(infile);

    return std::make_shared<DirectionalLight>(light_color, to_light);
}

std::shared_ptr<Light> SceneParser::readInPointLight(std::ifstream& infile) {
    std::string description;
    infile >> description;
    Eigen::Vector3d light_color = readInVector(infile);
    infile >> description;
    Eigen::Vector3d position = readInVector(infile);

    return std::make_shared<PointLight>(light_color, position);
}

std::shared_ptr<Light> SceneParser::readInAreaLight(std::ifstream& infile) {
    std::string description;

    infile >> description;
    Eigen::Vector3d light_color = readInVector(infile);

    infile >> description;
    Eigen::Vector3d center = readInVector(infile);

    infile >> description;
    Eigen::Vector3d aim = readInVector(infile);

    infile >> description;
    Eigen::Vector3d up = readInVector(infile);

    infile >> description;
    double width;
    infile >> width;

    infile >> description;
    double height;
    infile >> height;

    infile >> description;
    double resolution;
    infile >> resolution;

    return std::make_shared<AreaLight>(light_color, center, aim, up, width, height, resolution);
}

std::shared_ptr<BRDF> SceneParser::readInLambertian(std::ifstream &infile) {
    std::string description;
    infile >> description;

    std::shared_ptr<Texture> diffuseReflectance = readInTexture(infile);
    // diffuseReflectance = readInVector(infile);

    return std::make_shared<Lambertian>(diffuseReflectance);
}

std::shared_ptr<BRDF> SceneParser::readInPhong(std::ifstream& infile) {
    std::string description;

    infile >> description;
    std::shared_ptr<Texture> specularReflectance = readInTexture(infile);
    // specularReflectance = readInVector(infile);

    double shininess;
    infile >> description;
    infile >> shininess;

    return std::make_shared<Phong>(specularReflectance, shininess);
}

std::shared_ptr<BRDF> SceneParser::readInReflection(std::ifstream& infile) {
    std::string description;
    double jitterAmt;
    infile >> description;
    infile >> jitterAmt;

    return std::make_shared<Reflection>(jitterAmt);
}

std::shared_ptr<BRDF> SceneParser::readInTransmission(std::ifstream &infile) {
    std::string description;
    double ior, jitterAmt;
    infile >> description;
    infile >> ior;
    infile >> description;
    infile >> jitterAmt;

    return std::make_shared<Transmission>(ior, jitterAmt);
}

std::shared_ptr<BRDF> SceneParser::readInFresnelTransmission(std::ifstream &infile) {
    std::string description;

    double reflAmt, ior, reflJitter, transJitter;
    infile >> description;
    infile >> reflAmt;
    infile >> description;
    infile >> ior;
    infile >> description;
    infile >> reflJitter;
    infile >> description;
    infile >> transJitter;

    return std::make_shared<FresnelTransmission>(reflAmt, ior, reflJitter, transJitter);
}

std::shared_ptr<BRDF> SceneParser::readInBlinnPhong(std::ifstream& infile) {
    std::string description;
    infile >> description;

    std::shared_ptr<Texture> diffuseReflectance = readInTexture(infile);
    // Eigen::Vector3d diffuseReflectance = readInVector(infile);

    infile >> description;
    // Eigen::Vector3d specularReflectance = readInVector(infile);
    std::shared_ptr<Texture> specularReflectance = readInTexture(infile);

    double shininess;
    infile >> description;
    infile >> shininess;

    // infile >> description;
    // std::shared_ptr<Texture> ambientTexture = readInTexture(infile);
    // Eigen::Vector3d ambientColor = readInVector(infile);

    double diffuseK, specularK, ambientK;
    infile >> description;
    infile >> diffuseK;
    infile >> description;
    infile >> specularK;
    infile >> description;
    infile >> ambientK;

    return std::make_shared<BlinnPhong>(diffuseReflectance, specularReflectance, diffuseReflectance, shininess,
                                        diffuseK, specularK, ambientK);
}

std::shared_ptr<Emission> SceneParser::readInEmission(std::ifstream &infile) {
    std::shared_ptr<Texture> texture = readInTexture(infile);
    return std::make_shared<Emission>(texture);
}

std::shared_ptr<Material> SceneParser::readInMaterial(std::ifstream& infile) {

    std::vector<MaterialElement> required = {DIFFUSE_K, SPECULAR_K, AMBIENT_K,
                                  DIFFUSE, SPECULAR_COLOR,
                                  GLS_K, REFLECTION_K};
    std::vector<MaterialElement> included;

    double kd, ks, ka, kgls, refl, ior;
    std::shared_ptr<Texture> diffuse;
    Eigen::Vector3d specularColor;
    // variables with default values
    double rJitter = 0;
    double kRefraction = 0.0;
    double tJitter = 0;

    std::string identifier;
    std::string line;
    std::stringstream ss;

    while(true) {
        getline(infile, line);
        // if you find a blank line, stop reading
        if(line.empty()) {
            break;
        }

        // Load the line into a stream
        ss.clear();
        ss << line;
        ss >> identifier;

        // Find the element associated with the identifier in the input
        MaterialElement elem;
        bool valid = MaterialParser::findElem(identifier, elem);
        if(!valid) {
            throw std::invalid_argument("invalid material specifier: " + identifier);
        }

        included.push_back(elem);
        // Parse the element
        switch(elem) {
            case DIFFUSE_K : {
                ss >> kd;
                break;
            }
            case SPECULAR_K : {
                ss >> ks;
                break;
            }
            case AMBIENT_K : {
                ss >> ka;
                break;
            }
            case DIFFUSE_COLOR : {
                included.push_back(DIFFUSE);
                Eigen::Vector3d diffuseColor = readInVector(ss);
                diffuse = std::make_shared<ColorTexture>(diffuseColor);
                break;
            }
            case DIFFUSE_TEXTURE : {
                included.push_back(DIFFUSE);
                std::string filename;
                ss >> filename;
                diffuse = std::make_shared<ImageTexture>(filename);
                break;
            }
            case SPECULAR_COLOR : {
                specularColor = readInVector(ss);
                break;
            }
            case GLS_K : {
                ss >> kgls;
                break;
            }
            case REFLECTION_K : {
                ss >> refl;
                break;
            }
            case REFLECTION_JITTER : {
                ss >> rJitter;
                break;
            }
            case TRANSMISSION_K : {
                ss >> kRefraction;
                break;
            }
            case TRANSMISSION_JITTER : {
                ss >> tJitter;
                break;
            }
            case IOR : {
                ss >> ior;
                break;
            }
            default : {
                break;
            }
        };
    }

    if(kRefraction > 0.0) {
        required.push_back(IOR);
    }

    // Make sure all required elements were included
    std::string message;
    if(!hasAllRequired(required, included, message)) {
        throw std::invalid_argument("Missing required material information." + message);
    }

    // create material
    return std::make_shared<Material>(kd, ks, ka, kgls, diffuse, specularColor,
                                      refl, rJitter, tJitter, ior, kRefraction);
}

std::shared_ptr<Texture> SceneParser::readInTexture(std::ifstream& infile) {
    std::string description;
    infile >> description;
    if(description == "Color") {
        Eigen::Vector3d color = readInVector(infile);
        return std::make_shared<ColorTexture>(color);
    }
    else if(description == "Image") {
        std::string filename;
        infile >> filename;
        return std::make_shared<ImageTexture>(filename);
    }
    else {
        throw std::runtime_error("texture specification not recognized " + description);
    }
}

bool SceneParser::hasAllRequired(std::vector<MaterialElement> required, std::vector<MaterialElement> included,
                                 std::string& message) {
    for(MaterialElement r : required) {
        if(std::find(included.begin(), included.end(), r) == included.end()) {
            message = "Missing " + MaterialParser::elemToString[r];
            return false;
        }
    }
    return true;
}
