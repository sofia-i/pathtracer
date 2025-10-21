cmake_minimum_required(VERSION 3.29)
project(pathtracer)

set(CMAKE_CXX_STANDARD 20)

# Find Eigen path to include
find_path(EIGEN_DIR NAMES signature_of_eigen3_matrix_library
        PATHS
        /opt/homebrew/Cellar/eigen
        #        /usr/local/include/eigen3
        PATH_SUFFIXES eigen3 eigen)
if(NOT EXISTS ${EIGEN_DIR})
    message(FATAL_ERROR "Please check that the set Eigen directory is valid!")
endif()
include_directories(${EIGEN_DIR})

# openmp
find_package(OpenMP REQUIRED)
if (OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(pathtracer main.cpp
        src/BoundingVolume/BoundingBox.cpp
        src/BoundingVolume/BoundingVolumeHierarchy.cpp
        src/Camera.h
        src/Geometry/Cylinder.cpp
        src/Geometry/Extent.h
        src/Geometry/Geometry.h
        src/Geometry/Hittable.h
        src/Image.cpp
        src/ImageUtils.cpp
        src/Texture/Material.cpp
        src/Geometry/Mesh.cpp
        src/Pathtracer.cpp
        src/PpmWriter.hpp
        src/Ray.hpp
        src/Scene.hpp
        src/SceneParser.cpp
        src/Light/Skylight.h
        src/Geometry/Sphere.cpp
        src/Texture/Texture.cpp
        src/Geometry/Triangle.cpp
        src/MathUtils/Utils.h
        src/Texture/BRDF.h
        src/PngWriter.h
        src/MathUtils/Sampler.cpp
        src/MathUtils/Sampler.h
        src/Texture/Lambertian.h
        src/Texture/Phong.h
        src/Texture/Reflection.h
        src/Texture/Transmission.h
        src/Texture/BlinnPhong.h
        src/Texture/FresnelTransmission.h
        src/Texture/Emission.h
        src/Light/DirectionalLight.h
        src/Light/PointLight.h
        src/Light/AreaLight.h
        src/Light/EmissionLight.h
        # include/stb_image.h
        )

if (OpenMP_FOUND)
    target_link_libraries(pathtracer PUBLIC OpenMP::OpenMP_CXX)
endif()

# link_directories()
