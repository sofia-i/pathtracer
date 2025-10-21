//
//  Camera.h
//  raytracer_2
//
//  Created by Sofia Iannicelli on 2/15/23.
//

#ifndef Camera_h
#define Camera_h

#include <Eigen/Eigen>
#include <cmath>

// #define PI 3.14159265

class Camera {
public:
    Camera(Eigen::Vector3d lookAt, Eigen::Vector3d lookFrom, Eigen::Vector3d lookUp, double fov) :
            lookAt(lookAt), lookFrom(lookFrom), lookUp(lookUp),
            fovX(fov), fovY(fov) { }
    
    Camera(Eigen::Vector3d lookAt,Eigen::Vector3d lookFrom, Eigen::Vector3d lookUp, double fovX, double fovY) :
        lookAt(lookAt), lookFrom(lookFrom), lookUp(lookUp), fovX(fovX), fovY(fovY) { }

    Camera() : Camera(
            Eigen::Vector3d(0.0, 0.0, 0.0),  // at
            Eigen::Vector3d(0.0, 0.0, 1.0),  // from
            Eigen::Vector3d(0.0, 1.0, 0.0),  // up
            90.0) {}

    void setLookAt(Eigen::Vector3d newLookAt) { this->lookAt = newLookAt; }
    Eigen::Vector3d getLookAt() { return lookAt; }

    void setLookFrom(Eigen::Vector3d newLookFrom) { this->lookFrom = newLookFrom; }
    Eigen::Vector3d getLookFrom() { return lookFrom; }

    void setCameraLookUp(Eigen::Vector3d newLookUp) { this->lookUp = newLookUp; }
    Eigen::Vector3d getCameraLookUp() { return lookUp; }

    void setFovX(double newFovX) { this->fovX = newFovX; }
    void setFovY(double newFovY) { this->fovY = newFovY; }
    void updateFov(double ratioX, double ratioY) {
        fovX = (ratioX / ratioY) * fovY;
    }
    double getFovX() { return fovX; }
    double getFovY() { return fovY; }
    double getFovXRad() { return fovX * M_PI / 180; }
    double getFovYRad() { return fovY * M_PI / 180; }

    friend std::ostream& operator<<(std::ostream& os, Camera const &camera) {
        os << "Camera" << std::endl;
        os << "\t" << "Looking At " << camera.lookAt[0] << " " << camera.lookAt[1] << " " <<
                camera.lookAt[2] << std::endl;
        os << "\t" << "Looking From " << camera.lookFrom[0] << " " << camera.lookFrom[1] << " " <<
                camera.lookFrom[2] << std::endl;
        os << "\t" << "Look Up " << camera.lookUp[0] << " " << camera.lookUp[1] << " " <<
                camera.lookUp[2] << std::endl;
        os << "\t" << "Field of View: x " << camera.fovX << " y " << camera.fovY << std::endl;
        return os;
    }

private:
    Eigen::Vector3d lookAt;
    Eigen::Vector3d lookFrom;
    Eigen::Vector3d lookUp;
    double fovX;
    double fovY;
};


#endif /* Camera_h */
