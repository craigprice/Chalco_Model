
#include "ClassicalSpin3D.h"
#include <cmath>
#include <iostream>
#include <iomanip>

ClassicalSpin3D::ClassicalSpin3D(){
    x = 0;
    y = 0;
    z = 0;
    xOld = 0;
    yOld = 0;
    zOld = 0;
}
ClassicalSpin3D::ClassicalSpin3D(double x_, double y_, double z_){
    x = 0;
    y = 0;
    z = z_;
    xOld = 0;
    yOld = 0;
    zOld = 0;
}
void ClassicalSpin3D::print() const{
    std::cout << "Spin:" << std::endl;
    std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
    checkSpin();
}

void ClassicalSpin3D::doC2Rotation001() {
    x = 0;
    y = 0;
    z = z;
    xOld = 0;
    yOld = 0;
    zOld = zOld;
}

void ClassicalSpin3D::checkSpin() const{   
    if(
       (z != z)||
       ((fabs(sqrt(x * x + y * y + z * z) - 1)) >= 0.000001) ){
        std::cerr << "Error: Bad Spin values: " << "x: " << x << std::endl;
        std::cerr << "Error: Bad Spin values: " << "y: " << y << std::endl;
        std::cerr << "Error: Bad Spin values: " << "z: " << z << std::endl;
        std::cerr << "Error: Bad Spin values: " << "Magnitude: " << 
        std::setprecision(15) << std::setw(15) <<
        sqrt(x * x + y * y + z * z) << std::endl;
        exit(1);
    }
}

void ClassicalSpin3D::setRandomOrientation(){
    /**
     *Note that this is currently implemented using cartesian coordinates for
     *the spin vector. This is not necessary, it could be done in spherical, and
     *it may be faster to perform that way. This should be done as a comparison.
     *In fact I am positive that it would be faster to do this all in spherical
     *coordinates. However, currently, the "rotate" function (correctly) rotates
     *the spin vectors but it requires cartesian coordinates to work. Until I
     *find a rotation function that works in spherical, we are stuck with
     *cartesian.
     */
    x = 0;
    y = 0;
    
    if(drand48() < 0.5){
        z = 1;
    }else{
        z = -1;
    }
    
    checkSpin();
}

/**
 *First we construct a new spin perpendicular to the initial configuration by
 *adding PI/2 to the zenith angle. Then, this new vector is rotated, within the
 *same plane, anywhere in 2*PI radians. Lastly, the origial vector is rotated
 *around the perp spin by the angle "range"
 */
void ClassicalSpin3D::flip(double range){
    
    xOld = 0;
    yOld = 0;
    zOld = z;
    
    /*
    //Converts to spherical (in order to construct a perpendicular vector).
    //measured from the vertical.
    double theta = acos(z);
    double phi = atan2(y, x);
    
    double thetaPerp = theta + (PI / 2.0);
    double phiPerp = phi;
    if((thetaPerp > PI) && (thetaPerp < 2 * PI)){
        thetaPerp = PI - (thetaPerp - PI);
        phiPerp = phi + PI;
    }else if (thetaPerp >= 2 * PI){
        thetaPerp = thetaPerp - 2 * PI;
    }
    if(phiPerp > 2 * PI){
        phiPerp = phiPerp - 2 * PI;
    }
    
    //Converts back to cartesian and creates the perp spin.
    ClassicalSpin3D normalSpin(sin(thetaPerp) * cos(phiPerp),
                               sin(thetaPerp) * sin(phiPerp),
                               cos(thetaPerp));
    
    /**
     *Randomizes the direction that the perp spin is pionting. This will mean
     *that the original spin with end up rotating in a random direction
     *(within its range).
     *//*
    normalSpin.rotate(theta, phi, 2 * PI * drand48());
    
    //Once again convert to spherical in order to get the theta and phi values.
    theta = acos(normalSpin.getZ());
    phi = atan2(normalSpin.getY(), normalSpin.getX());
    
    //Rotation of the original spin within a range dictated by "range".
    this -> rotate(theta, phi, range * drand48());*/
    z = (-1) * z;
    x = 0;
    y = 0;
    
}
void ClassicalSpin3D::specifyRotation(double rotPhi, double rotTheta){
    /*xOld = x;
    yOld = y;
    zOld = z;
    
    //Converts to spherical (in order to construct a perpendicular vector).
    //measured from the vertical.
    double theta = acos(z);
    double phi = atan2(y, x);
    
    theta = theta + rotTheta;
    phi = phi + rotPhi;
    
    if((theta > PI) && (theta < 2 * PI)){
        theta = PI - (theta - PI);
        phi = phi + PI;//new
    }else if (theta >= 2 * PI){
        theta = theta - 2 * PI;
    }
    
    if(phi > 2 * PI){
        phi = phi - 2 * PI;
    }
    if(phi > 2 * PI){
        phi = phi - 2 * PI;
    }
    
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);
    z = cos(theta);
    */
    checkSpin();
    
}
void ClassicalSpin3D::reset(){
    x = 0;
    y = 0;
    z = zOld;
}
void ClassicalSpin3D::clear(){
    x = 0;
    y = 0;
    z = 0;
    xOld = 0;
    yOld = 0;
    zOld = 0;
}
/**Does rotation using quaternions. Requires that the spin has cartesian
 *components but that the axis that it rotates around is in spherical.
 *
 *Note: this function should be replaced with a function that only uses
 *spherical coordinates or only cartesian coordinates.
 */
void ClassicalSpin3D::rotate(double theta_, double phi_, double Beta) {
	/*double M[3][3];
	double q0, q1, q2, q3;
    double newS[3] = {0};
	q0 = cos(Beta / 2.0);
	q1 = sin(Beta / 2.0) * sin(theta_) * cos(phi_);
	q2 = sin(Beta / 2.0) * sin(theta_) * sin(phi_);
	q3 = sin(Beta / 2.0) * cos(theta_);
	M[0][0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
	M[0][1] = 2.0 * (q1 * q2 - q0 * q3);
	M[0][2] = 2.0 * (q1 * q3 + q0 * q2);
	M[1][0] = 2.0 * (q2 * q1 + q0 * q3);
	M[1][1] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
	M[1][2] = 2.0 * (q2 * q3 - q0 * q1);
	M[2][0] = 2.0 * (q3 * q1 - q0 * q2);
	M[2][1] = 2.0 * (q3 * q2 + q0 * q1);
	M[2][2] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
    
    /*
     for(int i = 0; i < 3; i ++){
     for(int j = 0; j < 3; j++){
     newS[i] += M[i][j] * components[j];
     }
     }
     */
    /*
    for(int i = 0; i < 3; i ++){
        for(int j = 0; j < 3; j++){
            if(i == 0 && j == 0){
                newS[i] += M[i][j] * x;
            }else if(i == 0 && j == 1){
                newS[i] += M[i][j] * y;
            }else if(i == 0 && j == 2){
                newS[i] += M[i][j] * z;
            }else if(i == 1 && j == 0){
                newS[i] += M[i][j] * x;
            }else if(i == 1 && j == 1){
                newS[i] += M[i][j] * y;
            }else if(i == 1 && j == 2){
                newS[i] += M[i][j] * z;
            }else if(i == 2 && j == 0){
                newS[i] += M[i][j] * x;
            }else if(i == 2 && j == 1){
                newS[i] += M[i][j] * y;
            }else if(i == 2 && j == 2){
                newS[i] += M[i][j] * z;
            }else{
                std::cerr << "Bad Rotation" << std::endl;
                exit(1);
            }
        }
    }
    
    //In Cartesian
    x = newS[0];
    y = newS[1];
    z = newS[2];*/
}

