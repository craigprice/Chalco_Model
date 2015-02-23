
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
    x = x_;
    y = y_;
    z = z_;
    xOld = 0;
    yOld = 0;
    zOld = 0;
}
void ClassicalSpin3D::print() const{
    std::cout << "Spin:" << std::endl;
    std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
}

void ClassicalSpin3D::doC2Rotation001() {
    x = (-1)*x;
    y = (-1)*y;
    z = z;
    xOld = (-1)*xOld;
    yOld = (-1)*yOld;
    zOld = zOld;
}

void ClassicalSpin3D::checkSpin(){   
    double size = sqrt(x*x + y*y + z*z);
    x = x / size;
    y = y / size;
    z = z / size;
    if((x > 1 || x < -1)||
       (y > 1 || y < -1)||
       (z > 1 || z < -1)||
       (x != x)||
       (y != y)||
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
    x = 2 * drand48() - 1;
    y = 2 * drand48() - 1;
    z = 2 * drand48() - 1;
    double size = sqrt(x * x + y * y + z * z);
    x = x / size;
    y = y / size;
    z = z / size;
    
}


void ClassicalSpin3D::specifyRotation(double rotPhi, double rotTheta){
    xOld = x;
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

    checkSpin();
    
}
void ClassicalSpin3D::reset(){
    x = xOld;
    y = yOld;
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
	double M[3][3];
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
    z = newS[2];
}

/**
 *First we construct a new spin perpendicular to the initial configuration by
 *adding PI/2 to the zenith angle. Then, this new vector is rotated, within the
 *same plane, anywhere in 2*PI radians. Lastly, the origial vector is rotated
 *around the perp spin by the angle "range"
 */
 void ClassicalSpin3D::flip(double range){
    
    xOld = x;
    yOld = y;
    zOld = z;
    
    double u[3];
    //Cross product with a sufficiently far away "g" to find a perp vector.
    //u[0] = g2s3 - g3s2;
    //u[1] = g3s1 - g1s3;
    //u[2] = g1s2 - g2s1;
    if(x > 0.57) {//g = (0,1,0)
        u[0] = z;
        u[1] = 0;
        u[2] = (-1) * x;
    }else if(y > 0.57){//g = (0,0,1)
        u[0] = (-1) * y;
        u[1] = x;
        u[2] = 0;
    }
    else{//g = (1,0,0)
        u[0] = 0;
        u[1] = (-1) * z;
        u[2] = y;
    }
    
    //Cross product to find third Orthogonal vector
    double v[3] = { u[1]*z - u[2]*y, u[2]*x - u[0]*z, u[0]*y - u[1]*x};
    
    //Normalize vectors
    double size = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    u[0] = u[0] / size;
    u[1] = u[1] / size;
    u[2] = u[2] / size;
    
    size = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] = v[0] / size;
    v[1] = v[1] / size;
    v[2] = v[2] / size;
    
    //spin, u, v are orthonormal.
    //See these two sites:
    //http://www.math.niu.edu/~rusin/known-math/96/sph.rand
    //http://objectmix.com/graphics/314285-random-points-spherical-cap.html
    //Except that the second link needs to be modified in a number of ways to
    //apply correctly to this situation. Most notably, before adding the radial
    //vector (from a vector in the plane perp to the spin) to the spin vector,
    //we need to shorten the spin vector such that the addition of the radial
    //vector lands on the surface of the sphere. We cannot simply normalize
    //at the end like they say.
    
    //Size of circle:
    if(range >= PI){range = PI;}
    double r = drand48() * 2.0 * PI;
    double distFromOrig = cos(range);//works both for pos dist and neg dist.
    double circDist = drand48() * (1 - distFromOrig) + distFromOrig;
    double circleRadius = sqrt(1 - circDist*circDist);
    
    x = x * circDist;
    y = y * circDist;
    z = z * circDist;
    
    x = x + circleRadius * (cos(r) * u[0] + sin(r) * v[0]);
    y = y + circleRadius * (cos(r) * u[1] + sin(r) * v[1]);
    z = z + circleRadius * (cos(r) * u[2] + sin(r) * v[2]);
    
    /*
     //Should not be needed.
     size = sqrt(x*x + y*y + z*z);
     x = x / size;
     y = y / size;
     z = z / size;
     */
    
}

