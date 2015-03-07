

/**
 ClassicalSpin3D.cpp
 Purpose: Spin object for the simulation of the lattice.
 
 @author Craig Price
 @version 1.0 2015/03/04
 */

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
    //(copied below)
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
    
    double cosr = cos(r);
    double sinr = sin(r);
    
    x = x * circDist;
    y = y * circDist;
    z = z * circDist;
    
    x = x + circleRadius * (cosr * u[0] + sinr * v[0]);
    y = y + circleRadius * (cosr * u[1] + sinr * v[1]);
    z = z + circleRadius * (cosr * u[2] + sinr * v[2]);
    
    /*
     //Should not be needed.
     size = sqrt(x*x + y*y + z*z);
     x = x / size;
     y = y / size;
     z = z / size;
     */
    
}

/*
 From: ags@seaman.cc.purdue.edu (Dave Seaman)
 Newsgroups: sci.math.num-analysis
 Subject: Re: N-dim spherical random number drawing
 Date: 20 Sep 1996 13:04:58 -0500
 
 In article <wmEgVDG00VUtMEhUwg@andrew.cmu.edu>,
 Shing-Te Li  <sl2x+@andrew.cmu.edu> wrote:
 >Does anyone know a good algorithm that can draw a vector of random
 >numbers uniformly from a N dimensional sphere?
 
 Here are four methods for generating uniformly-distributed points on a unit
 sphere:
 
 (1) The normal-deviate method.  Choose x, y, and z, each
 normally-distributed with mean 0 and variance 1.  (Actually, the
 variance does not matter, provided it remains fixed.)  Normalize
 the result to obtain a uniform density on the unit sphere.
 
 This method also generalizes nicely to n dimensions.  It works
 because the vector chosen (before normalization) has a density
 that depends only on the distance from the origin.  The method is
 mentioned in Knuth, _Seminumerical Algorithms_, 2nd. ed., pp.
 130-131.
 
 (2) The hypercube rejection method.  Choose x, y, and z, uniformly
 distributed on [-1,1].  Compute s = x^2 + y^2 + z^2.  If s > 1,
 reject the triplet and start over.  Once you have an acceptable
 vector, normalize it to obtain a point on the unit sphere.
 
 This method also generalizes to n dimensions, but as the dimension
 increases the probability of rejection becomes high and many random
 numbers are wasted.
 
 (3) The trig method.  This method works only in 3-space, but it is
 very fast.  It depends on the slightly counterintuitive fact (see
 proof below) that each of the three coordinates of a uniformly
 distributed point on S^2 is uniformly distributed on [-1,1] (but
 the three are not independent, obviously).  Therefore, it
 suffices to choose one axis (Z, say) and generate a uniformly
 distributed value on that axis.  This constrains the chosen point
 to lie on a circle parallel to the X-Y plane, and the obvious
 trig method may be used to obtain the remaining coordinates.
 
	(a) Choose z uniformly distributed in [-1,1].
	(b) Choose t uniformly distributed on [0, 2*pi).
	(c) Let r = sqrt(1-z^2).
	(d) Let x = r * cos(t).
	(e) Let y = r * sin(t).
 
 (4) A variation of method (3) that doesn't use trig may be found in
 Knuth, _loc. cit._.  The method is equivalent to the following:
 
	(a) Choose u and v uniformly distributed on [-1,1].
	(b) Let s = u^2 + v^2.
	(c) If s > 1, return to step (a).
	(d) Let a = 2 * sqrt(1-s).
	(e) The desired point is (a*u, a*v, 2*s-1).
 
 This method uses two-dimensional rejection, which gives a higher
 probability of acceptance compared to algorithm (2) in three
 dimensions.  I group this with the trig method because both
 depend on the fact that the projection on any axis is uniformly
 distributed on [-1,1], and therefore both methods are limited to
 use on S^2.
 
 I have found the trig method to be fastest in tests on an IBM RS/6000
 model 590, using the XLF 3.2 Fortran compiler, with the IMSL
 subroutines DRNUN and DRNNOA generating the random numbers.  The trig
 routines were accelerated by using the MASS library, and sqrt
 computations were performed by the hardware sqrt instruction available
 on the 590.  Timings for 1 million random points on S^2:
 
	(1) normal-deviate method                   9.60 sec.
	(2) 3-D rejection                           8.84 sec.
	(3) trig method                             2.71 sec.
	(4) polar with 2-D rejection                5.08 sec.
 
 The fact that algorithm (3) does not involve rejection turns out to be
 a huge advantage because it becomes possible to generate all the points
 in a single loop with no conditional statements, which optimizes very
 effectively and more than offsets the cost of using trig functions.
 The alternative with algorithms (2) and (4) is to generate the points
 one at a time and test each one for acceptance before proceeding to the
 next, which requires nested loops and many more subroutine calls to
 generate the random numbers.  Algorithm (1) does not explicitly use
 rejection, but the IMSL subroutine DRNNOA uses rejection in generating
 the normally-distributed numbers that are required.
 
 In higher dimensions, the normal-deviate method is preferred
 because the probability of rejection in the hypercube method
 grows very rapidly with the dimension.
 
 --------------------------------------------------------------------
 
 Here is a proof of correctness for the fact that the z-coordinate is
 uniformly distributed.  The proof also works for the x- and y-
 coordinates, but note that this works only in 3-space.
 
 The area of a surface of revolution in 3-space is given by
 
	A = 2 * pi * int_a^b f(x) * sqrt(1 + [f'(x}]^2) dx
 
 Consider a zone of a sphere of radius R.  Since we are integrating in
 the z direction, we have
 
	f(z) = sqrt(R^2 - z^2)
	f'(z) = -z / sqrt(R^2-z^2)
 
	1 + [f'(z)]^2 = r^2 / (R^2-z^2)
 
	A = 2 * pi * int_a^b sqrt(R^2-z^2) * R/sqrt(R^2-z^2) dz
 
 = 2 * pi * R int_a^b dz
 
 = 2 * pi * R * (b-a)
 
 = 2 * pi * R * h,
 
 where h = b-a is the vertical height of the zone.  Notice how the integrand
 reduces to a constant.  The density is therefore uniform.
 
 --
 Dave Seaman			dseaman@purdue.edu
 ++++ stop the execution of Mumia Abu-Jamal ++++
 ++++ if you agree copy these lines to your sig ++++
 ++++ see http://www.xs4all.nl/~tank/spg-l/sigaction.htm ++++
 ==============================================================================
 
 From: George Marsaglia <geo@stat.fsu.edu>
 Subject: Random points on a sphere.
 Date: Tue, 15 Jun 1999 18:03:07 -0400
 Newsgroups: sci.math.num-analysis,sci.math,sci.stat.math
 
 Questions on methods for random sampling from the surface
 of a sphere keep coming up, and seemingly get lost as
 new generations of computer users encounter the problem.
 
 I described the methods that I had used for years in
 "Sampling from the surface of a sphere", Ann Math Stat, v43, 1972,
 recommending the two that seemed the fastest:
 
 1) Choose standard normal variates x1,x2,x3,
 put R=1/sqrt(x1^2+x2^2+x3^2) and return the point
 (x1*R,x2*R,x3*R).
 
 2) Generate u and v, uniform in [-1,1]  until S=u^2+v^2 < 1,
 then return (1-2*S,u*r,v*r), with r=2*sqrt(1-S).
 
 The second method is fast and easiest to program---unless
 one already has a fast normal generator in his library.
 
 Method 1) seems by far the fastest (and easily applies to higher
 dimensions) if my latest normal generator is used.
 It produces normal variates at the rate of seven million
 per second in a 300MHz Pc.
 
 George Marsaglia
 */

/*
 Hi,
 
 if we have two sphere, S1 and S2, their intersection will be a circle
 C (suppose that the intersection of the two spheres exist and it is
 not a point). On S1 and S2 we will have a two spherical caps. C will
 be the border of the spherical cap.
 
 There is an efficient way to generate random points on one of the
 spherical cap ?
 
 One way that should work well, if the radius of the circle is small,
 is to generate a point p on the corresponding disk and to compute the
 intersection of the ray o-p with the sphere S1 (for example), where o
 is the center of S1.
 
 But if the radius of the circle is not small, the disk will not be a
 good approximation of the spherical cap, so the points will not be
 uniformly distributed on the cap...
 
 Thank you
 Luca
 Reply With Quote Reply With Quote
 12-11-2007 09:39 PM #2
 Default Re: Random points on a spherical cap
 <kingpin@freemail.it> wrote in message
 news:96bb8ce2-0d44-4b94-b5ec-9ee829e485e7@p69g2000hsa.googlegroups.com...
 
 To generate random points on the spherical cap for S1, let K1 be
 the sphere center and let R1 be the sphere radius. Let A be the
 circle center and let S be the circle radius. Compute
 W = (A-K1)/Length(A-K1). Choose U and V to be unit-length
 vectors so that {U,V,W} is an orthonormal set (mutually
 perpendicular, all unit-length).
 
 Select a random angle t in [0,2*pi) and a random length L in [0,S].
 The point P = A + L*(cos(t)*U + sin(t)*V) lies inside the circle.
 Define D = P - K1. The corresponding point on the spherical cap is
 Q = R1*D/Length(D).
 
 --
 Dave Eberly
 http://www.geometrictools.com
 
 */
