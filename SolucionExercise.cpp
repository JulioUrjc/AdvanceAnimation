#include <iostream>
using namespace std;
#include "Vec2.h"
#include "LinSys.h"
#include "Exercise.h"

// gravitational acceleration (9.81)
const float g = 9.81f;

// Linear system example
void LinearSystemExample(void)
{
   //Set up matrix A and vector b
   int size=3;
   MatrixMN A(size, size);
   A(0,0)=1.0f; A(0,1)=2.0f; A(0,2)=3.0f;
   A(1,0)=4.0f; A(1,1)=5.0f; A(1,2)=6.0f;
   A(2,0)=7.0f; A(2,1)=8.0f; A(2,2)=9.0f;
   Vector b(size);
   b(0)=10.0f; b(1)=11.0f; b(2)=12.0f;
   Vector x(size);

   //Solve linear system A*x=b
   GaussElimination(A, b, x);
}

void SetSubmatrix(MatrixMN& M, int row, int col, const Matrix2& subM)
{
   M(row, col) = subM.v[0][0];
   M(row, col+1) = subM.v[0][1];
   M(row+1, col) = subM.v[1][0];
   M(row+1, col+1) = subM.v[1][1];
}

void SetSubvector(Vector& v, int row, const Vec2& subv)
{
   v(row) = subv.x;
   v(row+1) = subv.y;
}

void GetSubvector(const Vector& v, int row, Vec2& subv)
{
   subv.x = v(row);
   subv.y = v(row+1);
}

class Spring1D
{
public:
   float p1;
   float p2;
   float L;
   float k;

   Spring1D(void){}
   Spring1D(float p1, float p2, float L, float k) : p1(p1), p2(p2), L(L), k(k) {}
   ~Spring1D(void){}

   float Force1(void) const
   {
      return -k*(p1-p2-L);
   }
   float Force2(void) const
   {
      return -k*(p2-p1+L);
   }

   float dF1dp1(void) const
   {
      return -k;
   }
   float dF1dp2(void) const
   {
      return k;
   }
   float dF2dp1(void) const
   {
      return k;
   }
   float dF2dp2(void) const
   {
      return -k;
   }

};


class Spring
{
public:
   Vec2 p1;
   Vec2 p2;
   float L0;
   float L;
   float k;
   Vec2 u1;
   Vec2 F;
   Matrix2 dFdp;

   Spring(void){}
   Spring(const Vec2& p1, const Vec2& p2, float L0, float k) : p1(p1), p2(p2), L0(L0), k(k)
   {
      L = (p1-p2).length();
      u1 = (1.0f/L) * (p1 - p2);
      F = (k*(L0-L)) * u1;

      //u1 = (p1 - p2)/l
      //l^2 = (p1 - p2)^T (p1 - p2)
      //2 l dl/dp1 = 2 (p1 - p2)^T
      //dl/dp1 = u1^T

      //l * u1 = (p1 - p2)
      //l * du1/dp1 + u1 * dl/dp1 = I
      //du1/dp1 = 1/l * (I - u1 * dl/dp1)
      //du1/dp1 = 1/l * (I - u1 * u1^T)

      //dF1/dp1 = - k * (l - l0) * du1/dp1 - k * u1 * dl/dp1
      //dF1/dp1 = - k * (l - l0) * 1/l * (I - u1 * u1^T) - k * u1 * u1^T
      //dF1/dp1 = - k * (1 - l0/l) * (I - u1 * u1^T) - k * u1 * u1^T
      //dF1/dp1 = - k * (I - u1 * u1^T) + k * l0/l * (I - u1 * u1^T) - k * u1 * u1^T
      //dF1/dp1 = - k * I + k * l0/l * (I - u1 * u1^T)
      //dF1/dp1 = - k * (1 - l0/l) * I - k * l0/l * u1 * u1^T

      dFdp = (k*(L0/L-1)) * Matrix2::IDENTITY - ((k*L0/L) * u1) * u1;
   }
   ~Spring(void){}

   const Vec2& Force1(void) const
   {
      return F;
   }
   Vec2 Force2(void) const
   {
      return (-1.0)*F;
   }

   const Matrix2& dF1dp1(void) const
   {
      return dFdp;
   }
   Matrix2 dF1dp2(void) const
   {
      return (-1.0)*dFdp;
   }
   Matrix2 dF2dp1(void) const
   {
      return (-1.0)*dFdp;
   }
   const Matrix2& dF2dp2(void) const
   {
      return dFdp;
   }

};


class Triangle
{
public:
   //3 points defined in counterclockwise order: p1, p2, p3
   Vec2 p1;
   Vec2 p2;
   Vec2 p3;
   float A;
   float A0;
   float k;
   Vec2 dAdp1T, dAdp2T, dAdp3T;
   Matrix2 d2AdpTdnext, d2AdpTdprev;

   Triangle(void){}
   Triangle(const Vec2& p1, const Vec2& p2, const Vec2& p3, float A0, float k) : p1(p1), p2(p2), p3(p3), A0(A0), k(k)
   {
      //A = 1/2 * (p2 - p1) x (p3 - p1) = 1/2 * (p3 - p2) x (p1 - p2) = 1/2 * (p1 - p3) x (p2 - p3)
      A = 0.5f * (p2 - p1).cross(p3 - p1);

      //A = 1/2 * (prev - next) x (p - next)
      //dAdp = 1/2 * (prev - next)*
      dAdp1T = 0.5f * Vec2::CrossProd(p3 - p2);
      dAdp2T = 0.5f * Vec2::CrossProd(p1 - p3);
      dAdp3T = 0.5f * Vec2::CrossProd(p2 - p1);

      //u* = (-y, x)
      //du*T/du = (0, -1; 1, 0)
      d2AdpTdprev = Matrix2(0.0f, -1.0f, 1.0f, 0.0f);
      d2AdpTdnext = Matrix2(0.0f, 1.0f, -1.0f, 0.0f);
      //d2AdpTdp = 0

      //E = 1/2 * k * (A - A0)^2

      //Fi = - k * (A - A0) * dA/dpi^T

      //dFi/dpj = - k * (A - A0) * d(dA/dpi^T)/dpj - k * dA/dpi^T * dA/dpj
   }
   ~Triangle(void){}

   Vec2 Force1(void) const
   {
      
      return (k*(A0-A)) * dAdp1T;
   }
   Vec2 Force2(void) const
   {
      return (k*(A0-A)) * dAdp2T;
   }
   Vec2 Force3(void) const
   {
      return (k*(A0-A)) * dAdp3T;
   }

   Matrix2 dF1dp1(void) const
   {
      return ((-k) * dAdp1T) * dAdp1T;
   }
   Matrix2 dF1dp2(void) const
   {
      return (k*(A0-A)) * d2AdpTdnext - (k * dAdp1T) * dAdp2T;
   }
   Matrix2 dF1dp3(void) const
   {
      return (k*(A0-A)) * d2AdpTdprev - (k * dAdp1T) * dAdp3T;
   }

   Matrix2 dF2dp1(void) const
   {
      return (k*(A0-A)) * d2AdpTdprev - (k * dAdp2T) * dAdp1T;
   }
   Matrix2 dF2dp2(void) const
   {
      return ((-k) * dAdp2T) * dAdp2T;
   }
   Matrix2 dF2dp3(void) const
   {
      return (k*(A0-A)) * d2AdpTdnext - (k * dAdp2T) * dAdp3T;
   }

   Matrix2 dF3dp1(void) const
   {
      return (k*(A0-A)) * d2AdpTdnext - (k * dAdp3T) * dAdp1T;
   }
   Matrix2 dF3dp2(void) const
   {
      return (k*(A0-A)) * d2AdpTdprev - (k * dAdp3T) * dAdp2T;
   }
   Matrix2 dF3dp3(void) const
   {
      return ((-k) * dAdp3T) * dAdp3T;
   }

};


void AdvanceEuler1(float k, float m, float d, float L, float kc, float dt,
                  float p1, float v1, float& p2, float& v2, bool collision)
{
   //Initialize forces
   float mg=-m*g;
   float f2=0.0f;

   //Add weights
   f2+=mg;

   //Add spring forces
   Spring1D sp12(p1, p2, L, k);
   f2+=sp12.Force2();

   //Add damping forces
   f2-=d*v2;

   //Test collisions
   if (collision)
   {
      if (p2 < 0.0f)
      {
         f2 -= kc*p2;
      }
   }

   //Integrate positions
   p2+=dt*v2;

   //Integrate velocities
   v2+=dt*(1.0f/m)*f2;
}


void AdvanceEulerSymplec1(float k, float m, float d, float L, float kc, float dt,
                     float p1, float v1, float& p2, float& v2, bool collision)
{
   //Initialize forces
   float mg=-m*g;
   float f2=0.0f;

   //Add weights
   f2+=mg;

   //Add spring forces
   Spring1D sp12(p1, p2, L, k);
   f2+=sp12.Force2();

   //Add damping forces
   f2-=d*v2;

   //Test collisions
   if (collision)
   {
      if (p2 < 0.0f)
      {
         f2 -= kc*p2;
      }
   }

   //Integrate velocities
   v2+=dt*(1.0f/m)*f2;

   //Integrate positions
   p2+=dt*v2;
}


void AdvanceEulerSymplec(float k, float m, float d, float L, float kA, float A0, float dt,
   const Vec2& p1, const Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3, Vec2& p4, Vec2& v4, bool springs, bool area)
{
   //Initialize forces
   Vec2 mg(0.0f, -m*g);
   Vec2 f2(0.0f, 0.0f);
   Vec2 f3(0.0f, 0.0f);
   Vec2 f4(0.0f, 0.0f);

   //Add weights
   f2+=mg;
   f3+=mg;
   f4+=mg;

   //Add spring forces
   if (springs)
   {
      Spring sp12(p1, p2, L, k);
      Spring sp23(p2, p3, L, k);
      Spring sp34(p3, p4, L, k);
      Spring sp41(p4, p1, L, k);
      f2+=sp23.Force1()+sp12.Force2();
      f3+=sp34.Force1()+sp23.Force2();
      f4+=sp41.Force1()+sp34.Force2();
   }

   //Add damping forces
   f2-=d*v2;
   f3-=d*v3;
   f4-=d*v4;

   //Add triangle forces
   if (area)
   {
      Triangle tri1(p1, p2, p3, A0/2.0f, kA);
      f2+=tri1.Force2();
      f3+=tri1.Force3();

      Triangle tri2(p3, p4, p1, A0/2.0f, kA);
      f3+=tri2.Force1();
      f4+=tri2.Force2();
   }

   //Integration
   v2+=dt*(1.0f/m)*f2;
   p2+=dt*v2;
   v3+=dt*(1.0f/m)*f3;
   p3+=dt*v3;
   v4+=dt*(1.0f/m)*f4;
   p4+=dt*v4;
}


void AdvanceMidPoint1(float k, float m, float d, float L, float kc, float dt,
                     float p1, float v1, float& p2, float& v2, bool collision)
{
   float mg=-m*g;

   //a0: a at t
   Spring1D sp12(p1, p2, L, k);
   float a0 = (1.0f/m)*(mg - d*v2 + sp12.Force2() - (collision && p2 < 0.0f) ? kc*p2 : 0.0f);

   //v2_: v at t+h/2
   float v2_ = v2 + (dt/2.0f)*a0;

   //p2_: p at t+h/2
   float p2_ = p2 + (dt/2.0f)*v2;

   //a2_: a at t+h/2 with p and v at t+h/2
   sp12=Spring1D(p1, p2_, L, k);
   float a2_ = (1.0f/m)*(mg - d*v2_ + sp12.Force2() - (collision && p2_ < 0.0f) ? kc*p2_ : 0.0f);

   //compute positions
   p2+=dt*v2_;

   //compute velocities
   v2+=dt*a2_;
}


void AdvanceBackEuler1(float k, float m, float d, float L, float kc, float dt,
                      float p1, float v1, float& p2, float& v2, bool collision)
{
   //Compute forces
   Spring1D sp12(p1, p2, L, k);
   float mg = -m*g;
   float F = mg - d*v2 + sp12.Force2() - ((collision && p2 < 0.0f) ? kc*p2 : 0.0f);

   //Set up system
   float A = m + dt*d - dt*dt*(sp12.dF2dp2() - ((collision && p2 < 0.0f) ? kc : 0.0f));
   float b = (m + dt*d) * v2 + dt * F;

   //Solve for velocity
   v2 = b/A;

   //Integrate position
   p2+=dt*v2;
}


void AdvanceBackEuler(float k, float m, float d, float L, float kA, float A0, float dt,
   const Vec2& p1, const Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3, Vec2& p4, Vec2& v4, bool springs, bool area)
{
   //Initialize forces
   Vec2 mg(0.0f, -m*g);
   Vec2 F2 = mg - d*v2;
   Vec2 F3 = mg - d*v3;
   Vec2 F4 = mg - d*v4;

   //Initialize Jacobians
   Matrix2 dF2dp2 = Matrix2::ZERO, dF2dp3 = Matrix2::ZERO;
   Matrix2 dF3dp3 = Matrix2::ZERO, dF3dp2 = Matrix2::ZERO, dF3dp4 = Matrix2::ZERO;
   Matrix2 dF4dp4 = Matrix2::ZERO, dF4dp3 = Matrix2::ZERO;

   //Springs
   if (springs)
   {
      Spring sp12(p1, p2, L, k);
      Spring sp23(p2, p3, L, k);
      Spring sp34(p3, p4, L, k);
      Spring sp41(p4, p1, L, k);

      F2 += sp23.Force1() + sp12.Force2();
      F3 += sp34.Force1() + sp23.Force2();
      F4 += sp41.Force1() + sp34.Force2();

      //Collect spring Jacobians
      dF2dp2 += sp23.dF1dp1() + sp12.dF2dp2();
      dF2dp3 += sp23.dF1dp2();
      dF3dp3 += sp34.dF1dp1() + sp23.dF2dp2();
      dF3dp2 += sp23.dF2dp1();
      dF3dp4 += sp34.dF1dp2();
      dF4dp4 += sp41.dF1dp1() + sp34.dF2dp2();
      dF4dp3 += sp34.dF2dp1();
   }

   //Area preservation
   if (area)
   {
      Triangle tri1(p1, p2, p3, A0/2.0f, kA);
      F2 += tri1.Force2();
      F3 += tri1.Force3();

      Triangle tri2(p3, p4, p1, A0/2.0f, kA);
      F3 += tri2.Force1();
      F4 += tri2.Force2();

      dF2dp2 += tri1.dF2dp2();
      dF2dp3 += tri1.dF2dp3();
      dF3dp3 += tri1.dF3dp3() + tri2.dF1dp1();
      dF3dp2 += tri1.dF3dp2();
      dF3dp4 += tri2.dF1dp2();
      dF4dp4 += tri2.dF2dp2();
      dF4dp3 += tri2.dF2dp1();
   }

   //Set up A matrix
   Matrix2 A11 = (m+dt*d)*Matrix2::IDENTITY - (dt*dt)*dF2dp2;
   Matrix2 A12 = (-dt*dt)*dF2dp3;
   Matrix2 A21 = (-dt*dt)*dF3dp2;
   Matrix2 A22 = (m+dt*d)*Matrix2::IDENTITY - (dt*dt)*dF3dp3;
   Matrix2 A23 = (-dt*dt)*dF3dp4;
   Matrix2 A32 = (-dt*dt)*dF4dp3;
   Matrix2 A33 = (m+dt*d)*Matrix2::IDENTITY - (dt*dt)*dF4dp4;

   //Set up b vector
   Vec2 b1 = (m+dt*d)*v2 + dt*F2;
   Vec2 b2 = (m+dt*d)*v3 + dt*F3;
   Vec2 b3 = (m+dt*d)*v4 + dt*F4;

   //Linear system. Solve for x in A*x=b
   //Solve for new velocities
   MatrixMN A(6, 6, 0.0f);
   SetSubmatrix(A, 0, 0, A11);
   SetSubmatrix(A, 0, 2, A12);
   SetSubmatrix(A, 2, 0, A21);
   SetSubmatrix(A, 2, 2, A22);
   SetSubmatrix(A, 2, 4, A23);
   SetSubmatrix(A, 4, 2, A32);
   SetSubmatrix(A, 4, 4, A33);
   Vector b(6);
   SetSubvector(b, 0, b1);
   SetSubvector(b, 2, b2);
   SetSubvector(b, 4, b3);
   Vector x(6);
   GaussElimination(A, b, x);
   GetSubvector(x, 0, v2);
   GetSubvector(x, 2, v3);
   GetSubvector(x, 4, v4);

   //Integrate positions
   p2+=dt*v2;
   p3+=dt*v3;
   p4+=dt*v4;

}


void AdvanceVerlet1(float k, float m, float d, float L, float kc, float dt,
                     float p1, float v1, float& p2, float& v2, float& p2old, bool collision)
{
   return; // CHECK FOR THE INITIALIZATION OF THE VERLET METHOD!!!!

   //Save old position
   float aux = p2;

   //Initialize forces
   float mg=-m*g;
   float f2=0.0f;

   //Add weights
   f2+=mg;

   //Add spring forces
   Spring1D sp12(p1, p2, L, k);
   f2+=sp12.Force2();

   //Add damping forces
   f2-=d*v2;

   //Test collisions
   if (collision)
   {
      if (p2 < 0.0f)
      {
         f2 -= kc*p2;
      }
   }

   //Integrate positions
   p2=2.0f*p2-p2old+dt*dt*(1.0f/m)*f2;

   //Compute velocities
   v2=(1.0f/2.0f/dt)*(p2-p2old);

   //Update old position
   p2old = aux;
}


// Exercise 1
// hanging mass point
void AdvanceTimeStep1(float k, float m, float d, float L, float kc, float dt, int meth, float p1, float v1, float& p2, float& v2, float& p2old, bool collision)
{
   switch(meth)
   {
   case 1:
      AdvanceEuler1(k, m, d, L, kc, dt, p1, v1, p2, v2, collision);
      break;

   case 2:
      AdvanceEulerSymplec1(k, m, d, L, kc, dt, p1, v1, p2, v2, collision);
      break;

   case 3:
      AdvanceMidPoint1(k, m, d, L, kc, dt, p1, v1, p2, v2, collision);
      break;

   case 4:
      AdvanceBackEuler1(k, m, d, L, kc, dt, p1, v1, p2, v2, collision);
      break;

   case 5:
      AdvanceVerlet1(k, m, d, L, kc, dt, p1, v1, p2, v2, p2old, collision);
      break;
   }
}


// Exercise 2
// Triangle
void AdvanceTimeStep2(float k, float m, float d, float L, float kA, float A, float dt, int method,
   const Vec2& p1, const Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3, Vec2& p4, Vec2& v4, bool springs, bool area)
{
   switch(method)
   {
   case 2:
      AdvanceEulerSymplec(k, m, d, L, kA, A, dt, p1, v1, p2, v2, p3, v3, p4, v4, springs, area);
      break;

   case 4:
      AdvanceBackEuler(k, m, d, L, kA, A, dt, p1, v1, p2, v2, p3, v3, p4, v4, springs, area);
      break;

   default:
      break;
   }
}


