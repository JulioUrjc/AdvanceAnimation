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

// Clase de Muelle en 1D
class Muelle1D{
public:
   float p1;
   float p2;
   float l;
   float k;

   Muelle1D(void){}
   Muelle1D(float p1, float p2, float l, float k) : p1(p1), p2(p2), l(l), k(k) {}
   ~Muelle1D(void){}

   float fuerza1(void) const{ return -k*(p1-p2-l); }
   float fuerza2(void) const{ return -k*(p2-p1+l); }
   float dF1dp1(void) const { return -k;}
   float dF1dp2(void) const { return k; }
   float dF2dp1(void) const { return k; }
   float dF2dp2(void) const { return -k;}
};

// Clase muelle 2D
class Muelle2D{
public:
   Vec2 p1;
   Vec2 p2;
   float L0;
   float L;
   float k;
   Vec2 u1;
   Vec2 F;
   Matrix2 dFdp;

   Muelle2D(void){}
   Muelle2D(const Vec2& p1, const Vec2& p2, float L0, float k) : p1(p1), p2(p2), L0(L0), k(k){
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
   ~Muelle2D(void){}

   const Vec2& fuerza1(void) const   {return F;}
   Vec2 fuerza2(void) const          {return (-1.0)*F;}
   const Matrix2& dF1dp1(void) const {return dFdp;}
   Matrix2 dF1dp2(void) const        {return (-1.0)*dFdp;}
   Matrix2 dF2dp1(void) const        {return (-1.0)*dFdp;}
   const Matrix2& dF2dp2(void) const {return dFdp;}
};

class Triangulo{
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

   Triangulo(void){}
   Triangulo(const Vec2& p1, const Vec2& p2, const Vec2& p3, float A0, float k) : p1(p1), p2(p2), p3(p3), A0(A0), k(k){
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
   ~Triangulo(void){}

   Vec2 force1(void) const{ return (k*(A0-A)) * dAdp1T; }
   Vec2 force2(void) const{ return (k*(A0-A)) * dAdp2T; }
   Vec2 force3(void) const{ return (k*(A0-A)) * dAdp3T; }

   Matrix2 dF1dp1(void) const{return ((-k) * dAdp1T) * dAdp1T;}
   Matrix2 dF1dp2(void) const{return (k*(A0-A)) * d2AdpTdnext - (k * dAdp1T) * dAdp2T;}
   Matrix2 dF1dp3(void) const{return (k*(A0-A)) * d2AdpTdprev - (k * dAdp1T) * dAdp3T;}
   Matrix2 dF2dp1(void) const{return (k*(A0-A)) * d2AdpTdprev - (k * dAdp2T) * dAdp1T;}
   Matrix2 dF2dp2(void) const{return ((-k) * dAdp2T) * dAdp2T;}
   Matrix2 dF2dp3(void) const{return (k*(A0-A)) * d2AdpTdnext - (k * dAdp2T) * dAdp3T;}
   Matrix2 dF3dp1(void) const{return (k*(A0-A)) * d2AdpTdnext - (k * dAdp3T) * dAdp1T;}
   Matrix2 dF3dp2(void) const{return (k*(A0-A)) * d2AdpTdprev - (k * dAdp3T) * dAdp2T;}
   Matrix2 dF3dp3(void) const{return ((-k) * dAdp3T) * dAdp3T;}
};

// Euler Explicito Ejercicio 1
void EulerExp(float k, float m, float d, float l, float dt, float p1, float v1, float& p2, float& v2, float kc, bool collision){

   //Inicializamos fuerzas
   float gravedad=-m*g;
   float f2=0.0f;

   //Incluimos gravedad
		f2+=gravedad;
   //Fuerzas del muelle 
		Muelle1D muelle12(p1, p2, l, k);
		f2+=muelle12.fuerza2();
   //Fuerza amortiguamiento
		f2-=d*v2;
   //Calculo de colision
		if (collision && (p2 < 0.0f)) f2-= kc*p2;
   //Integramos la posicion
		p2+=dt*v2;
   //Integramos la velocidad
		v2+=dt*(1.0f/m)*f2;
}

// Euler Simplectico Ejercicio 1
void EulerSimp(float k, float m, float d, float l, float dt, float p1, float v1, float& p2, float& v2, float kc, bool collision){
   
	//Inicializamos
   float gravedad=-m*g;
   float f2=0.0f;

   //Incluimos gravedad
		f2+=gravedad;
   //Fuerzas del muelle
		Muelle1D muelle12(p1, p2, l, k);
		f2+=muelle12.fuerza2();
   //Fuerza amortiguamiento
		f2-=d*v2;
   //Calculo de colision
		if (collision && (p2 < 0.0f)) f2-= kc*p2;
   //Integramos la velocidad
		v2+=dt*(1.0f/m)*f2;
   //Integramos la posicion
		p2+=dt*v2;
}

// Punto Medio Ejercicio 1
void PuntoMedio(float k, float m, float d, float l, float dt, float p1, float v1, float& p2, float& v2, float kc, bool collision){
   
	//Inicializamos
	float gravedad=-m*g;

   //aeleracion: aceleracion en t
		Muelle1D muelle12(p1, p2, l, k);
		float aceleracionActual = (1.0f/m)*(gravedad - d*v2 + muelle12.fuerza2());
   //Calculo de colision
		if (collision && (p2 < 0.0f)) aceleracionActual-= (1.0f/m)*(kc*p2);
   //v2medio: velocidad en t+h/2
		float v2medio = v2 + (dt/2.0f)*aceleracionActual;
   //p2medio: posicion en t+h/2
		float p2medio = p2 + (dt/2.0f)*v2;
   //a2medio: aceleracion en t+h/2 con p&v en t+h/2
		muelle12= Muelle1D(p1, p2medio, l, k);
		float a2medio = (1.0f/m)*(gravedad - d*v2medio + muelle12.fuerza2());
		if (collision && (p2medio < 0.0f)) a2medio-= (1.0f/m)*(kc*p2medio);

   //Calcular posicion
		p2+=dt*v2medio;
   //Calcular velocidad
		v2+=dt*a2medio;
}

// Euler Implicito Ejercicio 1
void EulerImp(float k, float m, float d, float l, float dt,float p1, float v1, float& p2, float& v2, float kc, bool collision){
   
	//Calculamos fuerzas
   Muelle1D muelle12(p1, p2, l, k);
   float gravedad = -m*g;
   float fuerza = gravedad - d*v2 + muelle12.fuerza2();
   //Calculo de colision
	if (collision && (p2 < 0.0f)) fuerza-= kc*p2;

   //Configuramos el sistema
		float m_impl = m+ dt*d - dt*dt*(muelle12.dF2dp2()- ((collision && p2 < 0.0f) ? kc : 0.0f));
		float b = (m+dt*d)*v2 + dt*fuerza;
   //Calculamos velocidad
		v2 = b / m_impl;
   //Integramos la posicion
		p2+=dt*v2;
}

// Verlet Ejercicio 1
void Verlet(float k, float m, float d, float l, float dt, float p1, float v1, float& p2, float& v2, float& p2old, float kc, bool collision){
   
	//Salvamos la posicion
		float oldPos = p2;
   //Inicializamos
		float gravedad=-m*g;
		float f2=0.0f;

   //Incluimos gravedad
		f2+=gravedad;
   //Fuerzas del muelle
		Muelle1D muelle12(p1, p2, l, k);
		f2+=muelle12.fuerza2();
   //Fuerza amortiguamiento
		f2-=d*v2;
   //Calculo de colision
		if (collision && (p2 < 0.0f)) f2-= kc*p2;
   //Integramos las posiciones
		p2=2.0f*p2 - p2old + dt*dt*(1.0f/m)*f2;
   //Calculamos la velocidad
		v2=(1.0f/2.0f/dt)*(p2-p2old);
   //Cargamos la posicion antigua
		p2old = oldPos;
}

// Exercise 1
// hanging mass point
void AdvanceTimeStep1(float k, float m, float d, float L, float kc, float dt, int meth, float p1, float v1, float& p2, float& v2, float& p2old, bool collision)
{
	switch(meth){
		case 1:
			EulerExp(k, m, d, L, dt, p1, v1, p2, v2, kc, collision);
			break;
		case 2:
			EulerSimp(k, m, d, L, dt, p1, v1, p2, v2, kc, collision);
			break;
		case 3:
			PuntoMedio(k, m, d, L, dt, p1, v1, p2, v2, kc, collision);
			break;
		case 4:
			EulerImp(k, m, d, L, dt, p1, v1, p2, v2, kc, collision);
			break;
		case 5:
			Verlet(k, m, d, L, dt, p1, v1, p2, v2, p2old, kc, collision);
			break;
		default:
			break;
   }
}

// Euler Simplectico 2D Ejercicio 2
void EulerSymplec2(float k, float m, float d, float L, float kA, float A0, float dt,
   const Vec2& p1, const Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3, Vec2& p4, Vec2& v4, bool springs, bool area)
{
   //Inicializacion
   Vec2 gravedad(0.0f, -m*g);
   //Vec2 f1(0.0f, 0.0f);
   Vec2 f2(0.0f, 0.0f);
   Vec2 f3(0.0f, 0.0f);
   Vec2 f4(0.0f, 0.0f);

   //Sumamos la gravedad
	   //f1+=gravedad;
	   f2+=gravedad;
	   f3+=gravedad;
	   f4+=gravedad;

   //Calculamos fuerzas de los muelles
   if (springs){
      Muelle2D muelle12(p1, p2, L, k);
      Muelle2D muelle23(p2, p3, L, k);
      Muelle2D muelle34(p3, p4, L, k);
	  Muelle2D muelle41(p4, p1, L, k);
      //f1+=muelle12.fuerza1()+muelle41.fuerza2();
      f2+=muelle23.fuerza1()+muelle12.fuerza2();
      f3+=muelle34.fuerza1()+muelle23.fuerza2();
	  f4+=muelle41.fuerza1()+muelle34.fuerza2();
   }

   //Sumamos el amortiguamiento
	   //f1-=d*v1;
	   f2-=d*v2;
	   f3-=d*v3;
	   f4-=d*v4;

   //Fuerzas para mantener el area
   if (area){
      Triangulo tri1(p1, p2, p3, A0/2, kA);
	  Triangulo tri2(p3, p4, p1, A0/2, kA);
      f2+=tri1.force2();
      f3+=tri1.force3() + tri2.force1();
      f4+=tri2.force2();
   }

   //Integramos
	   v2+=dt*(1.0f/m)*f2;
	   p2+=dt*v2;
	   v3+=dt*(1.0f/m)*f3;
	   p3+=dt*v3;
	   v4+=dt*(1.0f/m)*f4;
	   p4+=dt*v4;
}

void EulerImp2(float k, float m, float d, float l, float kA, float A0, float dt,
   const Vec2& p1, const Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3, Vec2& p4, Vec2& v4, bool springs, bool area){

		//Inicializar Fuerzas
		Vec2 gravedad(0.0f, -m*g);

		//Vec2 f1(0.0f, 0.0f);
		Vec2 f2(0.0f, 0.0f);
		Vec2 f3(0.0f, 0.0f);
		Vec2 f4(0.0f, 0.0f);

		//Calculamos gravedad
		//f1+=mg;
		f2 += gravedad;
		f3 += gravedad;
		f4 += gravedad;

		//Initialiciacion de las Jacobianas
		Matrix2 dF2dp2 = Matrix2::ZERO, dF2dp3 = Matrix2::ZERO, dF2dp4 = Matrix2::ZERO;
		Matrix2 dF3dp2 = Matrix2::ZERO, dF3dp3 = Matrix2::ZERO, dF3dp4 = Matrix2::ZERO;
		Matrix2 dF4dp2 = Matrix2::ZERO, dF4dp3 = Matrix2::ZERO, dF4dp4 = Matrix2::ZERO;

		//Calculo de los Muelles
		if (springs){
			Muelle2D muelle12(p1, p2, l, k);
			Muelle2D muelle23(p2, p3, l, k);
			Muelle2D muelle34(p3, p4, l, k);
			Muelle2D muelle41(p4, p1, l, k);

			//f1 += muelle12.fuerza1() + muelle31.fuerza2();
			f2 += muelle23.fuerza1() + muelle12.fuerza2();
			f3 += muelle34.fuerza1() + muelle23.fuerza2();
			f4 += muelle41.fuerza1() + muelle34.fuerza2();

			//Collect spring Jacobians
			dF2dp2 += muelle12.dF2dp2() + muelle23.dF1dp1();
			dF2dp3 += muelle23.dF1dp2();
			dF2dp4 += Matrix2::ZERO;
			dF3dp2 += muelle23.dF2dp1();
			dF3dp3 += muelle23.dF2dp2() + muelle34.dF1dp1();
			dF3dp4 += muelle34.dF1dp2();
			dF4dp2 += Matrix2::ZERO;
			dF4dp3 += muelle34.dF2dp1();
			dF4dp4 += muelle34.dF2dp2() + muelle41.dF1dp1();
		}

		if (area){
			Triangulo tri1(p1, p2, p3, A0 / 2, kA);
			Triangulo tri2(p3, p4, p1, A0 / 2, kA);
			//f1+=tri1.force1()+tri2.force3();
			f2 += tri1.force2();
			f3 += tri1.force3() + tri2.force1();
			f4 += tri2.force2();

			dF2dp2 += tri1.dF2dp2();
			dF2dp3 += tri1.dF2dp3();
			dF2dp4 += Matrix2::ZERO;
			dF3dp2 += tri1.dF3dp2();
			dF3dp3 += tri1.dF3dp3() + tri2.dF1dp1();
			dF3dp4 += tri2.dF1dp2();
			dF4dp2 += Matrix2::ZERO;
			dF4dp3 += tri2.dF2dp1();
			dF4dp4 += tri2.dF2dp2();
		}

		//Calculamos Amortiguamiento
		//f1-=d*v1;
		f2 -= d*v2;
		f3 -= d*v3;
		f4 -= d*v4;

		//Set up A matrix -> A*v=b
		//Axy = m - dt*dFxdvy - dt^2*dfxdpy

		Matrix2 A11 = (m + dt*d)*Matrix2::IDENTITY - (dt*dt)*dF2dp2;
		Matrix2 A12 = (-dt*dt)*dF2dp3;
		Matrix2 A13 = (-dt*dt)*dF2dp4;

		Matrix2 A21 = (-dt*dt)*dF3dp2;
		Matrix2 A22 = (m + dt*d)*Matrix2::IDENTITY - (dt*dt)*dF3dp3;
		Matrix2 A23 = (-dt*dt)*dF3dp4;

		Matrix2 A31 = (-dt*dt)*dF4dp2;
		Matrix2 A32 = (-dt*dt)*dF4dp3;
		Matrix2 A33 = (m + dt*d)*Matrix2::IDENTITY - (dt*dt)*dF4dp4;

		//Set up b vector
		//Vec2 b1 = (m + dt*d)*v1 + dt*f1;
		Vec2 b2 = (m + dt*d)*v2 + dt*f2;
		Vec2 b3 = (m + dt*d)*v3 + dt*f3;
		Vec2 b4 = (m + dt*d)*v4 + dt*f4;

		// Escribimos los valores en la matriz de forma escalonada
		MatrixMN A(6, 6);
		A(0, 0) = A11.v[0][0]; A(0, 1) = A11.v[0][1]; A(1, 0) = A11.v[1][0]; A(1, 1) = A11.v[1][1];
		A(0, 2) = A12.v[0][0]; A(0, 3) = A12.v[0][1]; A(1, 2) = A12.v[1][0]; A(1, 3) = A12.v[1][1];
		A(0, 4) = A13.v[0][0]; A(0, 5) = A13.v[0][1]; A(1, 4) = A13.v[1][0]; A(1, 5) = A13.v[1][1];

		A(2, 0) = A21.v[0][0]; A(2, 1) = A21.v[0][1]; A(3, 0) = A21.v[1][0]; A(3, 1) = A21.v[1][1];
		A(2, 2) = A22.v[0][0]; A(2, 3) = A22.v[0][1]; A(3, 2) = A22.v[1][0]; A(3, 3) = A22.v[1][1];
		A(2, 4) = A23.v[0][0]; A(2, 5) = A23.v[0][1]; A(3, 4) = A23.v[1][0]; A(3, 5) = A23.v[1][1];

		A(4, 0) = A31.v[0][0]; A(4, 1) = A31.v[0][1]; A(5, 0) = A31.v[1][0]; A(5, 1) = A31.v[1][1];
		A(4, 2) = A32.v[0][0]; A(4, 3) = A32.v[0][1]; A(5, 2) = A32.v[1][0]; A(5, 3) = A32.v[1][1];
		A(4, 4) = A33.v[0][0]; A(4, 5) = A33.v[0][1]; A(5, 4) = A33.v[1][0]; A(5, 5) = A33.v[1][1];

		Vector b(6);
		b(0) = b2.x; b(1) = b2.y;
		b(2) = b3.x; b(3) = b3.y;
		b(4) = b4.x; b(5) = b4.y;
		Vector x(6);
		GaussElimination(A, b, x);
		v2.x = x(0); v2.y = x(1);
		v3.x = x(2); v3.y = x(3);
		v4.x = x(4); v4.y = x(5);

		// Integramos la posicion
		p2 += dt*v2;
		p3 += dt*v3;
		p4 += dt*v4;
}

// Exercise 2
// square
void AdvanceTimeStep2(float k, float m, float d, float L, float kA, float A, float dt, int method,
   const Vec2& p1, const Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3, Vec2& p4, Vec2& v4, bool springs, bool area){
	switch(method){
		case 2:
			EulerSymplec2(k, m, d, L, kA, A, dt, p1, v1, p2, v2, p3, v3, p4, v4, springs, area);
			break;

		case 4:
			EulerImp2(k, m, d, L, kA, A, dt, p1, v1, p2, v2, p3, v3, p4, v4, springs, area);
			break;
		default:
			break;
   }
}


