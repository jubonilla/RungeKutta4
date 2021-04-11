#include <fstream>
#include <cmath>
#include <vector>

typedef std::vector<double> Vector;
const double G=9.8; //Aceleración gravitacional local
const double L=1.0; //Longitud de la cuerda
const double H=0.1; //Tamaño del intervalo
const double T0=0.0; //Tiempo inicial
const double TF=30.0; //Tiempo final
const double Theta0=1.0; //Theta inicial
const double dTheta0=0.0; //Velocidad angular inicial
const double  Omega0=G/L; //Frecuencia angular

void k_1_m ( Vector & k_i, const Vector & theta);
void k_2_m ( Vector & k_i, const Vector & theta, const Vector & k_ip);
void k_3_m ( Vector & k_i, const Vector & theta, const Vector & k_ip);
void k_4_m ( Vector & k_i, const Vector & theta, const Vector & k_ip);
void theta_m (double t0,double tf, Vector & k_1,Vector & k_2,Vector & K_3,Vector & k_4, Vector theta);

int main (void)
{
  Vector Theta (4);
  Theta = {Theta0, dTheta0, 0.0, 0.0};

  Vector K_1={0.0, 0.0}; 
  Vector K_2={0.0,0.0};
  Vector K_3={0.0,0.0};
  Vector K_4={0.0,0.0};
  //Vectores de dos elementos que corresponden a las dos ecuaciones diferenciales de primer orden sobre las que se aplica el método RK-4

  theta_m (T0, TF, K_1, K_2, K_3, K_4, Theta);

  return 0;
}

void k_1_m ( Vector & k_i, const Vector & theta){
  double a=theta[0]; 
  double y=theta[1];

 k_i[0]=y;
 k_i[1]=-Omega0*std::sin(a);
}

void k_2_m ( Vector & k_i, const Vector & theta, const Vector & k_ip){
  double a=theta[0]+((H/2)*k_ip[0]);
  double y=theta[1]+((H/2)*k_ip[1]);

 k_i[0]=y;
 k_i[1]=-Omega0*std::sin(a);
}

void k_3_m ( Vector & k_i, const Vector & theta, const Vector & k_ip){
  double a=theta[0]+((H/2)*k_ip[0]);
  double y=theta[1]+((H/2)*k_ip[1]);

 k_i[0]=y;
 k_i[1]=-Omega0*std::sin(a);
}

void k_4_m ( Vector & k_i, const Vector & theta, const Vector & k_ip){
  double a=theta[0]+(H*k_ip[0]);
  double y=theta[1]+(H*k_ip[1]);

   k_i[0]=y;
   k_i[1]=-Omega0*std::sin(a);
}

void theta_m (double t0,double tf,Vector & k_1,Vector & k_2,Vector & k_3,Vector & k_4, Vector theta)
{
  std::ofstream fout; 
  fout.open ("Datos.txt"); //Abre el archivo para guardar los datos
  int nsteps=(tf-t0)/H;
  for (int ii=0;ii<=nsteps;ii++){

    double t=t0+ii*H;

     k_1_m(k_1,theta);
     k_2_m(k_2,theta,k_1);
     k_3_m(k_3,theta,k_2);
     k_4_m(k_4,theta,k_3);

  theta[2]=theta[0]+((H/6)*(k_1[0]+(2*k_2[0])+(2*k_3[0])+k_4[0]));
  theta[3]=theta[1]+((H/6)*(k_1[1]+(2*k_2[1])+(2*k_3[1])+k_4[1]));
  fout<<t<<" "<<theta[0]<<" "<<theta[1]<<std::endl;
    theta[0]=theta[2];
    theta[1]=theta[3];
    }
 fout.close();
}