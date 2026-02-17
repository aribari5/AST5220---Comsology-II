#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{

  //=============================================================================
  // TODO: Compute OmegaR, OmegaNu, OmegaLambda, H0, ...
  //=============================================================================

  H0 = 100.0*h;                                                                 // Hubble parameter today in km/s/Mpc
  H0_SI = H0 * Constants.km / Constants.Mpc;                                    // Hubble parameter today in 1/s
  OmegaR = ( pow(Constants.pi,2)*pow(Constants.k_b*TCMB,4)*8*Constants.G ) 
  / ( 15.0*pow(Constants.c,5)*pow(Constants.hbar,3)*3.0*pow(H0_SI,2) );         // Radiation density today

  OmegaNu = Neff*(7.0/8.0)*pow(4.0/11.0,4.0/3.0)*OmegaR;                        // Neutrino density today

  OmegaLambda = 1.0 - OmegaB - OmegaCDM - OmegaR - OmegaNu - OmegaK;            // Dark energy density today
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  x_start = -10.0;
  x_end   = 0.0;
  npts    = 100;

  Vector x_array = Utils::linspace(x_start, x_end, npts);



  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
    //...
    //...

    detadx[0] = 0.0;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================
  // ...
  // ...
  // ...
  // ...

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{
  //=============================================================================
  // TODO: The Hubble parameter as a function of x = exp(a). 
  //=============================================================================

  double H = H0 * sqrt( 
    (OmegaB + OmegaCDM)*exp(-3.0*x) 
    + (OmegaR + OmegaNu)*exp(-4.0*x) 
    + OmegaLambda 
    + OmegaK*exp(-2.0*x) );

  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  //=============================================================================
  // TODO: The conformal Hubble parameter as a function of x = exp(a). Using Hp = a*H = exp(x)*H
  //=============================================================================

  double Hp = H_of_x(x) * exp(x);

  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  //=============================================================================
  // The derivative of the conformal Hubble parameter with respect to x. 
  //=============================================================================

  // Intermediate variables
  double H = H_of_x(x);         
  double Hp = Hp_of_x(x);
  


  double dHpdx = Hp + exp(x)*pow(H0,2.0)/(2*H)
                * ( -3.0*(OmegaB + OmegaCDM)*exp(-3.0*x) 
                    -4.0*(OmegaR + OmegaNu)*exp(-4.0*x) 
                    -2.0*OmegaK*exp(-2.0*x) );

  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  //=============================================================================
  // The double derivative of the conformal Hubble parameter with respect to x. 
  //=============================================================================


  // Intermediate variables
  double H = H_of_x(x);         
  double Hp = Hp_of_x(x);
  double dHpdx = dHpdx_of_x(x);
  double dHdx = pow(H0,2.0)/(2.0*H) *
                ( -3.0*(OmegaB + OmegaCDM)*exp(-3.0*x) 
                    -4.0*(OmegaR + OmegaNu)*exp(-4.0*x) 
                    -2.0*OmegaK*exp(-2.0*x) );







  double ddHpddx = dHpdx + pow(H0,2.0)/2.0 * 
                  (exp(x)*H-exp(x)*dHdx)/pow(H,2.0) *
                  ( -3.0*(OmegaB + OmegaCDM)*exp(-3.0*x) 
                    -4.0*(OmegaR + OmegaNu)*exp(-4.0*x) 
                    -2.0*OmegaK*exp(-2.0*x) ) +
                    exp(x)/H * 
                    ( 9.0*(OmegaB + OmegaCDM)*exp(-3.0*x) 
                    +16.0*(OmegaR + OmegaNu)*exp(-4.0*x) 
                    +4.0*OmegaK*exp(-2.0*x))
  return ddHpddx;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

