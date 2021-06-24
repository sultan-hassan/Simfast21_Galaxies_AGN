
/************************************************************************
SimFast21
User defined functions that are used in the simulation, such as fitting functions...
See arXiv:1510.04280 for more details
*************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>   /* header for complex numbers in c */

#include "Input_variables.h"
#include "auxiliary.h"

double Rion(float hmass, double redshift);
double Rrec(float overdensity, double redshift);
double G_H(double redshift);
double XHI(double ratio);
double sfr(float hmass, double z);
//double Qion(double z);
double L_nu (double l, double l912);
double Qion(double logmetal);
double MF_logmet(double redshift);
double ebG(double nu);
double Vcir(double hmass, double redshift);
double Mab(double lnu);

#define CM_PER_MPC 3.0857e24


/* units of s^-1 */
/***************   Rion *************/
/* ionising emissivity */
//#define c1      7.25511524e+39
//#define c1      4.2e39
#define c2      9.367608e+07
//#define c3      4.09039945e-01
//#define c3      0.44
#define c4      2.27044695e+00
#define V_norm  1.e+45

double Rion(float hmass, double redshift){
  double tmp_ion1 =  global_A_gal*hmass*pow((1.0 + redshift ),c4);
  double tmp_ion2 =  pow((hmass/c2),global_C_gal);
  double tmp_ion3 =  exp(- pow((c2/hmass),3.0));
  double tmp_ion4 =  tmp_ion1*tmp_ion2*tmp_ion3;
  double tmp_ion5 =  tmp_ion4*global_fesc_gal/V_norm;
  return tmp_ion5;
}



/***************   Rrec *************/
/* recombination rate */
#define r1 9.84815696
#define r2 1.76105133
#define r3 0.81722878
#define r4 5.06773134

double Rrec(float overdensity, double redshift){
  double rec1 =  1.e-24*r1*pow((1.0 + redshift ),r4);
  double rec2 =  pow((overdensity/r2),r3);
  double rec3 =  pow(rec2/(1.+rec2),4.);
  double DVV = pow(global_L*CM_PER_MPC/global_N_smooth/global_hubble/(1. + redshift),3.0);
  double rec55 = rec1*rec3*DVV/V_norm;
  return rec55;
}


/********* SFR ****************/
/* units: M/yr */

double sfr(float hmass, double z) {
  double logmetal = MF_logmet(z);
  return (Rion(hmass,z)/Qion(logmetal));

}


/***** metallicity Madau & Fragos 2016 (arxiv:1606.07887) ***/
double MF_logmet(double redshift){
  return 0.153 - 0.074*pow(redshift,1.34);
}


/********* Qion  Finlator + 2011 (arxiv:1106.4321)****************/
/* units: yr*s^-1*M^-1 */

double Qion(double logmetal){
  return pow(10,0.639* pow( - logmetal,1./3.) + 52.62 - 0.182);
}



/********* Qion ****************/
/* units: yr*s^-1*M^-1 */
//#define a 2.05707274e-02
//#define b 5.29944413e+01
//double Qion(double z) {
//  return pow(10.0, a*z + b - 45.0);
//}


/********** ratio between recombination rate coefficient (at T=10^4K) and interpolation function of Haardt & Madau (2012) uniform ionising background  ****************/
//double G_H(double redshift){
//  double gh_tmp1 =   1.86956756e-17*pow(redshift,5.)  - 9.05797228e-16*pow(redshift,4.) +  1.56916163e-14*pow(redshift,3.);
//  double gh_tmp2 =   -1.07046180e-13*pow(redshift,2.) + 1.15923648e-13*redshift + 1.04351586e-12;
//  double b_h =  4.19232273531e-13/(gh_tmp1 + gh_tmp2);
//  return b_h;
//e}
double G_H(double redshift){
  double gh_tmp1 =   - 4.37152818e-05*pow(redshift,5.)  + 1.64394635e-03*pow(redshift,4.) - 2.03044756e-02*pow(redshift,3.);
  double gh_tmp2 =   7.18335221e-02*pow(redshift,2.) - 8.85791388e-02*redshift -1.20887500e+01;
  double b_h =  4.19232273531e-13/pow(10.0, gh_tmp1 + gh_tmp2);
  return b_h;
}



/**************************** Popping et al. (2009) formula to compute the residual neutral fraction from ionising background ********/
double XHI(double ratio){
  double XHI_tmp1 = 2.*ratio + 1. - sqrt((2.*ratio + 1.)*(2.*ratio + 1.) - 4.*ratio*ratio  );
  double XHI_tmp3 = XHI_tmp1/(2.*ratio);
  if(ratio == 0.0) return 0.0; 
  else  return XHI_tmp3;
}



/* for xalpha calculation
/* generic function for number photons/baryon/freq.(Hz)  - use: A*nu^-alpha */
/* Assumes A/nu^0.9 SED with 20,000 Lyman photons/baryon between Lya and Ly_limit */
/* We can change A and index 0.9 (A is degenerate with SFR efficiency */
/* Frequency in Hz */
double ebG(double nu) {

  double A_Lya = 1979.878526;
  double alpha_Lya = 0.9;

  return A_Lya*pow(nu,-alpha_Lya);

}



/* AGN source model (Hassan et al. 2018, MNRAS, 473, 227, arXiv:1705.05398)*/

/* convert halo mass to circular velocity*/
double Vcir(double hmass, double redshift){
  double cvir =  0.014*pow(hmass*sqrt(global_omega_m*(1.0+redshift)*(1.0+redshift)*(1.0+redshift)+global_lambda),1.0/3.0);
  return cvir; // units of km/sec
}

/* convert halo circular velocity to black hole mass following local universe observations (Tremaine et al. 2002; Ferrarese 2002)  */
double BHmass(double vc){
  double tmp_mb = global_A_agn*pow(vc/159.4,5.0); // solar mass unit
  return tmp_mb;
}

/* Luminosity in B-band */
double LB(double mb){
  double tmp_lb = 5.7e3*pow(mb, 1.0 + global_C_agn); // units of solar luminosity in B-band
  return tmp_lb;
}

/* Luminosity at Lyman limit (912 A)*/

double L_912(double lb){
  double tmp_lnu = lb*pow(10.0,18.05); // units of erg s^-1 Hz^-1
  return tmp_lnu;
}


double Mab(double lnu){
  return -2.5*log10(lnu) + 51.60; // if lnu in erg s^-1 Hz^-1
}



/* Luminosity at any wavelength  */

double L_nu(double l, double l912){
  return l912*pow( 912/l ,-1.57); // units of erg s^-1 Hz^-1   
}


/* AGN ionising emissivity, integrating L_912 over an SED with a power law index -1.57 */

#define planck_const 6.626196e-27 // units of erg s

double Ragn(double l912){
  double tmp_ragn = global_fesc_agn*l912/(1.57*planck_const*V_norm);  // units of s^-1
  return tmp_ragn;
}


/* Quasar Halo Occupancy Distribution (QHOD) fit to Giallongo et al (2015) observations at z ~ 6  (AGN duty cycle)*/
double QHOD(double hmass_bin){

  double tmp_qhod = pow(hmass_bin/(2.18999687e+12), 9.03924875e-01) + 2.30988763e-02; //fraction of halos with active AGN in bins with mean halo mass of *hmass_bin*
  return tmp_qhod;
}




/* Cosmic Time at redshift z*/
double CosmicTime(double redshift){
  double tol=1.e-5;
  double etalast,t=1.;
  double f,fprime,aextemp;
  float etaold;
  float totMass,Lambda;
  float t0,Ho;
  int it;
  Lambda = global_lambda;
  totMass = global_omega_m;
  Ho = H0*0.7;
  it = 0;
  aextemp = 1./(1.+redshift);
  t = t0;
  do {
    f = sqrt(Lambda/totMass)*pow(aextemp,1.5) +
      sqrt(aextemp*aextemp*aextemp*Lambda/totMass+1) -
      exp(1.5*sqrt(Lambda)*Ho*t);
    fprime =-1.5*sqrt(Lambda)*Ho*exp(1.5*sqrt(Lambda)*Ho*t);
    etalast = t;
    t = t-f/fprime;
    if( (it++) > 10 ) break;
  } while( fabs(t-etalast)/etalast > tol );
  etaold = t;
  return t;
}
