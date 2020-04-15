/*

Program written by:

Winthrop Brown
MIT Lincoln Laboratory
781-981-6716
wbrown@ll.mit.edu
or
winthrop@alum.mit.edu

Last Modified: October, 2011.

This program uses the differential relativistic Thomson cross-section to calculate the scattered
x-ray spectrum and/or intensity profile of x-rays produced by an electron bunch represented by loaded
macro-particles with a Guassian profile laser pulse. Effect of laser bandwidth as well as perpendicular
k-vector distributions are included. It is assumed that the radiation is incoherent, hence the contributions
from each macro-paritcle are added linearly. Key words in input file "Compton.ini".
Initial particle coordinates are read from the file with format from "tape2" produced from PARMELA.

   Dimensionless Units used in simulation:

   t' =  t/t_unit

   x' = x/l_unit

   v' =  v*t_unit/l_unit

   p' =  gamma*v/c




  Coordinate Variables

   X[6] = [x',y',z',px',py',pz']


    Key words (for "compton.ini"):
    For up to date explanation, see example file of compton.ini.

  WAVELENGTH: Wavelength of laser (in meters).
  MAXFIELD:   Maximum normalized verctor potential (a0) at laser focus (z=0).
  TBASE:      Time base unit (in seconds)
  LBASE:      Length base unit (in meters)
  DT:         Time step in base time unit;
  OFFSET:     3 Arguments specifying offset in x,y,z (in meters) of the electron focus from the laser focus
  OFFSET2:    2 arguments specifying the angular offset in dx/dz and dy/dz (in radians).
  BEAM_WAIST: 1 argument specifying the minimum 1/e^2 intensity radius for the laser spot.
  DURATION:   Scale factor for electron bunch duration
  PULSE_WIDTH: Laser Pulse Duration in seconds.
  EX:         Signifies no more key words in file


  */

//#include "stdafx.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include "Include\Header\random2.h"
#include "Include\Header\matrix2.h"

using namespace std;

#define prec double
#define INFILE "Compton.ini" //simulation initialization parameters
#define INTENFILE "profile.txt"
#define INTENFILE2 "data_dump_profile.txt"

#define OUT_TIME_UNIT 1.e-12 //Time unit used in writing particle coordinates to file (ps)
#define OUT_TIME_UNIT_NAME "ps"
#define OUT_LENGTH_UNIT 1.e-3
#define OUT_LENGTH_UNIT_NAME "mm"

#define c 2.9979e8//speed of light;
#define LAMBDA_C 2.4263e-12    //compton wavelength (m)
#define R_E 2.8179e-15        //Classical Electron Radius (m)
#define REST 0.511 //rest mass of an electron in MeV
#define RATIO 1.7588e11 //electron charge to mass ratio (SI units)
#define CONV  0.01 //conversion from cm to m
#define PLANCK 6.62607e-34
#define E_CHARGE 1.6022e-19
#define NUMSKIP 1000  //number of particles to skip in outputing to screen
#define NUMSKIP_CHIRP 1000

#define NSTAT 18            //number of statistics recorded for each particle
#define WS 10

#define LEFT  ios::left, ios::adjustfield
#define PR setprecision(8)
#define PR2 setprecision(5)
#define PR3 setprecision(3)
#define PR4 setprecision(4)
#define W2 14
#define W3 40
#define INCR 1.e-6
#define ZERO 1.e-6
#define INFINITY 1.e8

#define MAX_PARTICLES 250000
#define MAX_XRAY_BIN 500
#define MAX_PHI_NUMBER 100
#define MAX_THETA_NUMBER 400
#define MAX_TIME_BIN 400

#define MAX_DATA_DUMP_FREQ_STEPS 500
#define DATA_DUMP_FILE_DEFAULT "data_dump.txt"

#define RANDOMIZE_START_SCALE 1.0


#define Templ "parm_in"

#define Gen_Part_File_Out "gen_tape.txt"

#define MIN_EXPONENT -150.0

char clear[200];
char TAB;

char CORFILE[50];
char OUTFILE[50];

/*********NEWER STUFF***********/

//The following three variables are utilized if the electron beam coordinates are to be generated within the program
//********************************************************
int GENERATE_PARTICLES = 0; //if 1, electron beam macro_particles are generated rather than read from a file
prec TWISS_6D[9];
/* emit(mm_mrad), alpha, beta(mm/mrad) emit(mm_mrad), alpha, beta(mm/mrad) emit(deg_MeV), alpha, beta(deg/MeV) */
prec AVG_6D[6]; // [mm mrad mm mrad deg MeV]
int NUM_GEN;  //number of particles to generate
//********************************************************

int ANGLE_SPECTRUM = 0;   //1 calculates spectrum for a range of theta angles at a given phi
                            //2 calculates spectrum for a range of phi angels at a given theta
prec ANGLE_SPECTRUM_START;
prec ANGLE_SPECTRUM_STOP;
prec ANGLE_SPECTRUM_FIXED;  //phi or theta fixed angle
int ANGLE_SPECTRUM_NUM = 0;




/********NEW STUFF************/
char WINDOW_FILE[50]; //file containing tabulated transmission data for window (two columns, keV, transmission)
int WINDOW_FILE_ON = 0;
char DETECTOR_FILE[50];
int DETECTOR_FILE_ON = 0;
int WINDOW_FILE_N = 0;  //number of tabulated points in WINDOW_FILE
int DETECTOR_FILE_N = 0; //number of tablulated points in DETECTOR_FILE
long Rand_Start_Init_Value = 0 ; //Initialization variable for starting randomization of integration start times
long Rand_Start_Init_Value_0 = 0;
int RANDOMIZE_START = 0; //Specifies the randomization of start times for macro particle integration (randomizes within one time step)
int RANDOMIZE_START_SPECIFIED = 0;

prec WINDOW_FILE_X[MAX_XRAY_BIN];
prec WINDOW_FILE_Y[MAX_XRAY_BIN];
prec DETECTOR_FILE_X[MAX_XRAY_BIN];
prec DETECTOR_FILE_Y[MAX_XRAY_BIN];


int NUMBER;     //number of particles in simulation
int PLANE_WAVE=1; //if 1, the calculation ignores the perpindicular distribution of k_perp
int SPECTRUM = 0; //if 1, then the spectrum calculation is performed
int INTENSITY = 0; //if 1, then intensity calculaion is performed
int FILTER_ON = 0;  //if 1, then the output spectrum is Filtered.
int WINDOW_ON = 0;  //if 1, then the output spectrum is modified by the presence of a window.
int FAST_INTENSITY =0 ;//if 1, then the intensity profile is calculated with analytic integration over frequencies
int CART_INTENSITY = 0; //if 1, then the inensity profile is calculated in cartesian coordinates


prec off[3];    //offset of particle coordinates in input file (distribution centers around reference particle)
prec off2[3];   //angular offset of dx/dz and dy/dz

prec EMAX;           //Maximum electric field at the focus;
prec PI;    //3.14
prec E_SCALE;        //Scale to convert normalized scattered frequency to keV.

prec TBASE = 1.0e-12;       //default time base is 1.0 picoseconds.
prec LBASE = 1.0e-6;        //default length base is 1.0 micro-meters
prec WAVELENGTH = 0.800e-6;  //default wavelength is 0.8 micrometers;
prec K_VECTOR;              //K vector
prec DT_BASE;              //c*TBASE/LBASE
prec ZR;                  //Raleigh Range of Laser Pulse
prec ZRx;
prec ZRy;
prec w0 = 10.e-6;         // 1/e^2 intensity radius of laser minimum spot size (default = 10 micrometers)
prec PULSE_WIDTH = 1.0e-12; //default 1/e^2 laser intensity pusle width is 1 ps.
prec BANDWIDTH;
prec BANDWIDTH_FACTOR=1.0;
prec AVG_E;
prec AVG_GAMMA;

prec M=1.;               //M^2 of laser focus
prec w0x;                 //1/e^2 waist in x and y
prec w0y;
prec THETA_0 = 0.;       //incidence angle of laser in x-z plane
prec sigma_theta_L_x = 0.;    //rms spread in incident laser photon angles (detemined from M^2 and Laser spot size).
prec sigma_theta_L_y = 0.;
prec A0 = 0.1;           // a0^2: peak normalized intensity at the laser focus (Not time averaged. This is the sqaure of the peak instantaneous field value)
prec PHI_P = 0.0;        //polarization rotation angle (0 -> polarized in x-z plane
prec OBSERVE_THETA=0;    //center observation angles in lab frame.
prec OBSERVE_PHI=0;
prec OBS_DISTANCE = 2.0; //Specifies position of detector plane. Has significance when the
                                //spot size of the interaction is considered.

prec FILTER_WIDTH=0.; //Width of frequency filter (in units of gamma^2)
prec FILTER_CENTER=0.; //Center of frequency filter (in units of gamma^2)
prec FILTER_AMP =1.;   //Peak transimission of filter

prec FULL_PULSE_LENGTH;  //Length incorporating all electrons in the bunch
prec MINIMUM_ELECTRON = INFINITY;
prec MAXIMUM_ELECTRON = -INFINITY;

prec Ws_Max=1000.;
prec Ws_Min=1.0;
int N_Ws = 100;

prec INTEGRAL_RANGE_SIGMA_TIME = 3.0;  //number of laser intensity time duration sigmas to integrate out in time integral
prec INTEGRAL_RANGE_SIGMA_THETA = 3.0; //number of k_perp distribution sigmans to integrate out to in k_perp integration
prec INTEGRAL_RANGE_SIGMA_W = 3.0;     //number of omega 1/e^2 widths to integrate over for non-analytic intensity profile calculations
int THETA_STEPS = 10;

prec N_photon; //total number of photons in incident laser pulse
prec N_electron = 1.0; //number of electrons per macro-particle
prec Laser_Energy = 0.0; //Laser Energy in mJ
int Laser_Energy_Specified = 0;

int T_STEPS = 31;           //Number of time steps to be taken.

int TIME_SCALE_SPECIFIED = 0;
int LENGTH_SCALE_SPECIFIED = 0;
int DT_SPECIFIED = 0;

prec N_X[MAX_TIME_BIN][MAX_XRAY_BIN]; //Stores the number of scattered x-rays at xray energy E_X
prec E_X[MAX_XRAY_BIN];
prec T_X[MAX_TIME_BIN];

int TIME_BIN_NUM = 11;

prec U_X[MAX_PHI_NUMBER][MAX_THETA_NUMBER]; //Stores the total scattered energy intensity at a given phi and theta
prec PHI[MAX_PHI_NUMBER];
prec THETA[MAX_THETA_NUMBER];

int PHI_NUMBER;             //Defines the number of observation points in phi and theta and the range of theta for
int THETA_NUMBER;           //determining the scattered intensity profile.
prec THETA_RANGE;
prec PHI_RANGE_MIN = 0.0;
prec PHI_RANGE_MAX = 1.5707963;


int THETA_NUMBER_X;        //Defines the number and range of observation points for the intensity profile
prec THETA_RANGE_X;        //calculataed in cartesian profile calculation
int THETA_NUMBER_Y;
prec THETA_RANGE_Y;

prec SOLID_THETA_MAX=0.;   //Specifies parameters for solid angle integration of spectrum
prec SOLID_THETA_MIN=0.;
prec SOLID_DTHETA=0.;
prec SOLID_DPHI=0.;
int SOLID_NTHETA;
int SOLID_NPHI;
int SOLID_ANGLE = 0;


int BRAGG=0;         //Specfies paramters for a Bragg reflection
prec BRAGG_SPACING;
prec BRAGG_ANGLE;    //Angle in rad
prec BRAGG_WIDTH = 1.0e-3; //1/e^2 attenuation fractional width

int DATA_DUMP = 0;        //Specifies that the file storing spectral information for each transverse position
char DATA_DUMP_FILE[100]; //for the cartesian intensity profile calculations is to be created.

int TOTAL_PHOTON=0;  //if 1, the program will calculated the total photon flux.
int CHIRP_LASER = 0;  //1 specifies the use of a chirped laser pulse.
int NONLINEAR_ON = 0; //specifies inclusion of non-linear pondermotive effects
prec NONLINEAR_THRESHOLD;
prec CHIRP_WIDTH = 1.0; //fractional increase of pulse width from bandwith limit, if negative => negative chirp

prec cc;               //Speed of light in base units: c*

prec FREQ = 2.7777778e9;   //Frequency used for in "tape2" input file for conversion of phase to time. 
						   //This default frequency corresponds to the case where 1 degree = 1 ps.

/***** Scales Electron Beam input to either change bunch lengh *****/
prec DURATION = 1.;
prec ENERGY_SCALE = 1.;
prec DIVERGENCE_SCALE = 1.0;
prec DE_MULT = 1.0;
prec SPOT_SIZE_SCALE = 1.0;

prec X[MAX_PARTICLES][6];   //electron particle coordinates [x,y,z,px,py,pz]
prec NG[MAX_PARTICLES];     //Time integrated intensity
prec MT[MAX_PARTICLES];     //average interaction time for each particle
int FIRST_TIME = 1;

prec R_E2;  //classical electron beam radius squared.
prec cR_E2_N;          //c*r_e^2 in base units

prec SPOT_SIZE_AVG_X = 0.;
prec SPOT_SIZE_NORM = 0.;
prec SPOT_SIZE_2_AVG_X = 0.;

prec SPOT_SIZE_AVG_Y = 0.;
prec SPOT_SIZE_2_AVG_Y = 0.;

prec SPOT_SIZE_X;
prec SPOT_SIZE_Y;

/************************* PROTOTYPES ****************************/
void rk4 (prec[], prec&, prec, prec, prec);
void DERIVATIVE (prec[], prec&, prec);
prec GAMMA (prec[]); //calculates gamma given a six dimensional coordinate vector [x,y,z,px,py,pz]
prec BETA (prec []);
prec CROSS (prec[], prec[], int); //performs a cross product of two vetors returning component int
void INITCOR (); //initializes particle coordinates
int INPUT();                           //reads in simulations initialization
void COPYARRAY(prec[],prec[],int);
void COPYARRAY_MAX(prec[],prec[],int);
void COPYARRAY_MIN(prec[],prec[],int);
int LaserFrame(prec[],prec[],prec);
prec Photon_Density(prec[], prec);
prec k_perp (prec, prec);
prec XP (prec x[]);
prec YP (prec x[]);
prec alpha_x (prec, prec[], prec, prec );
prec alpha_y (prec, prec[], prec, prec );
prec alpha_z (prec, prec[], prec, prec );
prec SINphiSINtheta (prec[], prec, prec);
prec COSphiSINtheta (prec[], prec, prec);
prec COStheta (prec[], prec, prec);
prec CrossSection(prec[], prec, prec, prec, prec);
prec Current_Position(prec[], int, prec);
prec H(prec[], prec, prec, prec, prec);
prec H_Time_Dependent(prec[], prec, prec, prec, prec, prec, prec);
prec OMEGA_FACTOR (prec[], prec, prec , prec, prec , prec,prec);
prec Integral_ng(prec[],prec, prec);
prec dN_dw_domega (prec[], prec , prec , prec , prec , prec, prec, prec, prec, prec );
prec fk_perp (prec, prec );
prec Integral_dN_y (prec [], prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, prec);
void rk4_dNy (prec[], prec& , prec, prec, prec , prec, prec , prec , prec, prec, prec, prec, prec);
void rk4_dNxy (prec[], prec& , prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, prec);
prec Integral_dN_x_y (prec[], prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, prec);
int Write_Spectrum (ofstream&);
int Get_Spectrum (prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, int, prec);
int Get_Fast_Profile_Plane (int, int, prec, int);
int Get_Profile (int, int, prec, prec, prec, prec, prec, prec, prec, prec, prec, int);
prec Get_Intensity (prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, int);
prec Get_Chirp_Intensity (prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, int);
int Write_Profile (ofstream&, int, int);
int FIND_TIME_BIN(prec, prec[],prec, prec);
prec FIND_TIME_FROM_BIN(int);
prec FILTER(prec, prec, prec);
prec FREQ_ATTEN_FACTOR (prec[],prec, prec, prec);
int Solid_Spectrum (prec, prec, int, int, prec, prec, prec, prec, prec, prec, prec, prec, int);
int Find_W_Range(prec[], prec, prec, prec, prec &, prec &);
prec WINDOW_FILTER (prec);
prec VEL_FACTOR(prec[], prec, prec);
int Get_Start_Stop_Time (prec[], prec& , prec &);
int Adjust_Start_Stop_Time (prec[], prec&, prec& );
int Get_Total_Photon (int,prec&);
int Alter_Observation(prec[], prec &, prec &, prec);
int Get_Chirp_Spectrum (prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, int, prec);
void DERIVATIVE_CHIRP (prec[], prec&,prec, prec, prec, prec, prec, prec, prec);
void rk4_chirp (prec[], prec&, prec, prec, prec, prec, prec, prec, prec , prec, prec);
prec OMEGA_FACTOR_CHRIP (prec[], prec, prec, prec, prec, prec,prec, prec, prec);
prec Integral_ngnw (prec[],prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, prec);
int Get_Cart_Profile (int, prec, int, prec, prec, prec, prec, prec, prec, prec, prec, prec, int);
int Write_Cart_Profile (ofstream&, int, int);
int READ_DATA (char[], prec[], prec[], int &, int);
int Get_Angle_Spectrum_All (prec, prec, prec, int, prec, prec, prec, prec, prec, prec, prec, prec, int, prec);
int Get_Angle_Spectrum_One (int, prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, int, prec);
int Write_Angle_Spectrum (prec, prec, prec, int);
int Find_W_Range_Chirp(prec[], prec, prec, prec, prec, prec, prec, prec, int, prec[], prec[]);
int Find_W_Range_Chirp_Kp(prec[], prec, prec, prec, prec, prec, prec, prec, prec, prec, int, prec[], prec[]);
prec BRAGG_REFLECT(prec ws, int xy, prec b_space, prec b_angle, prec b_width,prec, prec);
int WRITE_DATA_DUMP (ofstream &, prec[], prec, prec, int, prec, prec);
prec Get_Data_Dump_Intensity (prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, int, prec[]);
prec Get_Chirp_Data_Dump_Intensity (prec, prec, prec, prec, prec, prec, prec, prec, prec, prec, int, prec[]);
int Get_Time_Integraion(prec[], prec&, prec&, prec&, prec&, int);
int Get_Fast_Cart_Profile (int, prec, int, prec, prec, prec, prec, prec, prec, prec, int);
int Get_Time_Integraion_2(prec[], prec&, prec&, int);
void GENCOR(int,prec[], prec[]);
int Init_SpotSize_Vars ();
prec INTERPOLATE (prec, prec[], prec[], int);
int Randomize_Start_Time(prec&, prec&);
prec EXP(prec);

/******************************************************************************************/

int main()
{
    TAB = char(9);
    int i,j;
                             //   y^2,gamma*y'^2,gamma*yy',phi,phi^2,E,E^2,x,xp,y,yp]

    int ref = 1;             //denotes use of the reference particle
    PI = 2.*acos(0.0);
    R_E2 = R_E*R_E;

    PHI_RANGE_MAX = PI;

    int numpar;

    prec dW;                 //bin size of calculated spectrum
    prec start_theta;       //start angle for integration over initial laser states in y
    prec stop_theta;         //stop angle for integration over initial laser states in y
    prec d_theta;            //step size of integration in initial laser staes in y
    prec start_theta_x;
    prec stop_theta_x;
    prec d_theta_x;


    ofstream outfile;
//  ofstream fieldfile;
    ofstream intenfile;



	cout<<endl<<"Initializing variables . . ."<<endl;
    strcpy (CORFILE,"smalltape");
    strcpy (OUTFILE,"Compton.out");

    for(i=0; i<3; i++)
    {
        j=i+3;
        off[i] = 0.;
        off2[i] = 0.;
    }

    for (i=0; i<6; i++)
        for (j=0; j<MAX_PARTICLES; j++)
        {
            X[i][j]=0.;
        }//for

    for (i=0; i<MAX_PARTICLES; i++)
    {
        NG[i] = 0.;
    }

    for(j=0;j<MAX_TIME_BIN;j++)
    {
        T_X[j]=0.;
        for (i=0;i<MAX_XRAY_BIN;i++)
        {
            N_X[j][i]=0.;
            E_X[i]=0.;
        }
    }

    for (i=0;i<MAX_PHI_NUMBER;i++)
        for (j=0; j<MAX_THETA_NUMBER;j++)
        {
            U_X[i][j]=0.;
        }

    strcpy(DATA_DUMP_FILE,DATA_DUMP_FILE_DEFAULT);


    if(INPUT()) return(0);

//  cout<<"Time Base = "<<PR3<<TBASE/OUT_TIME_UNIT<<" "<<OUT_TIME_UNIT_NAME<<endl;
//  cout<<"Time Step = "<<PR3<<DT/OUT_TIME_UNIT*TBASE<<" "<<OUT_TIME_UNIT_NAME<<endl;

//***** Assume symmetric laser focus****//
    w0x=w0;
    w0y=w0;
//******************************//

    if (Laser_Energy_Specified)
    {
        N_photon = Laser_Energy/(PLANCK*c/WAVELENGTH*1.0e3);
        A0 = N_photon/((PI*PI*sqrt(PI)/sqrt(32.0))*(c*PULSE_WIDTH*w0x*w0y)/(LAMBDA_C*R_E*WAVELENGTH));
        cout<<"Normalized Intensity (A0) = "<<A0<<endl;
    }

    EMAX=2*PI*sqrt(A0)*TBASE*c/WAVELENGTH;

    E_SCALE = PLANCK*c/(WAVELENGTH*1000.*E_CHARGE);  //base frequency photon energy in keV

  /******Calculated key Simulation Parameters and convert to normalized units *******************/
    //N_photon = A0*(PI*PI*sqrt(PI)/sqrt(32.0))*(c*PULSE_WIDTH*w0x*w0y)/(LAMBDA_C*R_E*WAVELENGTH);
    N_photon = A0*(PI*PI*sqrt(PI)/sqrt(32.0))*(c*PULSE_WIDTH*w0x*w0y)/(LAMBDA_C*R_E*WAVELENGTH);
    Laser_Energy = PLANCK*c/WAVELENGTH*N_photon*1.0e3;//Laser Energy in mJ (for display only)
    BANDWIDTH = BANDWIDTH_FACTOR*WAVELENGTH/(PI*PULSE_WIDTH*c); //changed 11/06/03 from 2.0*WAVELENGTH/(PI*PULSE_WIDTH*c);
    WAVELENGTH = WAVELENGTH/LBASE;  //central wavelength of incident laser
    K_VECTOR= 2.0*PI/WAVELENGTH;    //k_vector magnitude
    DT_BASE = c*TBASE/LBASE;        //temporal integration time step
    cc = c*TBASE/LBASE;             //speed of light in normalized units
    w0 = w0/LBASE;                  //1/e^2 intensity waist of laser focus
    OBS_DISTANCE = OBS_DISTANCE/LBASE;  //Distance of Observation plane from interaction point
    ZR=PI*w0*w0/(M*WAVELENGTH);         //Rayleigh range of laser focus
    w0x=w0;                             //x and y waist are currently the same
    w0y=w0;
    ZRx=ZR;
    ZRy=ZR;
    PULSE_WIDTH = PULSE_WIDTH/TBASE*cc;          //Laser pulse 1/e^2 width in normalized length units
    sigma_theta_L_x = M*WAVELENGTH/(2.0*w0x*PI); //Width of perpendicular k spectrum
    sigma_theta_L_y = M*WAVELENGTH/(2.0*w0y*PI);
    cR_E2_N = R_E2/(LBASE*LBASE)*cc;

  //   cout<<"w0x = "<<w0x<<endl;
//   cout<<"ZR = "<<ZR<<endl;
//   cout<<"PULSE_WIDTH = "<<PULSE_WIDTH<<endl;
/***********************************************************************************************/

/****If NONLINEAR_ON = 1, then check to see that A0 is greater than or equal to the NONLINEAR_THRESHOLD. If not,
*****then set NONLINEAR_ON to 0. Also set CHIRP_LASER to 0 provided CHIRP_WIDTH does not indicate a chriped
*****laser pulse.*/

	if (NONLINEAR_ON & (A0 < NONLINEAR_THRESHOLD) )
	{
		NONLINEAR_ON = 0;
		if (fabs(CHIRP_WIDTH) <= (1.0+ZERO))
			CHIRP_LASER = 0;
	}

    cout<<"Number of Photons in Laser Pulse: "<<PR2<<N_photon<<endl;
    cout<<"Laser Energy = "<<PR2<<Laser_Energy<<" mJ"<<endl;

    for(i=0;i<3;i++)
        off[i]=off[i]/LBASE;

    if(GENERATE_PARTICLES)
    {
        cout<<"Generating macro_particles ..."<<endl;
        GENCOR(NUM_GEN,TWISS_6D,AVG_6D);
        strcpy(CORFILE,Gen_Part_File_Out);
        cout<<"Done"<<endl;
    }

    INITCOR();
    FILTER_CENTER=FILTER_CENTER*AVG_GAMMA*AVG_GAMMA;
    FILTER_WIDTH=FILTER_WIDTH*AVG_GAMMA*AVG_GAMMA;
    BRAGG_SPACING=BRAGG_SPACING/LBASE;


    cout<<"Raleigh Range = "<<PR3<<ZR*LBASE/OUT_LENGTH_UNIT<<" "<<OUT_LENGTH_UNIT_NAME<<endl;

    AVG_E = AVG_E/REST+1.;
    cout<<"Averge gamma = "<<AVG_E<<endl;

    numpar=NUMBER-1;

    N_electron = N_electron/(E_CHARGE*1.0e9*numpar);
    cout<<endl<<"N_e = "<<N_electron<<endl<<endl;

    cout<<"Number of particles in simulation: "<<numpar<<endl;
    cout<<"Integrating particles . . ."<<endl<<endl;

    /*******Define limits and step size of k_perp integration******************/
    start_theta = -INTEGRAL_RANGE_SIGMA_THETA*sigma_theta_L_y;
    stop_theta = INTEGRAL_RANGE_SIGMA_THETA*sigma_theta_L_y;
    d_theta = (stop_theta-start_theta)/(THETA_STEPS-1);

    start_theta_x = -INTEGRAL_RANGE_SIGMA_THETA*sigma_theta_L_y+THETA_0;
    stop_theta_x = INTEGRAL_RANGE_SIGMA_THETA*sigma_theta_L_y+THETA_0;
    d_theta_x = (stop_theta_x-start_theta_x)/(THETA_STEPS-1);
    /*************************************************************************/

    /*******Define limits and bin size of calculated x-ray spectrum******************/
    Ws_Max=Ws_Max*4.0*AVG_GAMMA*AVG_GAMMA;
    Ws_Min=Ws_Min*4.0*AVG_GAMMA*AVG_GAMMA;
    dW=(Ws_Max-Ws_Min)/prec(N_Ws-1);
    /********************************************************************************/


    if (WINDOW_FILE_ON)
        READ_DATA (WINDOW_FILE, WINDOW_FILE_X, WINDOW_FILE_Y, WINDOW_FILE_N, MAX_XRAY_BIN);
    if (DETECTOR_FILE_ON)
        READ_DATA (DETECTOR_FILE, DETECTOR_FILE_X, DETECTOR_FILE_Y, DETECTOR_FILE_N, MAX_XRAY_BIN);

    if(TOTAL_PHOTON)
    {
        ofstream dose_file;
        dose_file.open("Total_Dose.txt");
        prec total_dose = 0.;
        Get_Total_Photon(numpar,total_dose);
        SPOT_SIZE_X = SPOT_SIZE_2_AVG_X/SPOT_SIZE_NORM - SPOT_SIZE_AVG_X*SPOT_SIZE_AVG_X/(SPOT_SIZE_NORM*SPOT_SIZE_NORM);
        SPOT_SIZE_Y = SPOT_SIZE_2_AVG_Y/SPOT_SIZE_NORM - SPOT_SIZE_AVG_Y*SPOT_SIZE_AVG_Y/(SPOT_SIZE_NORM*SPOT_SIZE_NORM);
        if(SPOT_SIZE_X > 0)
        {
            SPOT_SIZE_X = (LBASE/1.0e-6)*sqrt(SPOT_SIZE_X);
            cout<<"Rms x-ray source size in X: "<<SPOT_SIZE_X<<" micro_meters."<<endl;
            dose_file<<"Rms_x-ray_source_size_X: "<<SPOT_SIZE_X<<" micro_meters"<<endl;
        }
        else
        {
            cout<<"Error computing x-ray source size in X."<<endl;
        }

        if(SPOT_SIZE_Y > 0)
        {
            SPOT_SIZE_Y = (LBASE/1.0e-6)*sqrt(SPOT_SIZE_Y);
            cout<<"Rms x-ray source size in Y: "<<SPOT_SIZE_Y<<" micro_meters."<<endl;
            dose_file<<"Rms_x-ray_source_size_Y: "<<SPOT_SIZE_Y<<" micro_meters"<<endl;
        }
        else
        {
            cout<<"Error computing x-ray source size in Y."<<endl;
        }
        dose_file<<"Total_Dose: "<<PR2<<total_dose<<" photons"<<endl;
        dose_file.close();
		
        Rand_Start_Init_Value = 0;
        //return(0);
    }//if

    outfile.open(OUTFILE);
    outfile.setf(LEFT);    

    if(SPECTRUM && !SOLID_ANGLE && !INTENSITY && !CART_INTENSITY && !ANGLE_SPECTRUM)
    {
        if(!CHIRP_LASER)
           Get_Spectrum(OBSERVE_PHI, OBSERVE_THETA, start_theta_x, stop_theta_x, d_theta_x, start_theta, stop_theta, d_theta, Ws_Min, dW,numpar,1.0);
        else if(Get_Chirp_Spectrum(OBSERVE_PHI, OBSERVE_THETA, start_theta_x, stop_theta_x, d_theta_x, start_theta, stop_theta, d_theta, Ws_Min, dW,numpar,1.0))
        {
                cout<<endl<<"WARNING: No capability for calculating the spectrum with a chriped laser pulse outside of the plane wave approximation."<<endl
                <<"Plane wave approximation has been invoked."<<endl<<endl;
        }
        Write_Spectrum (outfile);
    }
    else if (SPECTRUM && SOLID_ANGLE && !INTENSITY && !CART_INTENSITY && !ANGLE_SPECTRUM)
    {
        Solid_Spectrum(SOLID_THETA_MIN, SOLID_THETA_MAX, SOLID_NTHETA, SOLID_NPHI, start_theta_x, stop_theta_x, d_theta_x, start_theta, stop_theta, d_theta, Ws_Min, dW,numpar);
        Write_Spectrum (outfile);
    }
    else if ( INTENSITY && !CART_INTENSITY && !ANGLE_SPECTRUM )
    {
        intenfile.open(INTENFILE);
        if (FAST_INTENSITY)
        {
            Get_Fast_Profile_Plane (PHI_NUMBER, THETA_NUMBER, THETA_RANGE,numpar);
        }
        else
        {
            Get_Profile (PHI_NUMBER, THETA_NUMBER, THETA_RANGE,start_theta_x, stop_theta_x, d_theta_x, start_theta, stop_theta, d_theta, Ws_Min, dW,numpar);
        }
        Write_Profile (intenfile, PHI_NUMBER, THETA_NUMBER);
        intenfile.close();
    }
    else if (CART_INTENSITY && !ANGLE_SPECTRUM)
    {
        if(DATA_DUMP && !FAST_INTENSITY)
            intenfile.open(INTENFILE2);
        else
            intenfile.open(INTENFILE);
        if(FAST_INTENSITY)
        {
            Get_Fast_Cart_Profile(THETA_NUMBER_X, THETA_RANGE_X, THETA_NUMBER_Y, THETA_RANGE_Y, start_theta_x, stop_theta_x, d_theta_x, start_theta, stop_theta, d_theta, numpar);
        }
        else
        {
            Get_Cart_Profile (THETA_NUMBER_X, THETA_RANGE_X, THETA_NUMBER_Y, THETA_RANGE_Y, start_theta_x, stop_theta_x, d_theta_x, start_theta, stop_theta, d_theta, Ws_Min, dW,numpar);
        }
        Write_Cart_Profile (intenfile, THETA_NUMBER_X, THETA_NUMBER_Y);
        intenfile.close();
    }
    else if (ANGLE_SPECTRUM)
    {
        Get_Angle_Spectrum_All(ANGLE_SPECTRUM_START, ANGLE_SPECTRUM_STOP, ANGLE_SPECTRUM_FIXED, ANGLE_SPECTRUM_NUM, start_theta_x, stop_theta_x, d_theta_x, start_theta, stop_theta, d_theta, Ws_Min, dW,numpar,1.0);

        Write_Angle_Spectrum(ANGLE_SPECTRUM_START, ANGLE_SPECTRUM_STOP, ANGLE_SPECTRUM_FIXED, ANGLE_SPECTRUM_NUM);
    }
    else
    {
        outfile.close();
        return(0);
    }


    SPOT_SIZE_X = SPOT_SIZE_2_AVG_X/SPOT_SIZE_NORM - SPOT_SIZE_AVG_X*SPOT_SIZE_AVG_X/(SPOT_SIZE_NORM*SPOT_SIZE_NORM);
    SPOT_SIZE_Y = SPOT_SIZE_2_AVG_Y/SPOT_SIZE_NORM - SPOT_SIZE_AVG_Y*SPOT_SIZE_AVG_Y/(SPOT_SIZE_NORM*SPOT_SIZE_NORM);

    if(SPOT_SIZE_X > 0)
    {
        SPOT_SIZE_X = (LBASE/1.0e-6)*sqrt(SPOT_SIZE_X);
        cout<<"Rms x-ray source size in X: "<<SPOT_SIZE_X<<" micro_meters."<<endl;
    }
    else
    {
        cout<<"Error computing x-ray source size in X."<<endl;
    }

    if(SPOT_SIZE_Y > 0)
    {
        SPOT_SIZE_Y = (LBASE/1.0e-6)*sqrt(SPOT_SIZE_Y);
        cout<<"Rms x-ray source size in Y: "<<SPOT_SIZE_Y<<" micro_meters."<<endl;
    }
    else
    {
        cout<<"Error computing x-ray source size in Y."<<endl;
    }

    


    outfile.close();

    return(0);

}//main

/******************************************************************************************/
int Init_SpotSize_Vars()
{
        SPOT_SIZE_AVG_X = 0.;
        SPOT_SIZE_NORM = 0.;
        SPOT_SIZE_2_AVG_X = 0.;

        SPOT_SIZE_AVG_Y = 0.;
        SPOT_SIZE_2_AVG_Y = 0.;

return(0);
}

/******************************************************************************************/
int Write_Profile (ofstream& outfile, int n_phi, int n_theta)
/* Writes the calculated intensity profile to a file */
{
    int j,k;
    prec phi, theta;
    prec dphi,dtheta;
    prec Total_Energy = 0.;


    outfile.setf(LEFT);


    outfile<<"Intensity in units of keV/mrad^2"<<endl<<endl;

    outfile<<setw(W2)<<"THETA(mrad)";

    if(n_phi>1) dphi=PHI[1]-PHI[0];
    else dphi=PI;

    dtheta = THETA[1] - THETA[0];

    for (j=0;j<n_phi;j++)
    {
        phi = PHI[j]*180./PI;
        outfile<<setw(W2)<<phi;
    }
    outfile<<endl<<endl;
    for(k=0;k<n_theta;k++)
    {
        theta = THETA[k]*1000.0;
        outfile<<PR3<<setw(W2)<<theta;
        for(j=0;j<n_phi;j++)
        {
            outfile<<PR2<<setw(W2)<<E_SCALE*(U_X[j][k])/1.0e6;
            if(sqrt(THETA[k]*THETA[k])<ZERO)
            {
                Total_Energy=Total_Energy + E_SCALE*U_X[j][k]*dphi*sin(dtheta*0.5)*sin(dtheta*0.5);
            }
            else
            {
            Total_Energy=Total_Energy+E_SCALE*U_X[j][k]*dphi*dtheta*sqrt(sin(THETA[k])*sin(THETA[k]));
            }
        }//for j
        outfile<<endl;
    }//for k

    cout<<"*************Total Energy = "<<Total_Energy<<" keV ******************"<<endl;
    cout<<"*************Valid if phi_range = 180 degrees ******************"<<endl;

    return(0);

}//Write_Profile


/******************************************************************************************/
int Write_Cart_Profile (ofstream& outfile, int n_theta_x, int n_theta_y)
/* Writes the calculated intensity profile to a file */
{
    int j,k;
    prec theta_x, theta_y;
    prec Total_E;
    prec d_theta2;

    Total_E = 0.0;
    outfile.setf(LEFT);

    outfile<<setw(W2)<<"Theta_x(mrad)"<<setw(W2)<<"Theta_y(mrad)"<<"Inensity(keV/mrad^2)"<<" "<<n_theta_x<<" "<<n_theta_y;
    outfile<<endl<<endl;

    d_theta2 = (PHI[1] - PHI[0])*(THETA[1]-THETA[0])*1.0e6;

    for(j=0;j<n_theta_x;j++)
    {
        theta_x = PHI[j]*1000.0;
        for(k=0;k<n_theta_y;k++)
        {
            theta_y = THETA[k]*1000.0;
            outfile<<PR3<<setw(W2)<<theta_x<<setw(W2)<<theta_y<<E_SCALE*(U_X[j][k])/1.0e6;
            outfile<<endl;
            Total_E = Total_E + E_SCALE*(U_X[j][k])/1.0e6*d_theta2;
        }//for j
    }//for k

    cout<<"*************Total Energy = "<<Total_E<<" keV ******************"<<endl;

    return(0);

}//Write_Cart_Profile


/******************************************************************************************/
int Write_Spectrum (ofstream& outfile)
/* Writes the calculated x-ray spectrum to a file */
{
    int i,j;
    prec W,dt,dw,Total_Photons;
    prec Total_Nx[MAX_XRAY_BIN];
    prec IvTime[MAX_TIME_BIN];
    prec NvTime[MAX_TIME_BIN];
    ofstream outfile2,outfile3,outfile4;
    prec spot_size_x,spot_size_y;

    outfile3.open("Graph_Spectrum.txt");
    outfile3.setf(LEFT);

    outfile2.open("Int_Spectrum.txt");
    outfile4.open("Int_Time.txt");
    outfile2.setf(LEFT);
    outfile4.setf(LEFT);


    spot_size_x = SPOT_SIZE_2_AVG_X/SPOT_SIZE_NORM - SPOT_SIZE_AVG_X*SPOT_SIZE_AVG_X/(SPOT_SIZE_NORM*SPOT_SIZE_NORM);
    spot_size_y = SPOT_SIZE_2_AVG_Y/SPOT_SIZE_NORM - SPOT_SIZE_AVG_Y*SPOT_SIZE_AVG_Y/(SPOT_SIZE_NORM*SPOT_SIZE_NORM);


    if(spot_size_x > 0)
    {
        spot_size_x = (LBASE/1.0e-6)*sqrt(spot_size_x);
    }
    else
    {
        spot_size_x = -1;
        cout<<"Error computing x-ray source size in X."<<endl;
    }

    if(spot_size_y > 0)
    {
        spot_size_y = (LBASE/1.0e-6)*sqrt(spot_size_y);
    }
    else
    {
        spot_size_y = -1;
        cout<<"Error computing x-ray source size in Y."<<endl;
    }


    Total_Photons = 0.;
    dt=FIND_TIME_FROM_BIN(1)-FIND_TIME_FROM_BIN(0);
    dw=E_X[1]-E_X[0];
    for (i=0; i<MAX_XRAY_BIN; i++)
        Total_Nx[i]=0.;

    for (i=0; i<MAX_TIME_BIN; i++)
    {
        IvTime[i]=0.;
        NvTime[i]=0.;
    }//for

    if (!SOLID_ANGLE)
        outfile<<"             Spectral Density (photons/mrad^2/s/eV)"<<endl<<endl;
    else
        outfile<<"             Spectral Density (photons/s/eV)"<<endl<<endl;


    outfile<<setw(W2)<<"Energy(keV)"<<setw(W2)<<"Time(ps)";

    if (!SOLID_ANGLE)
        outfile3<<setw(W2)<<"Time(ps)"<<setw(W2)<<"Energy(keV)"<<"photons/mrad^2/s/eV"<<"  "<<TIME_BIN_NUM<<"  "<<N_Ws<<endl<<endl;
    else
        outfile3<<setw(W2)<<"Time(ps)"<<setw(W2)<<"Energy(keV)"<<"photons/s/eV"<<"  "<<TIME_BIN_NUM<<"  "<<N_Ws<<endl<<endl;

    for (j=0; j< TIME_BIN_NUM; j++)
    {
            outfile<<setw(W2)<<PR3<<FIND_TIME_FROM_BIN(j)*TBASE/OUT_TIME_UNIT;
    }//for i
    outfile<<endl<<endl;

    if (!SOLID_ANGLE)
        outfile2<<setw(W2)<<"Energy(keV)"<<"Spectral_Density(Photons/mrad^2/eV)"<<"  "<<spot_size_x<<"  "<<spot_size_y<<endl<<endl;
    else
        outfile2<<setw(W2)<<"Energy(keV)"<<"Spectral_Density(Photons/eV)"<<"  "<<spot_size_x<<"  "<<spot_size_y<<endl<<endl;

    for (i=0; i<N_Ws; i++)
    {
        W=E_SCALE*E_X[i];
        outfile<<setw(W2)<<W<<setw(W2)<<" ";
        outfile2<<setw(W2)<<W;
        for (j=0; j<TIME_BIN_NUM; j++)
        {
            if (!SOLID_ANGLE)
            {
                outfile<<setw(W2)<<PR2<<N_X[j][i]/E_SCALE/1.0e9/dt/TBASE;  //photons/mrad^2/s/eV
                Total_Nx[i]=Total_Nx[i]+N_X[j][i]/E_SCALE/1.0e9; //photons/mrad^2/eV
            }
            else
            {
                outfile<<setw(W2)<<PR2<<N_X[j][i]/E_SCALE/1.0e3/dt/TBASE;  //photons/s/eV
                Total_Nx[i]=Total_Nx[i]+ (N_X[j][i]/E_SCALE)/1.0e3; //photons/eV
            }
        }//for j
        if(SOLID_ANGLE) Total_Photons = Total_Photons+Total_Nx[i]*dw*E_SCALE*1.0e3;
        outfile2<<Total_Nx[i]<<endl;
        outfile<<endl;
    }//for i


    if (!SOLID_ANGLE)
        outfile4<<setw(W2)<<"Time(ps)"<<setw(W3)<<"Intensity(photons/s/mrad^2)"<<"Intensity(keV/s/mrad^2)"<<endl<<endl;
    else
    {
        outfile4<<setw(W2)<<"Time(ps)"<<setw(W3)<<"Intensity(photons/s)"<<"Intensity(keV/s)"<<endl<<endl;
        cout<<"Total Photons = "<<Total_Photons<<endl;
    }

    for (j=0; j<TIME_BIN_NUM; j++)
    {
        outfile4<<setw(W2)<<PR2<<FIND_TIME_FROM_BIN(j)*TBASE/OUT_TIME_UNIT;
        for (i=0; i<N_Ws; i++)
        {
            if (!SOLID_ANGLE)
            {
                NvTime[j]=NvTime[j]+N_X[j][i]*dw/dt/1.0e6/TBASE;  //photons/s/mrad^2
                IvTime[j]=IvTime[j]+(E_SCALE)*E_X[i]*N_X[j][i]*dw/dt/1.0e6/TBASE; //keV/time_unit/mrad^2
                W=E_SCALE*E_X[i];
                outfile3<<setw(W2)<<PR2<<FIND_TIME_FROM_BIN(j)*TBASE/OUT_TIME_UNIT<<setw(W2)<<W<<N_X[j][i]/E_SCALE/1.0e9/dt/TBASE<<endl;
            }
            else
            {
                NvTime[j]=NvTime[j]+N_X[j][i]*dw/dt/TBASE;  //photons/s
                IvTime[j]=IvTime[j]+(E_SCALE)*E_X[i]*N_X[j][i]*dw/dt/TBASE; //keV/time_unit/mrad^2
                W=E_SCALE*E_X[i];
                outfile3<<setw(W2)<<PR2<<FIND_TIME_FROM_BIN(j)*TBASE/OUT_TIME_UNIT<<setw(W2)<<W<<N_X[j][i]/E_SCALE/1.0e3/dt/TBASE<<endl;
            }
        }//for i
        outfile4<<PR2<<setw(W3)<<NvTime[j]<<IvTime[j]<<endl;
    }


    outfile2.close();
    outfile3.close();
    outfile4.close();

    return(0);

}//Write_Spectrum
/******************************************************************************************/
prec FILTER(prec freq, prec width, prec center)
/* This peformes an attenuation of the spectrum at a given wavelength
to simulate a filter function. Width is the FWHM of the distribution.
  */
{
    prec ans;
    prec dw2;
    prec width2;

    if(width<ZERO) return(1.);

    dw2=(freq-center)*(freq-center);

/*  width2 = 0.5*width*width;
    ans = exp(-dw2/width2);  */

    width2 = 0.25*width*width;
    ans = FILTER_AMP/(1.0+dw2/width2);

    return(ans);


}//FILTER
/******************************************************************************************/
int Get_Total_Photon (int numpar,prec& total_dose)
/* This procedure calculates the total number of scattered photons integrated
over all wavelenghths and solid angles. Assumes planewave, and time integration
assumes non-chipred x-ray pulse. Wavelength filters are not effective in this calculation.*/
{
    prec Xee[6];
    prec start_time,stop_time,dt,ing,time_zero;
    int i,t_bin;
    prec IvTime[MAX_TIME_BIN];
    prec time;

    ofstream out;

    out.open("total_flux.txt");
    out.setf(LEFT);


    dt=FIND_TIME_FROM_BIN(1)-FIND_TIME_FROM_BIN(0);
    prec total=0.0;

    for (i=0; i<MAX_TIME_BIN; i++)
        IvTime[i]=0.;

    for (i=0; i<numpar;i++)
    {
       COPYARRAY(X[i],Xee,6);


       if(Get_Start_Stop_Time(Xee,start_time,stop_time))
       {
           ing = 0.0;
       }
       else
            ing = Integral_ng(Xee,start_time, stop_time);

        SPOT_SIZE_AVG_X = SPOT_SIZE_AVG_X + ing*Xee[0];
        SPOT_SIZE_AVG_Y = SPOT_SIZE_AVG_Y + ing*Xee[1];
        SPOT_SIZE_2_AVG_X = SPOT_SIZE_2_AVG_X + ing*Xee[0]*Xee[0];
        SPOT_SIZE_2_AVG_Y = SPOT_SIZE_2_AVG_Y + ing*Xee[1]*Xee[1];
        SPOT_SIZE_NORM = SPOT_SIZE_NORM +ing;


       total = total + 8.0*PI/3.0*ing;
       //time_zero = Xee[2]/cc;
       time_zero = (start_time+stop_time)*0.5;
       t_bin=FIND_TIME_BIN(time_zero,Xee,0.0,0.0);
       IvTime[t_bin]=IvTime[t_bin] + 8.0*PI/3.0*ing/dt;
       if((i%NUMSKIP)==0)
       {
        cout<<i<<endl;
        cout<<"ing = "<<ing<<endl;
       }

    }//for

    for (i=0; i< TIME_BIN_NUM;i++)
    {
        time = FIND_TIME_FROM_BIN(i);
        out<<setw(15)<<PR2<<time<<IvTime[i]<<endl;
    }

    cout<<endl<<"Total Number of Scattered Photons: "<<PR2<<total<<endl;
    total_dose = total;
    out.close();

    return(0);

}//Get_Total_Photon
/******************************************************************************************/
int Get_Fast_Profile_Plane (int n_phi, int n_theta, prec theta_range, int numpar)
/* Calculates the angular intensity profile of the x-ray beam integrated over all wavelengths. Assumes the
laser is a plane wave (i.e., does not sum over initial k_perp's). n_phi is the number of phi observation angles
from 0 to 2*pi, n_theta is the number of theta observation angles from -theta_range to theta_range. The integration
over scattred wavelengths in this procedure is analytic, and does not include in filtering effects. */
{
    prec Xee[6];
    prec ing;
    int i,j,k;
    prec phi_obs, theta_obs;
    prec dphi;
    prec dtheta;
    prec phi_adj;
    prec theta_adj;
    prec mid_time;

/*  phi_obs = 0.;
    theta_obs = -theta_range;
    dphi = PI/prec(n_phi);
    dtheta = 2.0*theta_range/(n_theta);
    phi_obs =0.;
    theta_obs = -theta_range;*/

    dphi = (PHI_RANGE_MAX-PHI_RANGE_MIN)/prec(n_phi);
    dtheta = 2.0*theta_range/prec(n_theta-1);
    if ((PHI_RANGE_MAX-PHI_RANGE_MIN)< (PI-ZERO))
    {
        dphi = (PHI_RANGE_MAX-PHI_RANGE_MIN)/prec(n_phi-1);
    }
    phi_obs =PHI_RANGE_MIN;
    theta_obs = -theta_range;


    for (j=0;j<n_phi;j++)
        PHI[j] = phi_obs + j*dphi;

    for(k=0;k<n_theta;k++)
        THETA[k] = theta_obs + k*dtheta;

    for (i=0; i<numpar; i++)
    {
       COPYARRAY(X[i],Xee,6);

       Get_Time_Integraion_2(Xee, ing, mid_time, i);

       for (j=0;j<n_phi;j++)
       {
           for(k=0;k<n_theta;k++)
           {
                phi_adj = PHI[j];
                theta_adj = THETA[k];
                Alter_Observation(Xee, phi_adj, theta_adj, mid_time);
                U_X[j][k]=U_X[j][k] + CrossSection(Xee, phi_adj, theta_adj, THETA_0, 0)*ing*H(Xee,THETA_0, 0.0, phi_adj, theta_adj);
           }//for
       }//for
       if((i%NUMSKIP)==0)
       {
        cout<<i<<endl;
        cout<<"ing = "<<ing<<endl;
       }
    }//for


    return(0);

}//Get_Profile_Plane

/******************************************************************************************/
int Get_Profile (int n_phi, int n_theta, prec theta_range, prec start_thetax, prec stop_thetax, prec d_thetax, prec start_thetay, prec stop_thetay, prec d_thetay, prec Ws_Min, prec dW, int numpar)
/* Calculates the angular intensity profile of the x-ray beam integrated over all wavelengths.
 n_phi is the number of phi observation angles from 0 to 2*pi, n_theta is the number of theta observation
 angles from -theta_range to theta_range. */

{
    int j,k;

    prec phi_obs, theta_obs;
    prec dphi;
    prec dtheta;

    dphi = (PHI_RANGE_MAX-PHI_RANGE_MIN)/prec(n_phi);
    dtheta = 2.0*theta_range/prec(n_theta-1);
    if ((PHI_RANGE_MAX-PHI_RANGE_MIN)< (PI-ZERO))
    {
        dphi = (PHI_RANGE_MAX-PHI_RANGE_MIN)/prec(n_phi-1);
    }
    phi_obs =PHI_RANGE_MIN;
    theta_obs = -theta_range;

    for (j=0;j<n_phi;j++)
        PHI[j] = phi_obs + j*dphi;

    for(k=0;k<n_theta;k++)
        THETA[k] = theta_obs + k*dtheta;

    for (j=0;j<n_phi;j++)
    {
        for(k=0;k<n_theta;k++)
        {	
			if(!CHIRP_LASER)
				U_X[j][k]=U_X[j][k] + Get_Intensity (PHI[j], THETA[k], start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay, Ws_Min, dW,numpar);
			else
				U_X[j][k]=U_X[j][k] + Get_Chirp_Intensity (PHI[j], THETA[k], start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay, Ws_Min, dW,numpar);								
            cout<<"Phi = "<<PHI[j]<<"   Theta = "<<THETA[k]<<endl;
        }//for
    }//for

    return(0);

}

/******************************************************************************************/
int Get_Cart_Profile (int n_theta_x, prec theta_range_x, int n_theta_y, prec theta_range_y, prec start_thetax, prec stop_thetax, prec d_thetax, prec start_thetay, prec stop_thetay, prec d_thetay, prec Ws_Min, prec dW, int numpar)
/* Calculates the angular intensity profile of the x-ray beam integrated over all wavelengths.
 n_phi is the number of phi observation angles from 0 to 2*pi, n_theta is the number of theta observation
 angles from -theta_range to theta_range. */
{
    int j,k;

    prec theta_x_obs, theta_y_obs;
    prec dtheta_x;
    prec dtheta_y;
    prec theta_obs;
    prec phi_obs;

    phi_obs = 0.;
    theta_x_obs = -theta_range_x;
    theta_y_obs = -theta_range_y;
    dtheta_x = 2.0*theta_range_x/prec(n_theta_x-1);
    dtheta_y = 2.0*theta_range_y/prec(n_theta_y-1);

    prec DATA_DUMP_DATA[MAX_DATA_DUMP_FREQ_STEPS];

    ofstream dump;

    if (DATA_DUMP)
    {
        dump.open(DATA_DUMP_FILE);
        dump.setf(LEFT);
        dump<<setw(W2)<<"Theta_x(mrad)"<<setw(W2)<<"Theta_y(mrad)"<<"Inensity(keV/mrad^2/eV)"<<endl;
        dump<<setw(W2)<<n_theta_x<<setw(W2)<<n_theta_y<<N_Ws<<" "<<E_SCALE*Ws_Min<<" "<<E_SCALE*Ws_Max<<endl<<endl;
    }//if

    for (j=0;j<n_theta_x;j++)
        PHI[j] = theta_x_obs + j*dtheta_x;  //this stores theta_x values

    for(k=0;k<n_theta_y;k++)
        THETA[k] = theta_y_obs + k*dtheta_y; //this stores theta_y values

    for (j=0;j<n_theta_x;j++)
    {
        for(k=0;k<n_theta_y;k++)
        {
            phi_obs = atan2(THETA[k],PHI[j]);
            theta_obs = sqrt(THETA[k]*THETA[k] + PHI[j]*PHI[j]);
            if (DATA_DUMP==0)
				if(!CHIRP_LASER)
					U_X[j][k]=U_X[j][k] + Get_Intensity (phi_obs, theta_obs, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay, Ws_Min, dW,numpar);
				else
					U_X[j][k]=U_X[j][k] + Get_Chirp_Intensity (phi_obs, theta_obs, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay, Ws_Min, dW,numpar);				
            cout<<"Theta_x = "<<PHI[j]<<"   Theta_y = "<<THETA[k]<<endl;

            if (DATA_DUMP)
            {
				if(!CHIRP_LASER)
				{
					U_X[j][k] = Get_Data_Dump_Intensity (phi_obs, theta_obs, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay, Ws_Min, dW,numpar, DATA_DUMP_DATA);
					WRITE_DATA_DUMP (dump, DATA_DUMP_DATA, Ws_Min, Ws_Max, N_Ws, PHI[j], THETA[k]);
				}
				else
				{
					U_X[j][k] = Get_Chirp_Data_Dump_Intensity (phi_obs, theta_obs, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay, Ws_Min, dW,numpar, DATA_DUMP_DATA);
					WRITE_DATA_DUMP (dump, DATA_DUMP_DATA, Ws_Min, Ws_Max, N_Ws, PHI[j], THETA[k]);
				}
            }//if
        }//for
    }//for

    if (DATA_DUMP) dump.close();
    return(0);
}

/******************************************************************************************/
int Get_Fast_Cart_Profile (int n_theta_x, prec theta_range_x, int n_theta_y, prec theta_range_y, prec start_thetax, prec stop_thetax, prec d_thetax, prec start_thetay, prec stop_thetay, prec d_thetay, int numpar)
/* Calculates the angular intensity profile of the x-ray beam integrated over all wavelengths. Assumes the
laser is a plane wave (i.e., does not sum over initial k_perp's). n_theta_x is the number of theta observation
angles from -theta_range_x to theta_range_x. Likewise of n_theta_y. The integration over scattred wavelengths
in this procedure is analytic, and does not include in filtering effects. */
{
    prec Xee[6];
    prec ing;
    int i,j,k;
    prec phi_adj;
    prec theta_adj;
    prec mid_time;
    prec theta_x_obs, theta_y_obs;
    prec dtheta_x;
    prec dtheta_y;

    theta_x_obs = -theta_range_x;
    theta_y_obs = -theta_range_y;
    dtheta_x = 2.0*theta_range_x/prec(n_theta_x-1);
    dtheta_y = 2.0*theta_range_y/prec(n_theta_y-1);


    for (j=0;j<n_theta_x;j++)
        PHI[j] = theta_x_obs + j*dtheta_x;  //this stores theta_x values

    for(k=0;k<n_theta_y;k++)
        THETA[k] = theta_y_obs + k*dtheta_y; //this stores theta_y values

    for (i=0; i<numpar; i++)
    {
       COPYARRAY(X[i],Xee,6);

/*     start_time = -0.5*INTEGRAL_RANGE_SIGMA_TIME*PULSE_WIDTH/cc -0.5*Xee[2]/cc;
       stop_time = 0.5*INTEGRAL_RANGE_SIGMA_TIME*PULSE_WIDTH/cc -0.5*Xee[2]/cc;
       ing = Integral_ng(Xee,start_time, stop_time);*/

       Get_Time_Integraion_2(Xee, ing, mid_time, i);

       for (j=0;j<n_theta_x;j++)
       {
           for(k=0;k<n_theta_y;k++)
           {
                phi_adj = atan2(THETA[k],PHI[j]);
                theta_adj = sqrt(THETA[k]*THETA[k] + PHI[j]*PHI[j]);
                Alter_Observation(Xee, phi_adj, theta_adj, mid_time);
                U_X[j][k]=U_X[j][k] + CrossSection(Xee, phi_adj, theta_adj, THETA_0, 0)*ing*H(Xee,THETA_0, 0.0, phi_adj, theta_adj);
           }//for
       }//for
       if((i%NUMSKIP)==0)
       {
        cout<<i<<endl;
        cout<<"ing = "<<ing<<endl;
       }
    }//for

    return(0);

}//Get_Fast_Cart_Profile

/******************************************************************************************/

int WRITE_DATA_DUMP (ofstream & dump, prec data[], prec w_start, prec w_stop, int n, prec theta_x, prec theta_y)
{
    int i;


    theta_x = theta_x*1000.;
    theta_y = theta_y*1000.;

    dump<<PR2<<setw(W2)<<theta_x<<setw(W2)<<theta_y;

    for (i=0; i<(n-1); i++)
        dump<<PR2<<setw(W2)<<(data[i])/1.0e9;
    dump<<(data[n-1])/1.0e9<<endl;

    return(0);

}//WRITE_DATA_DUMP

/******************************************************************************************/
int Solid_Spectrum (prec t_min, prec t_max, int n_theta, int n_phi, prec start_thetax, prec stop_thetax, prec d_thetax, prec start_thetay, prec stop_thetay, prec d_thetay, prec Ws_Min, prec dW, int numpar)
/* Integrates the scattered frequency spectrum over a solid angle specified by t_min, t_max for the theta range over
2*PI in phi. */
{
    prec int_factor;
    prec d_theta, d_phi, theta, phi;
    int i,j;

    d_phi = 2*PI/prec(n_phi);
    d_theta = (t_max-t_min)/prec(n_theta);

/*  if ( (t_min<ZERO) && (t_min<d_theta) )
    {
        phi = 0.0;
        theta=0.0;
        int_factor = PI*sin(d_theta)*sin(d_theta);
        Get_Spectrum(phi, theta, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay, Ws_Min, dW,numpar,int_factor);
        t_min = t_min + d_theta;
        n_theta = n_theta - 1;
    }*/

    theta = t_min+d_theta*0.5;

    for(i=0; i<n_theta; i++)
    {
        int_factor = sin(theta)*d_phi*d_theta;
        for (j=0; j<n_phi; j++)
        {
           phi = prec(j)*d_phi;
           if(!CHIRP_LASER)
            Get_Spectrum(phi, theta, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay, Ws_Min, dW,numpar,int_factor);
           else
            Get_Chirp_Spectrum(phi, theta, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay, Ws_Min, dW,numpar,int_factor);

           cout<<"theta: "<<i<<"  phi: "<<j<<endl;
        }//for j
        theta = theta + d_theta;
    }// for i

//  if( CHIRP_LASER && !PLANE_WAVE) return(1);
    return(0);

}//Solid_Spectrum

/******************************************************************************************/
prec Get_Data_Dump_Intensity (prec phi0, prec theta0, prec start_thetax, prec stop_thetax, prec d_thetax, prec start_thetay, prec stop_thetay, prec d_thetay, prec Ws_Min, prec dW, int numpar, prec data[])
/* Calculates the scattered intensity at given observation direction by summing over the spectrum produced by
all the macro-particles in the simulation. If PLANE_WAVE = 1, assumes laser is a plane wave (i.e., does not sum over
initial k_perp's). If PLANE_WAVE=0, sums over initial k-vector states (drastically slows done the calculation). */
{
    prec W;
    prec dw;
    prec w_start;
    prec w_stop;
    prec mid_time;
    prec ing;
    prec Xee[6];
    prec phi,theta;
    prec total = 0.;

    int i,j;
    //dw = dW;
    //w_start = Ws_Min;

    /* Initialized data[] array */

    for (i = 0; i<N_Ws; i++)
        data[i] =0.;


    for (i=0; i<numpar; i++)
    {
       phi = phi0;
       theta=theta0;
       COPYARRAY(X[i],Xee,6);

       Get_Time_Integraion(Xee, ing, mid_time, phi, theta, i);

       Find_W_Range(Xee, phi, theta, THETA_0, w_start, dw);
       w_stop = w_start + prec(N_Ws-1)*dw;
        for (j=0;j<N_Ws;j++)
        {
            W=Ws_Min+prec(j)*dW;
            if( (W > w_start) && (W < w_stop)  )
            {
                if(PLANE_WAVE)
                    data[j] = data[j]+ W*CrossSection(Xee, phi, theta, THETA_0, 0)*OMEGA_FACTOR(Xee, THETA_0, 0, phi, theta, W, BANDWIDTH)*ing;
                else
                    data[j] =data[j] + W*Integral_dN_x_y(Xee,0,0,W,BANDWIDTH, ing, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay,0.,0.);
            }//if
        }//for j
    }// for i


    for (i=0; i<N_Ws; i++)
        total = total + dW*data[i];


    if(FIRST_TIME)
        FIRST_TIME = 0;

    return(total);

}//prec Get_Data_Dump_Intensity


/******************************************************************************************/
prec Get_Intensity (prec phi0, prec theta0, prec start_thetax, prec stop_thetax, prec d_thetax, prec start_thetay, prec stop_thetay, prec d_thetay, prec Ws_Min, prec dW, int numpar)
/* Calculates the scattered intensity at given observation direction by summing over the spectrum produced by
all the macro-particles in the simulation. If PLANE_WAVE = 1, sssumes laser is a plane wave (i.e., does not sum over
initial k_perp's). If PLANE_WAVE=0, sums over initial k-vector states (drastically slows done the calculation). */
{

    prec W;
    prec dw;
    prec w_start;
    prec U=0.;
    prec mid_time;
    prec ing;
    prec Xee[6];
    prec phi,theta;

    int i,j;

    //dw = dW;
    //w_start = Ws_Min;

    for (i=0; i<numpar; i++)
    {
       phi = phi0;
       theta=theta0;
       COPYARRAY(X[i],Xee,6);

       Get_Time_Integraion(Xee, ing, mid_time, phi, theta, i);

       Find_W_Range(Xee, phi, theta, THETA_0, w_start, dw);

       for (j=0;j<N_Ws;j++)
        {
            W=w_start+prec(j)*dw;
            if(PLANE_WAVE)
                U = U + dw*W*CrossSection(Xee, phi, theta, THETA_0, 0)*OMEGA_FACTOR(Xee, THETA_0, 0, phi, theta, W, BANDWIDTH)*ing;
            else
                U = U + dw*W*Integral_dN_x_y(Xee,0,0,W,BANDWIDTH, ing, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay,0.,0.);
        }//for j

       if((i%NUMSKIP)==0)
       {
         if(FIRST_TIME)
         {
            cout<<i<<endl;
            cout<<"ing = "<<ing<<endl;
         }
       }
    }// for i

    if(FIRST_TIME)
        FIRST_TIME = 0;

    return(U);

}//prec Get_Intensity

/******************************************************************************************/
prec Get_Chirp_Intensity (prec phi_0, prec theta_0, prec start_thetax, prec stop_thetax, prec d_thetax, prec start_thetay, prec stop_thetay, prec d_thetay, prec Ws_Min, prec dW, int numpar)
{
    prec start_time;
    prec stop_time;
    prec mid_time;
    prec dt;
    prec time_zero;
    prec Xee[6];
    prec ing;
    prec phi,theta;
    prec time;
    int i,j,k;
    prec W;
	prec U=0.;
    int t_bin;
    int t_steps;
    prec CS;
    prec wstart[MAX_TIME_BIN];
    prec wstop[MAX_TIME_BIN];


    for (i=0; i<numpar; i++)
    {
        phi=phi_0;
        theta=theta_0;
        COPYARRAY(X[i],Xee,6);

        if(Get_Start_Stop_Time(Xee,start_time,stop_time))
        {
           ing = 0.0;
        }
        else
        {
           mid_time = (start_time+stop_time)*0.5;
           Alter_Observation(Xee, phi, theta, mid_time);
           dt = (stop_time-start_time)/prec(T_STEPS);
           t_steps = int( (stop_time-start_time)/dt + 0.5);
           CS = CrossSection(Xee, phi, theta, THETA_0, 0);
           if (PLANE_WAVE)
            Find_W_Range_Chirp(Xee, THETA_0, 0.0, phi, theta, BANDWIDTH, start_time, dt, t_steps, wstart, wstop);
           else
            Find_W_Range_Chirp_Kp(Xee, start_thetax, stop_thetax, start_thetay, stop_thetay, phi, theta, BANDWIDTH, start_time, dt, t_steps, wstart, wstop);
           for (j=0;j<N_Ws;j++)
           {
                W=Ws_Min+prec(j)*dW;
                time = start_time;
                for (k=0;k<t_steps;k++)
                {
                   if ( (W >= wstart[k]) && (W<=wstop[k]) )
                   {
                    ing = Integral_ngnw(Xee,time,dt,THETA_0,0.0,phi,theta,W,BANDWIDTH, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay);
                    time_zero=time+0.5*dt;
                    t_bin=FIND_TIME_BIN(time_zero,Xee,phi,theta);
                    if (PLANE_WAVE)
                        U = U + dW*W*CS*ing;
                    else
                    {
                       // U = U + dw*W*ing;//Use this if CS is calculated in dN_dw_domega
						U = U + dW*W*CS*ing;
                    }
                   }
                   time = time+dt;
                }//for k
            }//for j
        }//else

       if((i%NUMSKIP_CHIRP)==0)
       {
          if(FIRST_TIME)
		  {
			cout<<i<<endl;
			cout<<"ing = "<<ing<<endl;
		  }
       }
    }//for i

  if(FIRST_TIME)
  {
	FIRST_TIME = 0;
	for (j=0;j<N_Ws;j++)
		{
			W=Ws_Min+prec(j)*dW;
			E_X[j]=W;
		}//for
  }

    return(U);

}//prec Get_Chirp_Intensity

/******************************************************************************************/
prec Get_Chirp_Data_Dump_Intensity (prec phi_0, prec theta_0, prec start_thetax, prec stop_thetax, prec d_thetax, prec start_thetay, prec stop_thetay, prec d_thetay, prec Ws_Min, prec dW, int numpar, prec data[])
{
    prec start_time;
    prec stop_time;
    prec mid_time;
    prec dt;
    prec Xee[6];
    prec ing;
    prec phi,theta;
    prec time;
    int i,j,k;
    prec W;
    int t_steps;
    prec CS;
    prec wstart[MAX_TIME_BIN];
    prec wstop[MAX_TIME_BIN];
	prec total = 0.;


/* Initialized data[] array */

    for (i = 0; i<N_Ws; i++)
        data[i] =0.;

    for (i=0; i<numpar; i++)
    {
        phi=phi_0;
        theta=theta_0;
        COPYARRAY(X[i],Xee,6);

        if(Get_Start_Stop_Time(Xee,start_time,stop_time))
        {
           ing = 0.0;
        }
        else
        {
           mid_time = (start_time+stop_time)*0.5;
           Alter_Observation(Xee, phi, theta, mid_time);
           dt = (stop_time-start_time)/prec(T_STEPS);
           t_steps = int( (stop_time-start_time)/dt + 0.5);
           CS = CrossSection(Xee, phi, theta, THETA_0, 0);
           if (PLANE_WAVE)
            Find_W_Range_Chirp(Xee, THETA_0, 0.0, phi, theta, BANDWIDTH, start_time, dt, t_steps, wstart, wstop);
           else
            Find_W_Range_Chirp_Kp(Xee, start_thetax, stop_thetax, start_thetay, stop_thetay, phi, theta, BANDWIDTH, start_time, dt, t_steps, wstart, wstop);
           for (j=0;j<N_Ws;j++)
           {
                W=Ws_Min+prec(j)*dW;
                time = start_time;
                for (k=0;k<t_steps;k++)
                {
                   if ( (W >= wstart[k]) && (W<=wstop[k]) )
                   {
                    ing = Integral_ngnw(Xee,time,dt,THETA_0,0.0,phi,theta,W,BANDWIDTH, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay);
                    if (PLANE_WAVE)						 
                        data[j] = data[j] + W*CS*ing;
                    else
                    {
                       // data[j] = data[j] + W*ing;//Use this if CS is calculated in dN_dw_domega
						data[j] = data[j] + W*CS*ing;
                    }
                   }
                   time = time+dt;
                }//for k
            }//for j
        }//else


    }//for i


  for (i=0; i<N_Ws; i++)
        total = total + dW*data[i];
	
	
  if(FIRST_TIME)
  {
	FIRST_TIME = 0;
	for (j=0;j<N_Ws;j++)
		{
			W=Ws_Min+prec(j)*dW;
			E_X[j]=W;
		}//for
  }
  

    return(total);

}//prec Get_Chirp_Data_Dump_Intensity


/******************************************************************************************/

int Get_Time_Integraion(prec Xee[], prec& ing, prec& mid_time, prec &phi, prec &theta, int npar)
/*
    This procdure is used to store and access the time integrated component of the overlap integral
    of the electron macro-particle with the laser pulse, and also the average time of the interaction
    "mid_time". The procedure first checks to see if the time integral has already be calculated
    for this macro-particle, designated by "npar".  If not, the integral if calculated, if so, then
    the results are accessed from the array "NG" and "MT" (for mid_time). The procedure for calculating
    the alteration of the observation angle due to finite spot size effects is also called within
    this procedure.
*/
{

    prec start_time, stop_time;

    if(FIRST_TIME)
    {
        if(Get_Start_Stop_Time(Xee,start_time,stop_time))
        {
           ing = 0.0;
           mid_time = 0.0;
           MT[npar] = mid_time;
           NG[npar] = ing;
           return(1);
        }
        else
        {
            ing = Integral_ng(Xee,start_time, stop_time);
            NG[npar] = ing;
            mid_time = (start_time+stop_time)*0.5;
            MT[npar] = mid_time;
            Alter_Observation(Xee, phi, theta, mid_time);
        }
        SPOT_SIZE_AVG_X = SPOT_SIZE_AVG_X + ing*Xee[0];
        SPOT_SIZE_AVG_Y = SPOT_SIZE_AVG_Y + ing*Xee[1];
        SPOT_SIZE_2_AVG_X = SPOT_SIZE_2_AVG_X + ing*Xee[0]*Xee[0];
        SPOT_SIZE_2_AVG_Y = SPOT_SIZE_2_AVG_Y + ing*Xee[1]*Xee[1];
        SPOT_SIZE_NORM = SPOT_SIZE_NORM +ing;
    }
    else
    {
        ing = NG[npar];
        mid_time = MT[npar];
        Alter_Observation(Xee, phi, theta, mid_time);
    }

    return(0);

}//Get_Time_Integration

/******************************************************************************************/
int Get_Time_Integraion_2(prec Xee[], prec& ing, prec& mid_time, int npar)
/*
    This procdure is used to calculate the time integrated component of the overlap integral
    of the electron macro-particle with the laser pulse, and also the average time of the interaction
    "mid_time".
*/
{

    prec start_time, stop_time;


    if(Get_Start_Stop_Time(Xee,start_time,stop_time))
    {
        ing = 0.0;
        mid_time = 0.0;
        return(1);
    }
    else
    {
        ing = Integral_ng(Xee,start_time, stop_time);
        mid_time = (start_time+stop_time)*0.5;
    }
    SPOT_SIZE_AVG_X = SPOT_SIZE_AVG_X + ing*Xee[0];
    SPOT_SIZE_AVG_Y = SPOT_SIZE_AVG_Y + ing*Xee[1];
    SPOT_SIZE_2_AVG_X = SPOT_SIZE_2_AVG_X + ing*Xee[0]*Xee[0];
    SPOT_SIZE_2_AVG_Y = SPOT_SIZE_2_AVG_Y + ing*Xee[1]*Xee[1];
    SPOT_SIZE_NORM = SPOT_SIZE_NORM +ing;


    return(0);

}//Get_Time_Integration2

/******************************************************************************************/

int Get_Start_Stop_Time (prec Xee[], prec & start, prec &stop)
/* Returns the start time and stop of the time integral for the electron represented by the
coordinates in Xee[x,y,z,px,py,pz]. Represents the laser pulse as a 3D ellipsoid with axis
lengths equal to 3X the 1/e^2 pulse length and radius at the focus, and finds the intersection
of the electron bunch with this ellipsoid. It then calls a procedure to check that the
ratio between the photon density at the start and middle of the integration is at least 10^-4. If
not, the exponential fall off factor is calculated, and the start and stop time of the simulation
are accordingly adjusted about the center integration time.*/
{
    prec wx,wx2;
    prec vz, vx,vx2,vz2,L2;
    prec xe0,ze0,xe02,ze02;
    prec L;
    prec cos_t,sin_t, sin_t2, sin_t4,sin_t3, cos_t2,cos_t3,cos_t4;;
    prec arg;
    prec MapleGenVar2,MapleGenVar1;
    prec zr_factor;

    vz = (Xee[5]/(GAMMA(Xee)));
    vx =  (Xee[3]/(GAMMA(Xee)));
    if(THETA_0 == 0)
        xe0 =0;
    else
        xe0=Xee[0];

    ze0=Xee[2];

    zr_factor = sqrt(0.25*PULSE_WIDTH*PULSE_WIDTH + ze0*ze0*0.25)/ZRx;
    //zr_factor = 0.5*PULSE_WIDTH/ZRx;
    if (zr_factor < 20.0)
    {
        wx =  INTEGRAL_RANGE_SIGMA_TIME*w0x*sqrt(1+zr_factor*zr_factor);
    }
    else
    {
        wx =  INTEGRAL_RANGE_SIGMA_TIME*w0x*20;
    }

    L = INTEGRAL_RANGE_SIGMA_TIME*PULSE_WIDTH;
    sin_t = sin(THETA_0);
    cos_t = cos(THETA_0);

    wx2=wx*wx;
    vx2=vx*vx;
    vz2=vz*vz;
    L2 = L*L;
    xe02=xe0*xe0;
    ze02=ze0*ze0;
    sin_t2=sin_t*sin_t;
    sin_t4=sin_t2*sin_t2;
    sin_t3=sin_t2*sin_t;

    cos_t2=cos_t*cos_t;
    cos_t4=cos_t2*cos_t2;
    cos_t3=cos_t2*cos_t;

    arg = -1.0;
    int count = 0;
    while ((arg<0.0) && (count<2))
    {
        arg = -wx2*L2*(-2.0*cos_t*xe0*sin_t*ze0-L2*cos_t2*vx2+cos_t2*xe02+sin_t2*ze0*
        ze0+2.0*L2*cos_t*vx*sin_t*vz+2.0*sin_t3*vx*ze02+cos_t2*
        cos_t2*vz2*xe02+2.0*cos_t3*vz*xe02+cos_t3*
        cos_t*vx2*ze02-2.0*sin_t2*xe0*cos_t*vx*ze0+2.0*cos_t2*vx2*
        sin_t2*ze02+2.0*cos_t2*ze02*vx*sin_t-2.0*cos_t3*ze0*
        xe0*vx-2.0*sin_t3*ze0*vz*xe0-2.0*cos_t4*ze0*vz*xe0*
        vx+2.0*cos_t2*xe02*sin_t2*vz2-4.0*cos_t2*xe0*sin_t2*vz
        *vx*ze0-2.0*cos_t2*xe0*sin_t*vz*ze0+2.0*cos_t*xe02*sin_t2*vz-2.0*
        sin_t4*xe0*vx*ze0*vz+sin_t4*vz2*xe02-2.0
        *wx2*cos_t*vx*sin_t*vz+sin_t4*vx2*ze02-L2*sin_t2*
        vz2-wx2*cos_t2*vz2-2.0*wx2*cos_t*vz-2.0*wx2*sin_t*vx-wx2*sin_t
        *sin_t*vx2-wx2);
        if(arg < 0.0)
        {
            wx = 2.0*wx;
            L = 2.0*L;
            wx2=wx*wx;
            L2 = L*L;
        }//if
        count++;
    }//if

    if(arg < 0.0) return(1);

    MapleGenVar2 = -wx2*sin_t2*xe0*vx-wx2*cos_t2*ze0*vz+L2*cos_t*xe0*sin_t*vz-wx2*cos_t*xe0*sin_t*vz-wx2*cos_t*vx*sin_t*ze0;
    MapleGenVar1 = MapleGenVar2-wx2*cos_t*ze0-L2*sin_t2*ze0*vz-wx2*sin_t*xe0+L2*cos_t*vx*sin_t*ze0-L2*cos_t2*xe0*vx;
    MapleGenVar2 = 1.0/(L2*cos_t2*vx2+L2*sin_t2*vz2+wx2*cos_t2*vz2+2.0*wx2*cos_t*vz+wx2*sin_t2*vx2+2.0*wx2*sin_t*vx+2.0*
    wx2*cos_t*vx*sin_t*vz-2.0*L2*cos_t*vx*sin_t*vz+wx2);

    start = (MapleGenVar1 - sqrt(arg))*MapleGenVar2/cc;
    stop = (MapleGenVar1 + sqrt(arg))*MapleGenVar2/cc;

    Adjust_Start_Stop_Time (Xee, start, stop);
    //cout.setf(LEFT);
    //cout<<"start stop = "<<PR2<<setw(W2)<<ze0/cc*0.5<<setw(W2)<<(start+stop)/2<<stop-start<<endl;

    if (RANDOMIZE_START)
	{
		Randomize_Start_Time(start,stop);
	}

    return(0);
}//Get_Start_Time

/******************************************************************************************/
int Randomize_Start_Time(prec &start,prec &stop)
/* Shifts start and stop integration time by a fraction of the time step size*/
{
	prec dt;    
    prec shift_amount;
	dt = (stop-start)/T_STEPS;
	if (Rand_Start_Init_Value ==0)
	{
		Rand_Start_Init_Value = Rand_Start_Init_Value_0;
	}
	
	shift_amount = (ran2(&Rand_Start_Init_Value) - 0.5)*RANDOMIZE_START_SCALE;
	shift_amount = dt*shift_amount;
	start = start + shift_amount;
	stop = stop + shift_amount; 

	return(0);
}

/******************************************************************************************/
int Adjust_Start_Stop_Time (prec Xee[], prec& start, prec& stop)
/* This procedure takes the calculated start and stop time of the time integral from
Get_Start_Stop_Time, and checks the contrast in the photon density between the start of
the integration and the middle of the integration. If the contrast is less than 10^4,
the time window of the integraion is expanded to accomodate. */
{
    prec a1,a2, r, t0,dt, alpha,dt2;

    t0 = 0.5*(start+stop);
    a1 = Photon_Density(Xee, t0);
    a2 = Photon_Density(Xee, start);

    if (a1 < 1.0e-150) return(0);		
		
    r = a2/a1;

    if (r < 0.0001) return(0);

    dt = stop-start;
    dt2=dt*dt;

    alpha = -log(r)/dt2;

    if (alpha<0.0) return(1);

    dt = sqrt(9.21/alpha);

    start = t0-dt*0.66;
    stop = t0+dt*0.66;

    return(0);

}//Adjust Start_Stop_Time

/******************************************************************************************/
int Find_W_Range(prec Xee[], prec phi, prec theta, prec theta_inc, prec &w_start, prec & dw)
/*
    Given the electron beam direction Xee, and the observation angle phi, and theta, and the laser
    incident angle, theta_inc, this procedure finds and returns the integration start frequency,
    w_start, and the integration step size, dw.
*/
{
    prec h;
    prec w_end;

    h = H(Xee, theta_inc, 0.0, phi, theta);
    w_start = h - INTEGRAL_RANGE_SIGMA_W*h*BANDWIDTH;
    w_end = h + INTEGRAL_RANGE_SIGMA_W*h*BANDWIDTH;

    if (w_start < 0.0) w_start = 0.0;

    dw=(w_end-w_start)/prec(N_Ws-1);

    return(0);

}//Find_W_Range


/******************************************************************************************/

int Get_Chirp_Spectrum (prec phi_0, prec theta_0, prec start_thetax, prec stop_thetax, prec d_thetax, prec start_thetay, prec stop_thetay, prec d_thetay, prec Ws_Min, prec dW, int numpar, prec int_factor)
{
    prec start_time;
    prec stop_time;
    prec mid_time;
    prec dt;
    prec time_zero;
    prec Xee[6];
    prec ing;
    prec phi,theta;
    prec time;
    int i,j,k;
    prec W;
    int t_bin;
    int t_steps;
    prec CS;
    prec wstart[MAX_TIME_BIN];
    prec wstop[MAX_TIME_BIN];


    for (i=0; i<numpar; i++)
    {
        phi=phi_0;
        theta=theta_0;
        COPYARRAY(X[i],Xee,6);

        if(Get_Start_Stop_Time(Xee,start_time,stop_time))
        {
           ing = 0.0;
        }
        else
        {
           mid_time = (start_time+stop_time)*0.5;
           Alter_Observation(Xee, phi, theta, mid_time);
           dt = (stop_time-start_time)/prec(T_STEPS);
           t_steps = int( (stop_time-start_time)/dt + 0.5);
           CS = int_factor*CrossSection(Xee, phi, theta, THETA_0, 0);
           if (PLANE_WAVE)
            Find_W_Range_Chirp(Xee, THETA_0, 0.0, phi, theta, BANDWIDTH, start_time, dt, t_steps, wstart, wstop);
           else
            Find_W_Range_Chirp_Kp(Xee, start_thetax, stop_thetax, start_thetay, stop_thetay, phi, theta, BANDWIDTH, start_time, dt, t_steps, wstart, wstop);
           for (j=0;j<N_Ws;j++)
           {
                W=Ws_Min+prec(j)*dW;
                time = start_time;
                for (k=0;k<t_steps;k++)
                {
                   if ( (W >= wstart[k]) && (W<=wstop[k]) )
                   {
                    ing = Integral_ngnw(Xee,time,dt,THETA_0,0.0,phi,theta,W,BANDWIDTH, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay);
                    time_zero=time+0.5*dt;
                    t_bin=FIND_TIME_BIN(time_zero,Xee,phi,theta);
                    if (PLANE_WAVE)
                        N_X[t_bin][j]=N_X[t_bin][j] + CS*ing;
                    else
                    {
                        //N_X[t_bin][j]=N_X[t_bin][j] + ing; Use this if CS is calculated in dN_dw_domega
                        N_X[t_bin][j]=N_X[t_bin][j] + CS*ing;
                    }
                   }
                   time = time+dt;
                }//for k
            }//for j
        }//else

       if((i%NUMSKIP_CHIRP)==0)
       {
          if(FIRST_TIME)
		  {
			cout<<i<<endl;
			cout<<"ing = "<<ing<<endl;
		  }
       }
    }//for i

  if(FIRST_TIME)
  {
	FIRST_TIME = 0;
	for (j=0;j<N_Ws;j++)
		{
			W=Ws_Min+prec(j)*dW;
			E_X[j]=W;
		}//for
  }

    return(0);

}//Get_Chrip Spectrum


/******************************************************************************************/
int Get_Spectrum (prec phi_0, prec theta_0, prec start_thetax, prec stop_thetax, prec d_thetax, prec start_thetay, prec stop_thetay, prec d_thetay, prec Ws_Min, prec dW, int numpar, prec int_factor)
/* Calculates the spectrum at a given observation direction by summing over the spectrum produced by
all the macro-particles in the simulation. If PLANE_WAVE = 1, assumes laser is a plane wave (i.e., does not sum over
initial k_perp's). If PLANE_WAVE=0, sums over initial k-vector states. */
{
    prec mid_time;
    prec time_zero;
    prec Xee[6];
    prec ing;
    prec phi,theta; //observation angles
    int i,j;
    prec W;
    int t_bin;
    
	if (FIRST_TIME)
		Init_SpotSize_Vars();
    for (i=0; i<numpar; i++)
    {
       COPYARRAY(X[i],Xee,6);
       phi = phi_0;
       theta=theta_0;

       Get_Time_Integraion(Xee, ing, mid_time, phi, theta, i);

       //time_zero = Xee[2]/cc;
       time_zero = mid_time;

       t_bin=FIND_TIME_BIN(time_zero,Xee,phi,theta);


       if(PLANE_WAVE)
       {
        for (j=0;j<N_Ws;j++)
        {
           W=Ws_Min+prec(j)*dW;
           N_X[t_bin][j]=N_X[t_bin][j] + int_factor*CrossSection(Xee, phi, theta, THETA_0, 0)*OMEGA_FACTOR(Xee, THETA_0, 0, phi, theta, W, BANDWIDTH)*ing;
        }//for
       }//if
       else
       {
        for (j=0;j<N_Ws;j++)
        {
           W=Ws_Min+prec(j)*dW;
           N_X[t_bin][j]=N_X[t_bin][j] + int_factor*Integral_dN_x_y(Xee,phi,theta,W,BANDWIDTH, ing, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay,0.,0.);
        }//for
       }//else

       if((i%NUMSKIP)==0)
       {
          if(FIRST_TIME)
          {
            cout<<i<<endl;
            cout<<"ing = "<<ing<<endl;
          }
       }
    }//for

    if(FIRST_TIME)
    {
        FIRST_TIME = 0;
        for (j=0;j<N_Ws;j++)
        {
            W=Ws_Min+prec(j)*dW;
            E_X[j]=W;
        }//for
    }

    return(0);

}//Get_Spectrum

/**************************************************************************************************/

int Get_Angle_Spectrum_All (prec angle_start, prec angle_stop, prec angle_fixed, int angle_num, prec start_thetax, prec stop_thetax, prec d_thetax, prec start_thetay, prec stop_thetay, prec d_thetay, prec Ws_Min, prec dW, int numpar, prec int_factor)
/*
Calls Get_Angle_Spectrum_One for each desired observation point specifed by the global variables
ANGLE_SPECTRUM, ANGLE_SPECTRUM_START, ANGLE_SPECTRUM_STOP, ANGLE_SPECTRUM_FIXED, and ANGLE_SPECTRUM_NUM.
*/
{
    prec d_angle;
    prec angle;
    prec phi, theta;
    prec W;
    int i;

    if(ANGLE_SPECTRUM==1)
    {
        phi = angle_fixed;
        angle = angle_start;
        d_angle = (angle_stop - angle_start)/prec(angle_num-1);
    }
    else if(ANGLE_SPECTRUM==2)
    {
        theta = angle_fixed;
        angle = angle_start;
        d_angle = (angle_stop - angle_start)/prec(angle_num);
    }
    else return(1);


    for (i=0; i<angle_num; i++)
    {

        T_X[i] = angle;

        if(ANGLE_SPECTRUM==1)
            Get_Angle_Spectrum_One (i, phi, angle, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay, Ws_Min, dW, numpar, int_factor);
        else
            Get_Angle_Spectrum_One (i, angle, theta, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay, Ws_Min, dW, numpar, int_factor);
        if(ANGLE_SPECTRUM==1)
        {
          cout<<"phi = "<<phi<<"   theta = "<<angle<<endl;
        }
        else
        {
          cout<<"phi = "<<angle<<"   theta = "<<theta<<endl;
        }
        angle = angle + d_angle;
    }//for


    for (i=0;i<N_Ws;i++)
    {
       W=Ws_Min+prec(i)*dW;
       E_X[i]=W;
    }//for

    return(0);

}//Get_Angle_Spectrum_All

/******************************************************************************************/
int Get_Angle_Spectrum_One (int angle_bin, prec phi_0, prec theta_0, prec start_thetax, prec stop_thetax, prec d_thetax, prec start_thetay, prec stop_thetay, prec d_thetay, prec Ws_Min, prec dW, int numpar, prec int_factor)
/* Calculates the spectrum at a given observation direction (phi_0, theta_0) by summing over the spectrum produced by
all the macro-particles in the simulation. If PLANE_WAVE = 1, assumes laser is a plane wave (i.e., does not sum over
initial k_perp's). If PLANE_WAVE=0, sums over initial k-vector states specfied by start_theta, stop_theta. This procedure
differs from Get_Spectrum in that it does stores in NX[angle][Ws] the spectrum at each observation
point rather than the spectrum at each observaton time.*/
{
    prec mid_time;
    prec Xee[6];
    prec ing;
    prec phi,theta;
    int i,j;
    prec W;
	
	//The following are used only if CHIRP_LASER = 1
		prec start_time, stop_time;
		int t_steps;
		int k;
		prec CS;
		prec time,dt;
		prec wstart[MAX_TIME_BIN];
		prec wstop[MAX_TIME_BIN];	

    for (i=0; i<numpar; i++)
    {
       COPYARRAY(X[i],Xee,6);
       phi = phi_0;
       theta=theta_0;

       if(!CHIRP_LASER)
	   {
			Get_Time_Integraion(Xee, ing, mid_time, phi, theta, i);

			if(PLANE_WAVE)
			{
				for (j=0;j<N_Ws;j++)
				{
					W=Ws_Min+prec(j)*dW;
					N_X[angle_bin][j]=N_X[angle_bin][j] + int_factor*CrossSection(Xee, phi, theta, THETA_0, 0)*OMEGA_FACTOR(Xee, THETA_0, 0, phi, theta, W, BANDWIDTH)*ing;
				}//for
			}//if
			else
			{
				for (j=0;j<N_Ws;j++)
				{
					W=Ws_Min+prec(j)*dW;
					N_X[angle_bin][j]=N_X[angle_bin][j] + int_factor*Integral_dN_x_y(Xee,phi,theta,W,BANDWIDTH, ing, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay,0.,0.);
				}//for
			}//else

			if((i%NUMSKIP)==0)
			{
				if(FIRST_TIME)
				{
					cout<<i<<endl;
					cout<<"ing = "<<ing<<endl;
				}
			}//if
       }//if(!CHIRP_LASER)
       else
	   {
			if(Get_Start_Stop_Time(Xee,start_time,stop_time))
			{
				ing = 0.0;
			}
			else
			{
				mid_time = (start_time+stop_time)*0.5;
				Alter_Observation(Xee, phi, theta, mid_time);
				dt = (stop_time-start_time)/prec(T_STEPS);
				t_steps = int( (stop_time-start_time)/dt + 0.5);
				CS = int_factor*CrossSection(Xee, phi, theta, THETA_0, 0);
				if (PLANE_WAVE)
					Find_W_Range_Chirp(Xee, THETA_0, 0.0, phi, theta, BANDWIDTH, start_time, dt, t_steps, wstart, wstop);
				else
					Find_W_Range_Chirp_Kp(Xee, start_thetax, stop_thetax, start_thetay, stop_thetay, phi, theta, BANDWIDTH, start_time, dt, t_steps, wstart, wstop);
				for (j=0;j<N_Ws;j++)
				{
					W=Ws_Min+prec(j)*dW;
					time = start_time;
					for (k=0;k<t_steps;k++)
					{
						if ( (W >= wstart[k]) && (W<=wstop[k]) )
						{
							ing = Integral_ngnw(Xee,time,dt,THETA_0,0.0,phi,theta,W,BANDWIDTH, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay);						
							if (PLANE_WAVE)
								N_X[angle_bin][j]=N_X[angle_bin][j] + CS*ing;
							else
							{
								//N_X[angle_bin][j]=N_X[angle_bin][j] + ing; Use this if CS is calculated in dN_dw_domega
								N_X[angle_bin][j]=N_X[angle_bin][j] + CS*ing;
							}
						}
						time = time+dt;
					}//for k
				}//for j
			}//else

			if((i%NUMSKIP_CHIRP)==0)
			{
				if(FIRST_TIME)
				{
					cout<<i<<endl;
					cout<<"ing = "<<ing<<endl;
				}
			}
		}//else //if(!CHIRP_LASER)
    }//for

    if(FIRST_TIME)
        FIRST_TIME = 0;

    return(0);

}//Get_Angle_Spectrum


/*********************************************************************************/

int Write_Angle_Spectrum (prec angle_start, prec angle_stop, prec angle_fixed, int angle_num)
/* Writes the calculated x-ray spectrum to a file */
{
    int i,j;
    prec W,dw;
    prec angle, d_angle, theta, phi;

    ofstream outfile3;
    ofstream outfile4;

    outfile3.open("Graph_Angle_Spectrum.txt");
    outfile3.setf(LEFT);

    outfile4.open("Graph_Angle_Spectrum2.txt");
    outfile4.setf(LEFT);


    if(ANGLE_SPECTRUM==1)
    {
        phi = angle_fixed;
        angle = angle_start;
        d_angle = (angle_stop - angle_start)/prec(angle_num-1);
    }
    else if(ANGLE_SPECTRUM==2)
    {
        theta = angle_fixed;
        angle = angle_start;
        d_angle = (angle_stop - angle_start)/prec(angle_num);
    }

    dw=E_X[1]-E_X[0];

    if (ANGLE_SPECTRUM==1)
    {
        outfile3<<setw(W2)<<"Theta(mrad)"<<setw(W2)<<"Energy(keV)"<<"keV/mrad^2/eV"<<"  "<<angle_num<<"  "<<N_Ws<<endl<<endl;
        outfile4<<setw(W2)<<"Theta(mrad)"<<setw(W2)<<"Energy(keV)"<<"photons/mrad^2/eV"<<"  "<<angle_num<<"  "<<N_Ws<<endl<<endl;
    }
    else
    {
        outfile3<<setw(W2)<<"Phi(degrees)"<<setw(W2)<<"Energy(keV)"<<"keV/mrad^2/eV"<<"  "<<angle_num<<"  "<<N_Ws<<endl<<endl;
        outfile4<<setw(W2)<<"Phi(degrees)"<<setw(W2)<<"Energy(keV)"<<"photons/mrad^2/eV"<<"  "<<angle_num<<"  "<<N_Ws<<endl<<endl;
    }


    for (j=0; j<angle_num; j++)
    {
        if (ANGLE_SPECTRUM==1)
        {
            angle = T_X[j]*1000.0;
        }
        else
        {
            angle = T_X[j]/PI*180.;
        }
        for (i=0; i<N_Ws; i++)
        {
            //  NvTime[j]=NvTime[j]+N_X[j][i]*dw/1.0e6;  //photons/mrad^2
            //  IvTime[j]=IvTime[j]+(E_SCALE)*E_X[i]*N_X[j][i]*dw/1.0e6; //keV/mrad^2
                W=E_SCALE*E_X[i];
                outfile3<<setw(W2)<<PR2<<angle<<setw(W2)<<W<<N_X[j][i]*W/E_SCALE/1.0e9<<endl;
                outfile4<<setw(W2)<<PR2<<angle<<setw(W2)<<W<<N_X[j][i]/E_SCALE/1.0e9<<endl;
        }//for i
    }//for j

    outfile3.close();
    outfile4.close();

    return(0);

}//Write_Angle_Spectrum

/******************************************************************************************/
int Alter_Observation(prec Xee[], prec &phi, prec &theta, prec time)
/* This procedure alters the lab observation direction (phi, theta) to account for the
location of the radiating particle at time given by "time". This becomes significant
when the transverse location is not extremely small compared to the observation distance
from the interaction specifed in the global varibale "OBS_DISTANCE".
*/
{
    prec x0[3];
    prec temp;

    if(theta<0.)
    {
        theta = -theta;
        phi = phi+PI;
    }

    x0[0]=OBS_DISTANCE*tan(theta)*cos(phi);
    x0[1]=OBS_DISTANCE*tan(theta)*sin(phi);
    x0[2]=OBS_DISTANCE;

    x0[0]=x0[0]- Current_Position(Xee, 0, time);
    x0[1]=x0[1] - Current_Position(Xee, 1, time);
    x0[2]=x0[2] - Current_Position(Xee, 2, time);

    phi = atan2(x0[1],x0[0]);


    temp=sqrt(x0[0]*x0[0]+x0[1]*x0[1]);
    theta = atan2(temp,x0[2]);

    return (0);
}//Alter_Observation

/******************************************************************************************/
prec Integral_dN_x_y (prec Xe[], prec phi, prec theta, prec ws, prec dw, prec i_ng, prec t_startx, prec t_stopx, prec dtx, prec t_starty, prec t_stopy, prec dty, prec t_start, prec dt)
/*
Integrates the scattered photon density per solid angle per unit frequency over all incident laser photon x directions and y directions
for a given electron coordinates Xe[x,y,z,px,pz,py], scattered frequency ws,
laser bandwidth dw, and obsvervation angles phi and theta. t_startx and t_stopx define the integral range in x, while
t_starty and t_stopy define the integral range in y.*/
{
   int i,stepsx;
   stepsx = int( (t_stopx - t_startx)/dtx + 0.5);
   prec ans,d_ans, tlx;
   ans=0.;
   d_ans = 0.;
   tlx=t_startx;
   for(i=0; i<stepsx; i++)
   {
        d_ans = Integral_dN_y(Xe, tlx, phi, theta, ws, dw, i_ng, t_starty, t_stopy, dty,t_start,dt);
        rk4_dNxy(Xe,ans,d_ans,dtx,tlx,phi, theta, ws, dw, i_ng,t_starty,t_stopy,dty,t_start,dt);
        tlx = tlx + dtx;
   }


   return(ans);

}//Integral_dN_x_y

/******************************************************************************************/

prec Integral_dN_y (prec Xe[], prec tlx, prec phi, prec theta, prec ws, prec dw, prec i_ng, prec t_start, prec t_stop, prec dt, prec time_start, prec d_time)
/* Integrates the scattered photon density per solid angle per unit frequency over all incident laser photon y directions
for a given for a given incident laser direction in x (tlx), a given electron coordinates Xe[x,y,z,px,pz,py], scattered frequency ws,
laser bandwidth dw, and obsvervation angles phi and theta. t_start and t_stop define the integral range.*/
{
   int i,steps;
   steps = int( (t_stop - t_start)/dt + 0.5);
   prec ans,d_ans, tly;
   ans=0.;
   d_ans = 0.;
   tly=t_start;
   for(i=0; i<steps; i++)
   {
        d_ans = dN_dw_domega (Xe, tlx, tly, phi, theta, ws, dw, i_ng, time_start, d_time);
        rk4_dNy(Xe,ans,d_ans,dt,tly,tlx,phi, theta, ws, dw, i_ng, time_start, d_time);
        tly = tly + dt;
   }

   return(ans);


}//Integral_dN_y

/******************************************************************************************/

void rk4_dNxy (prec Xe[], prec& Y, prec dY, prec h, prec time, prec phi, prec theta, prec ws, prec dw, prec i_ng, prec t_starty, prec t_stopy, prec dty, prec t_start, prec dt)
/* This procedure advances the value Y by one time step using the runge-kutta 4 method. Xe(x,y,z,px,py,pz)
contains the initial coordinates of the electron, which in required in the DERIVATIVE procedure.*/
{
   prec hh,h6;
   prec t = time;

   hh = 0.5*h;
   h6 = h/6.0;

   prec K1;
   prec K2;
   prec K3;
   prec K4;

   K1 = dY;
   t = time + h*0.5;
   K2 = Integral_dN_y(Xe, t, phi, theta, ws, dw, i_ng, t_starty, t_stopy, dty,t_start,dt);
   K3 = K2;
   t = time + h;
   K4 = Integral_dN_y(Xe, t, phi, theta, ws, dw, i_ng, t_starty, t_stopy, dty,t_start,dt);

   Y = Y + h6*(K1 +2*(K2 + K3) +K4);
}


/******************************************************************************************/

void rk4_dNy (prec Xe[], prec& Y, prec dY, prec h, prec time, prec tlx, prec phi, prec theta, prec ws, prec dw, prec i_ng, prec t_start, prec dt)
/* This procedure advances the value Y by one time step using the runge-kutta 4 method. Xe(x,y,z,px,py,pz)
contains the initial coordinates of the electron, which in required in the DERIVATIVE procedure.*/
{
   prec hh,h6;
   prec t = time;

   hh = 0.5*h;
   h6 = h/6.0;

   prec K1;
   prec K2;
   prec K3;
   prec K4;

   K1 = dY;
   t = time + h*0.5;
   K2 = dN_dw_domega (Xe, tlx, t, phi, theta, ws, dw, i_ng,t_start,dt);
   K3 = K2;
   t = time + h;
   K4 = dN_dw_domega (Xe, tlx, t, phi, theta, ws, dw, i_ng,t_start,dt);

   Y = Y + h6*(K1 +2*(K2 + K3) +K4);
}//rk4_dNy


/******************************************************************************************/

prec dN_dw_domega (prec Xe[], prec tlx, prec tly, prec phi, prec theta, prec ws, prec dw, prec i_ng, prec t_start, prec dt)
/* Calculates the scattered photon density per solid angle per unit frequency for a given incident
laser direction (tlx, tly), a given electron coordinates Xe[x,y,z,px,pz,py], scattered frequency ws,
laser bandwidth dw, and obsvervation angles phi and theta. */
{
    prec ans1,ans2,ans3;
    prec ans, dn_g,n_g;

    dn_g =0.;
    n_g =0.;


    ans3 = fk_perp(tlx, tly);
    if (CHIRP_LASER==0)
    {
        ans1 = CrossSection(Xe, phi, theta, tlx, tly);
        ans2 = OMEGA_FACTOR(Xe, tlx, tly, phi, theta, ws, dw);
        ans = ans1*ans2*ans3*i_ng;
    }//if
    else
    {
        //ans1 = CrossSection(Xe, phi, theta, tlx, tly);
        /* To save computation time, CrossSection is only calculated once for each
           macro_particle */
        ans1 = 1.0;
        DERIVATIVE_CHIRP(Xe,dn_g,t_start,tlx,tly,phi,theta,ws,dw);
        rk4_chirp(Xe,n_g,dn_g,dt,t_start,tlx,tly,phi,theta,ws,dw);
        ans = ans1*ans3*n_g;
    }//else

//	vel_factor = VEL_FACTOR(Xe,tlx,tly);
//	ans = ans*vel_factor;

    return(ans);

}//dN_dw_domega

/******************************************************************************************/
prec Integral_ngnw (prec Xe[],prec t_start, prec dt, prec tlx, prec tly, prec phi, prec theta, prec ws, prec dw0, prec start_thetax, prec stop_thetax, prec d_thetax, prec start_thetay, prec stop_thetay, prec d_thetay)
/* Integrates the overlap integral of the laser photon density multiplied by the initial frequency
distribution function of a give scattered frequency ws. with an electron with initial
position specified in Xe[x,y,z,px,py,pz] from time t_start to t_start + dt. tlx,y specifies the
initial photon direction, phi and theta the observation direction*/
{

   prec vel_factor;

   prec n_g=0.;
   prec dn_g=0.;

   if (PLANE_WAVE)
   {
    DERIVATIVE_CHIRP(Xe,dn_g,t_start,tlx,tly,phi,theta,ws,dw0);
    rk4_chirp(Xe,n_g,dn_g,dt,t_start,tlx,tly,phi,theta,ws,dw0);
   }//if
   else
   {
      n_g = Integral_dN_x_y(Xe,phi,theta,ws,dw0, 1.0, start_thetax, stop_thetax, d_thetax, start_thetay, stop_thetay, d_thetay, t_start, dt);
   }

   vel_factor = VEL_FACTOR(Xe,THETA_0,0.0);
   n_g = n_g*vel_factor;

   return(n_g);

}//Integral_ng

/******************************************************************************************/

void rk4_chirp (prec Xe[], prec& Y, prec dY, prec h, prec time, prec tlx, prec tly, prec phi, prec theta, prec ws, prec dw0)
/* This procedure advances the value Y by one time step using the runge-kutta 4 method. Xe(x,y,z,px,py,pz)
contains the initial coordinates of the electron, which in required in the DERIVATIVE procedure.*/
{
   prec hh,h6;
   prec t = time;

   hh = 0.5*h;
   h6 = h/6.0;

   prec K1;
   prec K2;
   prec K3;
   prec K4;


   K1 = dY;
   t = time + h*0.5;
   DERIVATIVE_CHIRP(Xe,K2,t,tlx,tly,phi,theta,ws,dw0);
   K3=K2;
   t = time + h;
   DERIVATIVE_CHIRP(Xe,K4,t,tlx,tly,phi,theta,ws,dw0);
   Y = Y + h6*(K1 +2.0*(K2 + K3) +K4);
}//rk4_chirp

/******************************************************************************************/
prec  Integral_ng(prec Xe[],prec t_start, prec t_stop)
/* Integrates the overlap integral of the laser photon density with an electron with initial
position specified in Xe[x,y,z,px,py,pz] from time t=T_START to t=T_STOP*/
{

   int i;
   int t_steps;
   prec time=0.;
   prec vel_factor;
   time = t_start;
   prec dt;

  // ofstream out_test;
  // out_test.open("test.txt");
  // out_test.setf(LEFT);


   dt = (t_stop-t_start)/prec(T_STEPS);
   t_steps = int( (t_stop - t_start)/dt + 0.5);

   prec n_g=0.;
   prec dn_g=0.;

   for(i=0; i<t_steps; i++)
   {
        DERIVATIVE(Xe,dn_g,time);
        rk4(Xe,n_g,dn_g,dt,time);
        time = time + dt;
  //      out_test<<setw(20)<<time<<Photon_Density(Xe, time)<<endl;
   }

   vel_factor = VEL_FACTOR(Xe,THETA_0,0.0);
   n_g = n_g*vel_factor;

   //out_test.close();

   return(n_g);

}//Integral_ng

/******************************************************************************************/

void rk4 (prec Xe[], prec& Y, prec dY, prec h, prec time)
/* This procedure advances the value Y by one time step using the runge-kutta 4 method. Xe(x,y,z,px,py,pz)
contains the initial coordinates of the electron, which in required in the DERIVATIVE procedure.*/
{
   prec hh,h6;
   prec t = time;

   hh = 0.5*h;
   h6 = h/6.0;

   prec K1;
   prec K2;
   prec K3;
   prec K4;


   K1 = dY;
   t = time + h*0.5;
   DERIVATIVE(Xe,K2,t);
   K3=K2;
   t = time + h;
   DERIVATIVE(Xe,K4,t);
   Y = Y + h6*(K1 +2.0*(K2 + K3) +K4);
}

/***********************************************************************************/

void DERIVATIVE (prec Xe[], prec& dY,prec t)
{

    dY = Photon_Density(Xe, t);

} //DERIVATIVE

/***********************************************************************************/

void DERIVATIVE_CHIRP (prec Xe[], prec& dY,prec t, prec tlx, prec tly, prec phi, prec theta, prec ws, prec dw0)
{
    prec photon_density;
	
	photon_density = Photon_Density(Xe, t);
    dY = photon_density*OMEGA_FACTOR_CHRIP (Xe, tlx, tly, phi, theta, ws,dw0, t,photon_density);
//  dY = Photon_Density(Xe, t)*OMEGA_FACTOR (Xe, tlx, tly, phi, theta, ws,dw0);

} //DERIVATIVE

/***********************************************************************************/

prec GAMMA (prec x[])
/* Given the coordinate vector x[x,y,z,px,py,pz] this procedure calculates the relativistic
   gamma factor. */
{
    int i;
    prec value = 0.;


    for (i=3; i<6; i++)
        value = value + x[i]*x[i];

    value = sqrt(value + 1);

    return (value);

} //GAMMA

/***********************************************************************************/

prec BETA (prec x[])
/* Given the coordinate vector x[x,y,z,px,py,pz] this procedure calculates the relativistic
   beta factor. */
{
    prec value = 0.;

    value = GAMMA(x);
    value = sqrt(1.0-1.0/(value*value));

    return (value);

} //BETA

/***********************************************************************************/

prec XP (prec x[])
/* Given the coordinate vector x[x,y,z,px,py,pz] this procedure calculates and return dx/dz */
{

    return (x[3]/x[5]);

} //XP

/***********************************************************************************/

prec YP (prec x[])
/* Given the coordinate vector x[x,y,z,px,py,pz] this procedure calculates and return dy/dz */
{

    return (x[4]/x[5]);

} //XP

/***********************************************************************************/

prec CROSS (prec x[], prec y[], int k)
// Returns the k component of the cross product of  x X y.
{
    // int i,j;
    prec value;

    value = 0.;
    if(k==0)
        value=x[1]*y[2] - x[2]*y[1];
    else if(k==1)
        value=x[2]*y[0] - x[0]*y[2];
    else
        value=x[0]*y[1] - x[1]*y[0];

    return(value);

}//CROSS

/***********************************************************************************/
int Find_W_Range_Chirp_Kp(prec Xe[], prec tlx_start, prec tlx_stop, prec tly_start, prec tly_stop, prec phi, prec theta, prec dw0, prec start_t, prec dt, int t_steps, prec w_start[], prec w_stop[])
{
    prec w_start_temp[MAX_TIME_BIN];
    prec w_stop_temp[MAX_TIME_BIN];

    Find_W_Range_Chirp(Xe, tlx_start, tly_start, phi, theta, dw0, start_t, dt, t_steps, w_start_temp, w_stop_temp);
    COPYARRAY(w_start_temp,w_start, t_steps);
    COPYARRAY(w_stop_temp,w_stop, t_steps);
    Find_W_Range_Chirp(Xe, tlx_stop, tly_start, phi, theta, dw0, start_t, dt, t_steps, w_start_temp, w_stop_temp);
    COPYARRAY_MIN(w_start_temp,w_start, t_steps);
    COPYARRAY_MAX(w_stop_temp,w_stop, t_steps);
    Find_W_Range_Chirp(Xe, tlx_start, tly_stop, phi, theta, dw0, start_t, dt, t_steps, w_start_temp, w_stop_temp);
    COPYARRAY_MIN(w_start_temp,w_start, t_steps);
    COPYARRAY_MAX(w_stop_temp,w_stop, t_steps);
    Find_W_Range_Chirp(Xe, tlx_stop, tly_stop, phi, theta, dw0, start_t, dt, t_steps, w_start_temp, w_stop_temp);
    COPYARRAY_MIN(w_start_temp,w_start, t_steps);
    COPYARRAY_MAX(w_stop_temp,w_stop, t_steps);

    return(0);

}//Find_W_Range_Kp
/***********************************************************************************/

int Find_W_Range_Chirp(prec Xe[], prec tlx, prec tly, prec phi, prec theta, prec dw0, prec start_t, prec dt, int t_steps, prec w_start[], prec w_stop[])
{
    prec h;
    prec Xl[3];
    prec z;
    prec f;
    prec dw;
    prec cw;
    prec DW;
    int i;
    prec time;

    cw = CHIRP_WIDTH*CHIRP_WIDTH;
    if(cw<1.0)
    {
        dw = dw0;
        DW = 0.0;
    }
    else
    {
        //dw = dw0/sqrt(cw);
        dw=dw0;
        //DW= dw0*sqrt(1.0- 1.0/(cw));
        DW= dw0*sqrt(cw)*sqrt(1.0- 1.0/(cw));
    }

    
	if (!NONLINEAR_ON)
		h = H(Xe,tlx,tly,phi,theta);

    time = start_t+0.5*dt;
    for (i=0; i<t_steps; i++)
    {
        LaserFrame(Xe,Xl,time);
        z = Xl[2] + cc*time;
		if (NONLINEAR_ON)
		{
			h = H_Time_Dependent(Xe,tlx,tly,phi,theta,time,Photon_Density(Xe,time));
		}

        if(CHIRP_WIDTH<0.)
        f=1.0 + z/PULSE_WIDTH*DW;
        else
        f=1.0 - z/PULSE_WIDTH*DW;
        f=f*h;
        w_start[i] = f - INTEGRAL_RANGE_SIGMA_W*h*dw;
        w_stop[i] = f + INTEGRAL_RANGE_SIGMA_W*h*dw;
        time = time+dt;
    }
    return(0);
}


/***********************************************************************************/

prec OMEGA_FACTOR_CHRIP (prec Xe[], prec tlx, prec tly, prec phi, prec theta, prec ws,prec dw0, prec t, prec photon_density)
/* Returns the value of the scale factor of the intnensity scattered at frequency ws
at observation angle phi and theta assuming initial photon direction tlx and tly and
initial electron beam coordinates Xe[x,y,z,px,py,pz]. dw is the total bandwidth of the laser
pulse, and CHIRP_WIDTH (global) is the ratio of the actual pulse length to the bandwidth limited
pulse length. Note, |chirp_factor| >= 1.  chirp_factor<0 => negative chrip. The global variable
PULSE_WIDTH is the 1/e^2 intensity pulse width (non-bandwidth limited of the laser pulse. */
{
    prec h, ans;
    prec Xl[3];
    LaserFrame(Xe,Xl,t);
    prec z;
    prec f;
    prec dw;
    prec cw;
    prec DW;

    cw = CHIRP_WIDTH*CHIRP_WIDTH;
    if(cw<1.0)
    {
        dw = dw0;
        DW = 0.0;
    }
    else
    {
        //dw = dw0/sqrt(cw);
        dw=dw0;
        //DW= dw0*sqrt(1.0- 1.0/(cw));
        DW= dw0*sqrt(cw)*sqrt(1.0- 1.0/(cw));
    }

    z = Xl[2] + cc*t;

    if(CHIRP_WIDTH<0.)
        f=1.0 + z/PULSE_WIDTH*DW;
    else
        f=1.0 - z/PULSE_WIDTH*DW;
    
	if (NONLINEAR_ON)
		h = H_Time_Dependent(Xe,tlx,tly,phi,theta,t,photon_density);
	else		   
		h = H(Xe,tlx,tly,phi,theta);

    ans = sqrt(2.0/PI) * EXP(-2.0*(ws/h-f)*(ws/h-f)/(dw*dw) )/(h*dw)*FREQ_ATTEN_FACTOR(Xe,ws,phi,theta);

    return(ans);

}

/***********************************************************************************/

prec OMEGA_FACTOR (prec Xe[], prec tlx, prec tly, prec phi, prec theta, prec ws,prec dw)
/* Returns the value of the scale factor of the intnensity scattered at frequency ws
at observation angle phi and theta assuming initial photon direction tlx and tly and
initial electron beam coordinates Xe[x,y,z,px,py,pz]. */
{
    prec h, ans;
    h = H(Xe,tlx,tly,phi,theta);
    ans = sqrt(2.0/PI) * EXP(-2.0*(ws/h-1.0)*(ws/h-1.0)/(dw*dw) )/(h*dw)*FREQ_ATTEN_FACTOR(Xe,ws,phi,theta);

    return(ans);

}

/***********************************************************************************/
prec FREQ_ATTEN_FACTOR (prec Xe[],prec ws, prec phi, prec theta)
/* Attenuates the scattered photon inensity at a given frequency and observation
direction (phi, theta) for a function specified in this procedure. */
{
    prec T =1.0;
    prec E;
    E=ws*E_SCALE;
    if(WINDOW_ON)
        T = WINDOW_FILTER(ws);
    if(WINDOW_FILE_ON)
        T= T*INTERPOLATE(E, WINDOW_FILE_X, WINDOW_FILE_Y, WINDOW_FILE_N);
    if (DETECTOR_FILE_ON)
        T = T*INTERPOLATE(E, DETECTOR_FILE_X, DETECTOR_FILE_Y, DETECTOR_FILE_N);

    if(FILTER_ON)
        T = T*FILTER(ws,FILTER_WIDTH,FILTER_CENTER);

    if (BRAGG>0)
        T = T*BRAGG_REFLECT(ws, BRAGG, BRAGG_SPACING, BRAGG_ANGLE, BRAGG_WIDTH,theta,phi);

    return(T);


}//FREQ_ATTEN_FACTOR
/***********************************************************************************/
prec BRAGG_REFLECT(prec ws, int xy, prec b_space, prec b_angle, prec b_width,prec theta_obs, prec phi_obs)
{
    prec w_bragg, l_bragg;
    prec tx,ty;
    prec w,dw,T;



        tx = theta_obs*cos(phi_obs)+b_angle;
        ty = theta_obs*sin(phi_obs);
        l_bragg = b_space*sin(tx)/cos(ty);
        w_bragg = WAVELENGTH/l_bragg;
        dw= 1.0/(b_width*w_bragg*0.5);
        w=dw*(ws-w_bragg);
        w=w*w;
        T = 1.0/(1.0+w);

        return(T);

}//BRAGG_REFLECT

/***********************************************************************************/
prec WINDOW_FILTER (prec ws)
/* Returns the attenuation factor for transmission through a window for a photon frequency
specified by ws. */
{
    prec E;
    prec ans1,ans2;

    E = ws*E_SCALE;
    ans1 = 24.2*pow(E/10.,-3.2) + 0.19*pow(E/60.,-0.4);
    ans2 = EXP(-ans1*2.51*1.658);

    return(ans2);

}//WINDOW_FILTER

/***********************************************************************************/

prec H(prec Xe[], prec tlx, prec tly, prec phi, prec theta)
/* Returns the ratio of the scattered photon energy to the incident photon energy given the
electron coordinates Xe[x,y,z,px,py,pz], the photon direction: tlx, tly, and the observation
angles phi and theta */
{
    prec beta,xpe,ype;
    prec tx, ty,ans;

    xpe=XP(Xe);
    ype=YP(Xe);
    beta=BETA(Xe);
    tx=tlx-xpe;
    ty=tly-ype;

    //ans=(1.0+beta*cos(tx)*cos(ty))/(1.0-beta*COStheta(Xe, phi,theta));
    ans=(1.0+beta*(cos(tx)*cos(ty) + sin(tly)*sin(ype)*(1.0-cos(tx))))/(1.0-beta*COStheta(Xe, phi,theta));

    return(ans);
}

/***********************************************************************************/

prec H_Time_Dependent(prec Xe[], prec tlx, prec tly, prec phi, prec theta, prec t, prec photon_density)
/* Returns the ratio of the scattered photon energy to the incident photon energy given the
electron coordinates Xe[x,y,z,px,py,pz], the photon direction: tlx, tly, and the observation
angles phi and theta at specified time t. This called when nonlinear effects are to be included in the
calculatoin of the scattered x-ray wavelenght. This routine takes into account pondermotive shifting of 
the scattered wavelength */
{
    prec beta,xpe,ype;
    prec tx, ty,ans,a02;

    xpe=XP(Xe);
    ype=YP(Xe);
    beta=BETA(Xe);
    tx=tlx-xpe;
    ty=tly-ype;
    
	//Calculate the localized normalized laser intensity
//	a02 = Photon_Density(Xe, t);
	
	//Calculate sqaure of peak instantaneous vector potential based on the photon density (assumes linear polarization)
	a02 = (photon_density/(PI*N_electron))*(2*LAMBDA_C*R_E*WAVELENGTH)/(LBASE*LBASE);

    //Calculate the ratio between incident and scattered photon energy 
    ans=(1.0+beta*(cos(tx)*cos(ty) + sin(tly)*sin(ype)*(1.0-cos(tx))))/(1.0-beta*COStheta(Xe, phi,theta));

	//Apply shift due to non-linear pondermotive effects (scaled based on time averaged squared vector potential)
	ans = ans/(1+a02*0.5);

    return(ans);
}


/***********************************************************************************/
prec fk_perp (prec theta_x, prec theta_y)
/* returns the value of the k perpendicular distribution fucntion for the photon
angles theta_x and theta_y corresponding to the perpendicular k vector.
*/
{

    prec ans;

    ans = 1.0/(2.0*PI*sigma_theta_L_x*sigma_theta_L_y)*EXP(-0.5*(theta_x - THETA_0)*(theta_x - THETA_0)/(sigma_theta_L_x*sigma_theta_L_x)
           -0.5*theta_y*theta_y/(sigma_theta_L_y*sigma_theta_L_y));

    return(ans);


}//fk_perp

/***********************************************************************************/

prec CrossSection(prec Xe[], prec phi, prec theta, prec tlx, prec tly)
/* Given the electron beam coordinates Xe, and the laser incident angles tlx, and tly, this
function returns the value of the lab frame differntial cross section at the lab observation
angles, phi, and theta.
*/
{
    prec ans, ax, ay, az, beta, gamma;
    prec cospsint, sinpsint, cost, bcos_1;
    prec ans1, ans2, ans3, ans4, ans5, ans6, ans7;
//  prec vel_factor; //note, this geometric factor is actually not part of the differential cross-section, but
                     //but is necessary to calculate the scattered photon rate.

    ax = alpha_x(PHI_P, Xe, tlx, tly);
    ay = alpha_y(PHI_P, Xe, tlx, tly);
    az = alpha_z(PHI_P, Xe, tlx, tly);

    beta = BETA(Xe);
    gamma = GAMMA(Xe);

    cospsint = COSphiSINtheta(Xe, phi, theta);
    sinpsint = SINphiSINtheta(Xe, phi, theta);
    cost     = COStheta(Xe, phi, theta);
    bcos_1 = (1-beta*cost)*(1-beta*cost);

//  vel_factor = VEL_FACTOR(Xe,tlx,tly);
/* Velocity factor is currently calculated during the time integral procedure. For the most
    part, this avoids redundancy. But if the program is run without the plane wave approximation
    there will be some minor, albeit, most likely insignificant differences.*/

    ans1 = ax*ax*(1 - cospsint*cospsint/(gamma*gamma*bcos_1));
    ans2 = ay*ay*(1- sinpsint*sinpsint/(gamma*gamma*bcos_1));
    ans3 = az*az*(1 - (cost-beta)*(cost-beta)/bcos_1);
    ans4 = -2.0*ax*ay*(cospsint*sinpsint)/(gamma*gamma*bcos_1);
    ans5 = -2.0*ax*az*(cost-beta)*cospsint/(gamma*bcos_1);
    ans6 = -2.0*ay*az*(cost-beta)*sinpsint/(gamma*bcos_1);
    ans7 = (1-beta*beta)/bcos_1;

    ans = (ans1+ans2+ans3+ans4+ans5+ans6)*ans7;//*vel_factor;

    return(ans);

}//CrossSection


/***********************************************************************************/
prec VEL_FACTOR(prec Xe[], prec tlx, prec tly)
/* Given the electron beam coordinates Xe, and the laser beam direction tlx and tly, this
procedures returns the veloctiy factor that determines the rate of scattered photons depending
on the differential cross section. This factor ranges fromm c(1+beta)*re^2 for a head on collison, to
c*(1-beta)*re^2 for a collinear collision.
*/
{
    prec gamma,betax,betay,betaz;
    prec kb;
    prec ans;

    gamma = GAMMA(Xe);
    gamma = 1./gamma;
    betax = Xe[3]*gamma;
    betay = Xe[4]*gamma;
    betaz = Xe[5]*gamma;

    kb = betax*sin(tlx) + (betay*sin(tly) + betaz*cos(tly))*cos(tlx);

    ans = cR_E2_N*(1.0 + kb);


return(ans);

}//VEL_FACTOR

/***********************************************************************************/
prec alpha_x (prec phi_p, prec Xe[], prec tlx, prec tly)
/* Given the rotation of the laser polarization vector about the laser k-vector (phi_p = 0
corresponds to the vector being in the x-z plane), the electron beam incident angles (tex, and tey), and
the laser incident angle (tlx, and tly), this function returns the x-component of the polarization vector in
the electron rest frame. */
/* Note, in the case below, the photon angles about the central laser direction are ignored. */
{
    prec alpha;
    prec tx;
    prec ty;
    prec beta;
    prec cosx;
    prec cosy;
    prec sinx;
    prec siny;

    tx = THETA_0-XP(Xe);
    ty = -YP(Xe);
    beta = BETA(Xe);

    cosx=cos(tx);
    cosy=cos(ty);
    sinx=sin(tx);
    siny=sin(ty);


    alpha = (cos(phi_p)*(cosx+beta*cosy) + beta*sin(phi_p)*sinx*siny)/(1.0+beta*cosx*cosy);

    return(alpha);

}//alpha_x

/***********************************************************************************/
prec alpha_y (prec phi_p, prec Xe[], prec tlx, prec tly)
/* Given the rotation of the laser polarization vector about the laser k-vector (phi_p = 0
corresponds to the vector being in the x-z plane), the electron beam incident angles (tex, and tey), and
the laser incident angle (tlx, and tly), this function returns the y-component of the polarization vector in
the electron rest frame. */
/* Note, in the case below, the photon angles about the central laser direction are ignored. */
{
    prec alpha;
    prec tx;
    prec ty;
    prec beta;
    prec cosx;
    prec cosy;
    prec sinx;
    prec siny;


    tx = THETA_0-XP(Xe);
    ty = -YP(Xe);
    beta = BETA(Xe);

    cosx=cos(tx);
    cosy=cos(ty);
    sinx=sin(tx);
    siny=sin(ty);

    alpha = (-cos(phi_p)*sinx*siny + sin(phi_p)*(cosy+beta*cosx) )/(1+beta*cosx*cosy);

    return(alpha);

}//alpha_y

/***********************************************************************************/

prec alpha_z (prec phi_p, prec Xe[], prec tlx, prec tly)
/* Given the rotation of the laser polarization vector about the laser k-vector (phi_p = 0
corresponds to the vector being in the x-z plane), the electron beam incident angles (tex, and tey), and
the laser incident angle (tlx, and tly), this function returns the z-component of the polarization vector in
the electron rest frame. */
/* Note, in the case below, the photon angles about the central laser direction are ignored. */
{
    prec alpha;
    prec tx;
    prec ty;
    prec beta;
    prec gamma;
    prec cosx;
    prec cosy;
    prec sinx;
    prec siny;


    tx = THETA_0-XP(Xe);
    ty = -YP(Xe);
    beta = BETA(Xe);
    gamma = GAMMA(Xe);

    cosx=cos(tx);
    cosy=cos(ty);
    sinx=sin(tx);
    siny=sin(ty);

    alpha = (-cos(phi_p)*sinx*cosy + sin(phi_p)*siny )/(gamma*(1+beta*cosx*cosy));

    return(alpha);

}//alpha_z

/***********************************************************************************/

prec COSphiSINtheta (prec Xe[], prec phi, prec theta)
/* Converts cos(phi)sin(theta), where theta and phi are defined in the coordinate system
where +z lies along the electron beam direction, to an expression in terms of the lab
frame phi and theta. Basically, alters the value of cos(phi)sin(theta) in the lab frame
cross section to account for the electron beam direction.
*/
{
    prec txe;
    prec tye;
    prec ans;

    prec sin_theta, cos_theta, sin_xe, cos_xe;
//  prec sin_ye, cos_ye;

    txe = XP(Xe);
    tye = YP(Xe);

    sin_theta = sin(theta);
    cos_theta = cos(theta);
    sin_xe = sin(txe);
//  sin_ye = sin(tye);
    cos_xe = cos(txe);
//  cos_ye = cos(tye);

//    ans = cos(phi)*sin_theta*cos_xe - sin(phi)*sin_theta*sin_xe*sin_ye- cos_theta*sin_xe*cos_ye;

    ans = cos(phi)*sin_theta*cos_xe - cos_theta*sin_xe;

    return(ans);

}//COSphiSINtheta

/***********************************************************************************/

prec SINphiSINtheta (prec Xe[], prec phi, prec theta)
/* Converts sin(phi)sin(theta), where theta and phi are defined in the coordinate system
where +z lies along the electron beam direction, to an expression in terms of the lab
frame phi and theta. Basically, alters the value of sin(phi)sin(theta) in the lab frame
cross section to account for the electron beam direction.
*/
{
    prec txe;
    prec tye;
    prec ans;

    txe = XP(Xe);
    tye = YP(Xe);

    //ans = sin(phi)*sin(theta)*cos(tye)-cos(theta)*sin(tye);
    ans = sin(phi)*sin(theta)*cos(tye)-cos(theta)*sin(tye)*cos(txe)-cos(phi)*sin(theta)*sin(tye)*cos(txe);

    return(ans);

}//SINphiSINtheta


/***********************************************************************************/

prec COStheta (prec Xe[], prec phi, prec theta)
/* Converts cos(theta), where theta and phi are defined in the coordinate system
where +z lies along the electron beam direction, to an expression in terms of the lab
frame phi and theta. Basically, alters the value of cos(theta) in the lab frame
cross section to account for the electron beam direction.
*/
{
    prec txe;
    prec tye;
    prec ans;

    prec sin_theta, cos_theta, sin_xe, sin_ye, cos_xe, cos_ye;

    txe = XP(Xe);
    tye = YP(Xe);

    sin_theta = sin(theta);
    cos_theta = cos(theta);
    sin_xe = sin(txe);
    sin_ye = sin(tye);
    cos_xe = cos(txe);
    cos_ye = cos(tye);

    //ans = cos_theta*cos_xe*cos_ye + sin_theta*sin(phi)*cos_xe*sin_ye + cos(phi)*sin_theta*sin_xe;
    ans = cos_theta*cos_xe*cos_ye + sin_theta*sin(phi)*sin_ye + cos(phi)*sin_theta*sin_xe*cos_ye;

    return(ans);

}//COStheta

/***********************************************************************************/

prec Photon_Density(prec Xe[], prec t)
/* Calculates the photon density * number of electrons per macro particle
  at the lab frame electron position specified in Xe[] at time t */
{

    prec z2x,z2y,x2,y2,ZR2x,ZR2y,w02x,w02y,z,dz,dose;
    prec Xl[3];
    prec norm;

    LaserFrame(Xe,Xl,t);

    z = Xl[2];
    x2 = Xl[0]*Xl[0];
    y2 = Xl[1]*Xl[1];
    ZR2x = ZRx*ZRx;
    ZR2y = ZRy*ZRy;
    z2x=z*z/ZR2x;
    z2y=z*z/ZR2y;
    w02x = w0x*w0x;
    w02y = w0y*w0y;

    //norm = N_electron*N_photon/( (sqrt(PI*0.5)*PI*0.5)*w0x*w0y*(PULSE_WIDTH/cc)*sqrt(1.0 + z2x)*sqrt(1.0 +z2y) );
    norm = N_electron*N_photon/( (sqrt(PI*0.5)*PI*0.5)*w0x*w0y*(PULSE_WIDTH)*sqrt(1.0 + z2x)*sqrt(1.0 +z2y) );
    dz = (z+cc*t)*(z+cc*t)/(PULSE_WIDTH*PULSE_WIDTH);
    dose = norm*EXP( -2.0*x2/(w02x*(1+z2x)) )*EXP( -2.0*y2/(w02y*(1+z2y)) )*EXP(-2.0*dz);
    return(dose);

}//Photon_Density


/***********************************************************************************/


int LaserFrame(prec Xe[],prec Xl[],prec t)
/* Converts the current lab frame position of the electron (at time t) to a frame whose z axis corresponds
to the negative direction of the laser pulse.
  */
{

    prec x[3];
    prec gamma;
    int i;

    gamma = GAMMA(Xe);
    gamma = 1.0/gamma;

    for (i=0;i<3;i++)
        x[i] = Current_Position(Xe, i, t);

    Xl[0]= x[0]*cos(THETA_0)-x[2]*sin(THETA_0);
    Xl[1]= x[1];
    Xl[2]= x[2]*cos(THETA_0)+x[0]*sin(THETA_0);

    return(0);


}
/***********************************************************************************/

prec Current_Position(prec Xe[], int index, prec t)
/* Given the initial electron position stored X[x,y,z,px,py,pz], this procedure
returns the postion at time t. index: 0->x, 1->y, 2->z*/
{
    prec v;
    prec x;
    prec gamma;

    gamma = GAMMA(Xe);

    v= cc*Xe[index+3]/gamma;

    x=Xe[index]+v*t;

    return(x);

}

/******************************************************************************************/

prec FIND_TIME_FROM_BIN(int bin)
{
    prec dt;
    prec ans;

    dt = FULL_PULSE_LENGTH/prec(TIME_BIN_NUM);

    ans = MINIMUM_ELECTRON + dt*0.5 +prec(bin)*dt - OBS_DISTANCE/cc;

    return (ans);

}//FIND_TIME_FROM_BIN

/******************************************************************************************/

int FIND_TIME_BIN(prec time, prec Xe[],prec phi, prec theta)
/*  Calculates the length from the electron position at time (time) for the observation
point specified by the phi, theta (in rad) and the global variable OBS_DISTANCE */
{
    prec dt,length;
    int ans,i;
    prec x[3];
    prec d[3];

    for(i=0;i<3;i++)
        x[i] = Current_Position(Xe, i, time);

    d[0] = OBS_DISTANCE*tan(theta)*cos(phi);
    d[1] = OBS_DISTANCE*tan(theta)*sin(phi);
    d[2] = OBS_DISTANCE;

    length=0.;
    for(i=0;i<3;i++)
    {
        length = length + (x[i]-d[i])*(x[i]-d[i]);
    }//for
    length = sqrt(length)/cc + time;

    dt = FULL_PULSE_LENGTH/prec(TIME_BIN_NUM);

    ans = int((length- MINIMUM_ELECTRON)/dt);

    return (ans);

}//FIND_TIME_BIN

/***********************************************************************************/

void INITCOR ()
 /* This functions sets the initial particle coordinates and injection times */
{
    ifstream infile;

    char ch;
    char clear[100];
    int nbuf,i;
    prec avg[3];
    avg[2]=0;
    avg[1]=0;
    avg[0]=0;
    prec gamma;
    prec temp;

    AVG_E =0.;
    AVG_GAMMA = 0.;


        infile.open(CORFILE);

        infile>>ch;
        infile>>ch;
        infile>>ch;
        infile>>ch;
        infile>>ch;
        infile>>ch;
        infile>>ch;

        infile>>nbuf;
        //cout<<"Number of macro-particles = "<<nbuf<<endl;

        NUMBER = nbuf + 1;
        infile.getline(clear,99);
        infile.getline(clear,99);
        infile.getline(clear,99);

        if (nbuf>MAX_PARTICLES)
        {
            nbuf = MAX_PARTICLES;
            cout<<"WARNING: Total number of macro-particles exceeds the allocated maximum."<<endl;
            cout<<"Only the first "<<MAX_PARTICLES<<" will be used."<<endl;
        }

        for(i=0;i<nbuf;i++)
        {
            infile>>X[i][0]>>X[i][3]>>X[i][1]>>X[i][4]>>X[i][2]>>X[i][5];
            X[i][2] = DURATION*X[i][2]*PI/180./(2.0*PI*FREQ*TBASE);
            X[i][5] = X[i][5]*ENERGY_SCALE;
            AVG_E = AVG_E + X[i][5]/prec(nbuf);
            infile.getline(clear,99);
            avg[0]=0.;
            avg[1]=0.;
            avg[2]=avg[2]+X[i][2]/nbuf;
//          if(X[i][2]>MAXIMUM_ELECTRON) MAXIMUM_ELECTRON = X[i][2];
//          if(X[i][2]<MINIMUM_ELECTRON) MINIMUM_ELECTRON = X[i][2];
        }//for
//      FULL_PULSE_LENGTH=(MAXIMUM_ELECTRON-MINIMUM_ELECTRON)*cc;
        for (i=0;i<nbuf;i++)
        {
            X[i][5] = (X[i][5]-AVG_E)*DE_MULT + AVG_E;
        }

    infile.close();

    for (i=0;i<nbuf;i++)
    {
        X[i][0]=CONV*X[i][0]*SPOT_SIZE_SCALE/LBASE;
        X[i][1]=CONV*X[i][1]*SPOT_SIZE_SCALE/LBASE;
        X[i][3]=X[i][3]*1.e-3*DIVERGENCE_SCALE + off2[0];
        X[i][4]=X[i][4]*1.e-3*DIVERGENCE_SCALE + off2[1];
        X[i][5]= sqrt(((X[i][5]/REST + 1)*(X[i][5]/REST + 1) -1 )/(X[i][4]*X[i][4] + X[i][3]*X[i][3] + 1));
        X[i][3]=X[i][5]*X[i][3];
        X[i][4]=X[i][5]*X[i][4];
        X[i][2]= X[i][2]  - avg[2];
        gamma=GAMMA(X[i]);
        AVG_GAMMA=AVG_GAMMA+gamma/prec(nbuf);
        X[i][0]= X[i][0] + off[0];
        X[i][1]= X[i][1] + off[1]; //- (-T_START+X[i][2])*(cc*X[i][4]/gamma)
        X[i][2]= -1.0*(X[i][2])*(cc*X[i][5]/gamma) + off[2];
        /* Calculate the Minimum and maximum position of the electrons in normailzed z units at t=0*/
        if(X[i][2]>MAXIMUM_ELECTRON) MAXIMUM_ELECTRON = X[i][2];
        if(X[i][2]<MINIMUM_ELECTRON) MINIMUM_ELECTRON = X[i][2];
        /*Propagate electron transvserve position (takes into account small beta functions)*/
        X[i][0] = X[i][0] + (X[i][3]/X[i][5])*(X[i][2]-off[2]);
        X[i][1] = X[i][1] + (X[i][4]/X[i][5])*(X[i][2]-off[2]);

    }
    FULL_PULSE_LENGTH=(MAXIMUM_ELECTRON-MINIMUM_ELECTRON);
    /* Find time bin range for time resolved intensity and x-ray spectra calculations */
    /* Assume speed of light to detector position on axis, and broaden by 20% to account for any
       broadening due to offaxis x-rays and slow electrons.*/
    temp = MAXIMUM_ELECTRON;
    MAXIMUM_ELECTRON =  (OBS_DISTANCE - MINIMUM_ELECTRON)/cc;
    MINIMUM_ELECTRON = (OBS_DISTANCE - temp)/cc;
    FULL_PULSE_LENGTH=(MAXIMUM_ELECTRON-MINIMUM_ELECTRON);
    temp = MINIMUM_ELECTRON;
    MAXIMUM_ELECTRON = MAXIMUM_ELECTRON + 0.2*FULL_PULSE_LENGTH;
    MINIMUM_ELECTRON = MINIMUM_ELECTRON - 0.2*FULL_PULSE_LENGTH;
    FULL_PULSE_LENGTH=(MAXIMUM_ELECTRON-MINIMUM_ELECTRON);

}//INITCOR

/******************************************************************************************/
void GENCOR(int num_part,prec twiss[], prec avg[])
/*
    Generates a particle coordinate file of the form meant to be read by the procedure INITCOR. num_part
    specifies the number of particles to be generated. twiss contains the 3 emittances and 6 twiss parameters in x, y, and z such that
    twiss[0] = emit_x = unnormalized rms emittance in x (mm-mrad)
    twiss[1] = alpha_x
    twiss[2] = beta_x (mm/mrad)
    twiss[6] = emit_z = rms emittance in z (deg-MeV)
    twiss[7] = alpha_z
    twiss[8] = beta_z (deg/MeV)
    avg contains the average value in all six dimensions (mm mrad mm mrad deg MeV).
    Each particles is coordinate is written to the output file specified by string Gen_Part_File_Out with each column
    specified with
    x(cm)   x'(mrad)   y(cm)  y'(mrad)  t(deg) E(MeV)
*/

{
    int i,j;
    prec** cov_6D;
    long nrow = 6;
    long init_value;
    prec yy[6];
    cov_6D = dmatrix(nrow, nrow);
    init_matrix(cov_6D,nrow,nrow);

    prec gamma_x, gamma_y, gamma_z;
    gamma_x = (1.0+twiss[1]*twiss[1])/twiss[2];
    gamma_y = (1.0+twiss[4]*twiss[4])/twiss[5];
    gamma_z = (1.0+twiss[7]*twiss[7])/twiss[8];

    cov_6D[0][0] = twiss[0]*twiss[2];
    cov_6D[1][1] = twiss[0]*gamma_x;
    cov_6D[1][0] = -1.0*twiss[1];
    cov_6D[0][1] = -1.0*twiss[1];

    cov_6D[2][2] = twiss[3]*twiss[5];
    cov_6D[3][3] = twiss[3]*gamma_y;
    cov_6D[3][2] = -1.0*twiss[4];
    cov_6D[2][3] = -1.0*twiss[4];

    cov_6D[4][4] = twiss[6]*twiss[8];
    cov_6D[5][5] = twiss[6]*gamma_z;
    cov_6D[5][4] = -1.0*twiss[7];
    cov_6D[4][5] = -1.0*twiss[7];


    ofstream out;
    out.open(Gen_Part_File_Out);
    out.setf(LEFT);
    out<<" numbuf= "<<num_part<<endl;
    out<<"       x          xp          y          yp         phi         w    particle #"<<endl;
    out<<endl;

    init_value = Get_ran2_init();

    for (i=0;i<num_part;i++)
    {
        randnCovA(cov_6D,yy,6,&init_value);
        //Convert from mm to cm
        yy[0] = yy[0]*0.1 + avg[0]*0.1;
        yy[1] = yy[1] + avg[1];
        //Convert from mm to cm
        yy[2] = yy[2]*0.1 + avg[2]*0.1;
        yy[3] = yy[3] + avg[3];
        yy[4] = yy[4] + avg[4];
        yy[5] = yy[5] + avg[5];

        for (j=0;j<6;j++)
        {
            out<<setprecision(8)<<setw(20)<<yy[j];
        }
        out<<i<<endl;
    }
    out.close();
}

/****************************************************************************************************/

void COPYARRAY_MIN(prec X[],prec Y[] , int ord)
/* Copys the array X of order ord to Y for each i if X[i]<Y[i] */
{
    int i;
    for(i=0;i<ord;i++)
    {
        if(X[i]<Y[i])
            Y[i] = X[i];
    }
}//COPYARRAY_MIN

/******************************************************************************************/

void COPYARRAY_MAX(prec X[],prec Y[] , int ord)
/* Copys the array X of order ord to Y for each i if X[i]>Y[i] */
{
    int i;
    for(i=0;i<ord;i++)
    {
        if(X[i]>Y[i])
            Y[i] = X[i];
    }
}//COPYARRAY_MAX


/******************************************************************************************/

void COPYARRAY(prec X[],prec Y[] , int ord)
/* Copys the array X of order ord to Y */
{
    int i;
    for(i=0;i<ord;i++)
        Y[i] = X[i];
}//COPYARRAY

/******************************************************************************************/

prec INTERPOLATE (prec x, prec xa[], prec ya[], int n)
/* Given the arrays xa[0..n-1] and ya[0..n-1], which tabulate a function (with
the xa's in order), and given a value fo x, this routine returns the linearly interpolated
value y(x).  */
{
    int klo,khi,k;
    prec h,y;

    klo = 0;
    khi = n-1;
    while (khi-klo > 1)
    {
        k = (khi+klo) >> 1;
        if (xa[k] > x) khi = k;
        else klo =k;
    }
    h = xa[khi]-xa[klo];

    y = ya[klo] + (ya[khi]-ya[klo])/h*(x-xa[klo]);

    return(y);

}//INTERPOLATE


/******************************************************************************************/

int READ_DATA (char filename[], prec xx[], prec yy[], int &count, int max)
/* Reads in data stored in two columns in the file named "filename", and
stored the first column in X and the second column in Y. Make sure no empty
lines at the end of the file. */
{
    int i=0;
    char ch[300];
    int j;

    ifstream infile;
    infile.open(filename);
    while((!infile.eof())&&(i<max))
    {
        infile >> xx[i] >> yy[i];
        infile.getline(ch,299);
        i++;
    }

    count = i;

    for (j=i; j<max; j++)
    {
        xx[j] = INFINITY;
        yy[j] = 0.;
    }


    return(0);
    infile.close();
}

/******************************************************************************************/


prec EXP(prec input)
{
	if (input < MIN_EXPONENT)
	{
		return(0.0);
	}
	else
	{
		return(exp(input));
	}
	
}


/******************************************************************************************/

int INPUT()
/* This procedure reads in data from the input file "Gun3D.ini" to initialize
   key simulation parameters.
*/
{
    ifstream infile;
    infile.open (INFILE);

    char key[50];
    char clear[200];
    int count=0;
    prec temp_var;
    int i;


while ((strcmp(key,"EX")!=0) && (count < 1000))
   {
       infile>>key;
       if( (strcmp(key,"wavelength") == 0)||(strcmp(key,"WAVELENGTH") == 0))
       {
           infile>>WAVELENGTH;
       }//if

       if( (strcmp(key,"frequency") == 0)||(strcmp(key,"FREQUENCY") == 0))
       {
           infile>>FREQ;
       }//if
       else if( (strcmp(key,"input_file") == 0)||(strcmp(key,"INPUT_FILE") == 0))
       {
           infile>>CORFILE;

       }
       else if( (strcmp(key,"out_file") == 0)||(strcmp(key,"OUT_FILE") == 0))
       {
           infile>>OUTFILE;

       }
       else if( (strcmp(key,"LASER_INCIDENCE") == 0)||(strcmp(key,"laser_incidence") == 0))
       {
           infile>>THETA_0;
           THETA_0=THETA_0/180.*PI;
       }
       else if( (strcmp(key,"maxfield") == 0)||(strcmp(key,"MAXFIELD") == 0))
       {
           infile>>A0;
       }
       else if( (strcmp(key,"laser_energy") == 0)||(strcmp(key,"LASER_ENERGY") == 0))
       {
           infile>>Laser_Energy;
           Laser_Energy = Laser_Energy*1000;
//         cout<<"Laesr energy = "<<Laser_Energy<<endl;
           Laser_Energy_Specified = 1;
       }
       else if( (strcmp(key,"FAST_INTENSITY") == 0)||(strcmp(key,"fast_intensity") == 0))
       {
           if(!CHIRP_LASER)
				FAST_INTENSITY = 1;
		   else
			   cout<<"WARNING: FAST_INTENSITY option is not compatible with other inputs, and is therefore disabled."<<endl; 

       }
       else if( (strcmp(key,"INTENSITY") == 0)||(strcmp(key,"intensity") == 0))
       {
           infile>>PHI_NUMBER>>THETA_NUMBER>>THETA_RANGE;
           INTENSITY =1;
           SPECTRUM = 0;
       }
       else if( (strcmp(key,"TIME_STEPS") == 0)||(strcmp(key,"time_steps") == 0))
       {
            infile>>T_STEPS;
       }
       else if( (strcmp(key,"CHARGE") == 0)||(strcmp(key,"charge") == 0))
       {
            infile>>N_electron;
       }
       else if( (strcmp(key,"TIME_RANGE") == 0)||(strcmp(key,"time_range") == 0))
       {
            infile>>INTEGRAL_RANGE_SIGMA_TIME;
       }
	   else if( (strcmp(key,"RANDOMIZE_START") == 0)||(strcmp(key,"randomize_start") == 0)) 
	   {
		   RANDOMIZE_START = 1;
		   Rand_Start_Init_Value_0 = Get_ran2_init();
           Rand_Start_Init_Value = Rand_Start_Init_Value_0;
	   }
       else if( (strcmp(key,"offset") == 0)||(strcmp(key,"OFFSET") == 0))
       {
            infile>>off[0]>>off[1]>>off[2];
       }
       else if( (strcmp(key,"offset2") == 0)||(strcmp(key,"OFFSET2") == 0))
       {
            infile>>off2[0]>>off2[1];
       }
       else if( (strcmp(key,"beam_waist") == 0)||(strcmp(key,"BEAM_WAIST") == 0))
       {
            infile>>w0;
       }
       else if( (strcmp(key,"OBSERVATION") == 0)||(strcmp(key,"observation") == 0))
       {
            infile>>OBSERVE_PHI>>OBSERVE_THETA;
       }
       else if( (strcmp(key,"POLARIZATION") == 0)||(strcmp(key,"polarization") == 0))
       {
            infile>>PHI_P;
            PHI_P = PHI_P/180.*PI;
       }
       else if( (strcmp(key,"PLANE_WAVE") == 0)||(strcmp(key,"plane_wave") == 0))
       {
            infile>>PLANE_WAVE;
       }
       else if( (strcmp(key,"duration") == 0)||(strcmp(key,"DURATION") == 0))
       {
            infile>>DURATION;
       }
       else if( (strcmp(key,"BANDWIDTH_FACTOR") == 0)||(strcmp(key,"bandwidth_factor") == 0))
       {
            infile>>BANDWIDTH_FACTOR;
       }
       else if( (strcmp(key,"energy_scale") == 0)||(strcmp(key,"ENERGY_SCALE") == 0))
       {
            infile>>ENERGY_SCALE;
       }
       else if( (strcmp(key,"DIVERGENCE_SCALE") == 0)||(strcmp(key,"divergence_scale") == 0))
       {
            infile>>DIVERGENCE_SCALE;
       }
       else if( (strcmp(key,"SPOT_SIZE_SCALE") == 0)||(strcmp(key,"spot_size_scale") == 0))
       {
            infile>>SPOT_SIZE_SCALE;
       }
       else if( (strcmp(key,"de_mult") == 0)||(strcmp(key,"DE_MULT") == 0))
       {
            infile>>DE_MULT;
       }
       else if( (strcmp(key,"SPECTRUM_RANGE") == 0)||(strcmp(key,"spectrum_range") == 0))
       {
           infile>>Ws_Min>>Ws_Max>>N_Ws;
           SPECTRUM = 1;
       }
       else if( (strcmp(key,"tbase") == 0)||(strcmp(key,"TBASE") == 0))
       {
            infile>>TBASE;
            TIME_SCALE_SPECIFIED = 1;
       }
       else if( (strcmp(key,"TIME_BIN_NUMBER") == 0)||(strcmp(key,"time_bin_number") == 0))
       {
            infile>>TIME_BIN_NUM;
            if (TIME_BIN_NUM > MAX_TIME_BIN)
                TIME_BIN_NUM=MAX_TIME_BIN;
       }
       else if( (strcmp(key,"SOLID_ANGLE") == 0)||(strcmp(key,"solid_angle") == 0))
       {
           infile>>SOLID_THETA_MIN>>SOLID_THETA_MAX>>SOLID_NTHETA>>SOLID_NPHI;
           SOLID_ANGLE=1;
       }
       else if( (strcmp(key,"FILTER") == 0)||(strcmp(key,"filter") == 0))
       {
           infile>>FILTER_CENTER>>FILTER_WIDTH;
           FILTER_ON = 1;
       }
       else if( (strcmp(key,"FILTER_AMP") == 0)||(strcmp(key,"filter_amp") == 0))
       {
           infile>>FILTER_AMP;
       }
       else if( (strcmp(key,"WINDOW") == 0)||(strcmp(key,"window") == 0))
       {
            WINDOW_ON = 1;
       }
       else if( (strcmp(key,"WINDOW_FILE") == 0)||(strcmp(key,"window_file") == 0))
       {
           infile>>WINDOW_FILE;
           WINDOW_FILE_ON = 1;
       }
       else if( (strcmp(key,"DETECTOR_FILE") == 0)||(strcmp(key,"detector_file") == 0))
       {
           infile>>DETECTOR_FILE;
           DETECTOR_FILE_ON = 1;
       }
       else if( (strcmp(key,"lbase") == 0)||(strcmp(key,"LBASE") == 0))
       {
            infile>>LBASE;
            LENGTH_SCALE_SPECIFIED = 1;
       }
       else if( (strcmp(key,"pulse_width") == 0)||(strcmp(key,"PULSE_WIDTH") == 0))
       {
            infile>>PULSE_WIDTH;
       }
       else if ( (strcmp(key,"CART_INTENSITY")==0)|| (strcmp(key,"cart_intensity") == 0))
       {
           CART_INTENSITY=1;
           infile>>THETA_NUMBER_X>>THETA_RANGE_X>>THETA_NUMBER_Y>>THETA_RANGE_Y;
           if (THETA_NUMBER_X > MAX_PHI_NUMBER)
               THETA_NUMBER_X = MAX_PHI_NUMBER;
           if (THETA_NUMBER_Y > MAX_THETA_NUMBER)
               THETA_NUMBER_Y = MAX_THETA_NUMBER;
       }
       else if ( (strcmp(key,"DATA_DUMP")==0)|| (strcmp(key,"data_dump") == 0))
       {
           DATA_DUMP = 1;
//         infile>>DATA_DUMP_FILE;
       }
       else if( (strcmp(key,"TOTAL_PHOTON") == 0)||(strcmp(key,"total_photon") == 0))
       {
            TOTAL_PHOTON = 1;
       }
       else if( (strcmp(key,"ANGLE_SPECTRUM") == 0)||(strcmp(key,"angle_spectrum") == 0))
       {

           infile>>ANGLE_SPECTRUM>>ANGLE_SPECTRUM_FIXED>>ANGLE_SPECTRUM_START>>ANGLE_SPECTRUM_STOP>>ANGLE_SPECTRUM_NUM;
           if ( (ANGLE_SPECTRUM > 2) || (ANGLE_SPECTRUM < 0 ) ) ANGLE_SPECTRUM = 1;
           if (ANGLE_SPECTRUM_NUM > MAX_TIME_BIN) ANGLE_SPECTRUM_NUM = MAX_TIME_BIN;
           if (ANGLE_SPECTRUM_FIXED > PI ) ANGLE_SPECTRUM_FIXED = ANGLE_SPECTRUM_FIXED/180.0*PI;
       }
       else if( (strcmp(key,"DETECT_POSITION") == 0)||(strcmp(key,"detect_position") == 0))
       {
           infile>>OBS_DISTANCE;
       }
       else if ( (strcmp(key,"CHIRP_LASER")==0) || (strcmp(key,"chirp_laser") == 0))
       {
			infile>>CHIRP_WIDTH;
			CHIRP_LASER = 1;
			TOTAL_PHOTON = 1;
			if(FAST_INTENSITY)
			{
				FAST_INTENSITY = 0;		
				cout<<"WARNING: FAST_INTENSITY option is not compatible with other inputs, and is therefore disabled."<<endl; 
			}
       }
	   else if ( (strcmp(key,"NONLINEAR")==0) || (strcmp(key,"nonlinear") == 0))
	   {
		   if(!(infile>>NONLINEAR_THRESHOLD))
		   {
			   cout<<"ERROR: NONLINEAR must be followed by a nonlinear calculation threshold greater than or equal to zero."<<endl;
			   return(1);
		   }
			NONLINEAR_ON = 1;
			CHIRP_LASER = 1;
			TOTAL_PHOTON = 1;
			if(FAST_INTENSITY)
			{
				FAST_INTENSITY = 0;		
				cout<<"WARNING: FAST_INTENSITY option is not compatible with other inputs, and has therefore been disabled."<<endl; 
			}
	   }
       else if ( (strcmp(key,"THETA_STEPS")==0) || (strcmp(key,"theta_steps") == 0))
       {
           infile>>THETA_STEPS;
       }
       else if ( (strcmp(key,"BRAGG")==0) || (strcmp(key,"bragg") == 0))
       {
           infile>>BRAGG>>BRAGG_SPACING>>BRAGG_ANGLE>>BRAGG_WIDTH;
           BRAGG_ANGLE = BRAGG_ANGLE/180.*PI;
       }
       else if ( (strcmp(key,"PHI_RANGE")==0) || (strcmp(key,"phi_range") == 0))
       {

           infile>>PHI_RANGE_MIN>>PHI_RANGE_MAX;
           if (PHI_RANGE_MAX<PHI_RANGE_MIN)
           {
               temp_var = PHI_RANGE_MIN;
               PHI_RANGE_MIN = PHI_RANGE_MAX;
               PHI_RANGE_MAX=temp_var;
           }
           PHI_RANGE_MIN = PHI_RANGE_MIN/180.*PI;
           PHI_RANGE_MAX = PHI_RANGE_MAX/180.*PI;
       }
       else if ( (strcmp(key,"GEN_PART")==0) || (strcmp(key,"gen_part") == 0))
       {
           GENERATE_PARTICLES = 1;
           infile>>NUM_GEN>>TWISS_6D[0]>>TWISS_6D[1]>>TWISS_6D[2]>>TWISS_6D[3]>>TWISS_6D[4]>>TWISS_6D[5]
               >>TWISS_6D[6]>>TWISS_6D[7]>>TWISS_6D[8]>>AVG_6D[5];
           for (i=0;i<5;i++)
               AVG_6D[i]=0.;
           if (NUM_GEN>MAX_PARTICLES)
               NUM_GEN = MAX_PARTICLES;
           //Convert normalized emittance to unnormailzed emittance
           temp_var = AVG_6D[5]/REST + 1;
           temp_var = sqrt(temp_var*temp_var - 1);
           TWISS_6D[0] = TWISS_6D[0]/temp_var;
           TWISS_6D[3] = TWISS_6D[3]/temp_var;

       }
       else if ( (strcmp(key,"GEN_PART_BASIC")==0) || (strcmp(key,"gen_part_basic") == 0) )
       {
           GENERATE_PARTICLES = 1;
           infile>>NUM_GEN>>TWISS_6D[0]>>TWISS_6D[2]>>TWISS_6D[3]>>TWISS_6D[5]
               >>TWISS_6D[7]>>TWISS_6D[8]>>AVG_6D[5];
           for (i=0;i<5;i++)
               AVG_6D[i]=0.;
           if (NUM_GEN>MAX_PARTICLES)
               NUM_GEN = MAX_PARTICLES;
           //Convert normalized emittance to unnormailzed emittance
           temp_var = AVG_6D[5]/REST + 1;
           temp_var = sqrt(temp_var*temp_var - 1);


           TWISS_6D[0] = TWISS_6D[0]/temp_var;  //emit_x
           TWISS_6D[3] = TWISS_6D[3]/temp_var;  //emit_y
           TWISS_6D[1] = 0.;                   //alpha_x
           TWISS_6D[4] = 0.;                   //alpha_y
           TWISS_6D[2] = (TWISS_6D[2]*TWISS_6D[2])/TWISS_6D[0]; //beta_x
           TWISS_6D[5] = (TWISS_6D[5]*TWISS_6D[5])/TWISS_6D[3]; //beta_y
           TWISS_6D[8] = TWISS_6D[8]*AVG_6D[5]; //energy spread (MeV)
           TWISS_6D[6] = TWISS_6D[8]*TWISS_6D[7];  //emit_z(deg-MeV)
           TWISS_6D[8] = (TWISS_6D[7]*TWISS_6D[7])/TWISS_6D[6];
           TWISS_6D[7] = 0;
       }


       infile.getline(clear,199);
       count++;

    }//while

return(0);
}//INPUT
