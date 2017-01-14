/***************************************************************************/
/*  TWO-DIMENSIONAL STRUCTURED O-MESH GENERATOR AROUND UNIT CHORD ELLIPSE  */
/*  INPUTS NFAR AS A MULTIPLE OF CHORD LENGTH (TYPICALLY <50C)    				 */
/*  CREATED BY J.O.MACLEOD                                        				 */
/***************************************************************************/

#include <iostream>
#include <ios>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <fstream>
#include <float.h>
#include <sstream>
#include <string>
#include <string.h>
#include <stdio.h>

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::string;
/*** Global Variables ***/
bool correct;
int nfar, imax, jmax, geometrycheck, NACAseries, reflex, it=0;
double tc=1, cellheight=0, cellratio=0, p, m, r, k1;
std::stringstream NACA;
string NACAnumber;
/* Solution arrays */
double *x, *z, *st;

/*** Function Prototypes ***/
/* Standard */

void Initial_Data(void);         // READ IN INPUT DATA
void find_Boundary(void);   		 // FIND THE BOUNDARY VALUES
void calcnormals(void);					 // CALCULATE NORMALS AT INNER AND OUTER BOUNDARY
void find_cellheight(void);			 // FIND ALPHA ARRAY TO HAVE CONSTANT CELL HEIGHT
void find_cellratio(void);			 // FIND ALPHA ARRAY TO HAVE CONSTANT CELL RATIO
int TFI(int);	    						   // USE TRANSFINITE INTERPOLATION TO GENERATE INNER POINTS
void Mesh_smoothing(void);       // LAPLACE SMOOTHING TO REDUCE THE SHARP EDGES OF TFI
int ind(int, int);							 // FUNCTION STEP ALONG THE 1D ARRAYS OF X AND Z.
void tecplotwrite(void);				 // WRITE TO TXT FILE TO BE DISPLAYED IN TECPLOT

/*********** END OF FUNCTION PROTOTYPES *************/

/********************************/
/*           MAIN CODE          */
/********************************/

int main(void)
{
	int narraylength, i;

	Initial_Data();

/* Create solution arrays */
	narraylength=(jmax+1)*(imax+it+1);
	x = (double *) calloc(narraylength, sizeof(double));
	z = (double *) calloc(narraylength, sizeof(double));
	st = (double *) calloc(imax+it, sizeof(double));
	for (i=0; i<=imax+it; i++ ) {
		st[i]=10;
	}

/* Find the boundary values, both inner and outer */
	find_Boundary();

/* Iterate alpha array for cell height definition */
	if (cellheight!=0) {
	find_cellheight();
  cout << endl << "Mesh Converged. Generating Mesh..." << endl;
	//Perform TFI for the whole grid
	TFI(jmax);
	}

/* Iterate alpha array for cell aspect ratio definition */
	if (cellratio!=0) {
	find_cellratio();
	cout << endl << "Mesh Converged. Generating Mesh..." << endl;
	// Generate the whole grid
	TFI(jmax);
	}

  //cout << endl << "Smoothing mesh..." << endl;
	//Mesh_smoothing();

	cout << endl << "Writing to file..." << endl;
	tecplotwrite();

	cout << endl << "Mesh complete, and available, in CpAjm13368_output.dat" << endl;
	cout << "  " <<endl;
	cout << "---------------------------------------------------------------------------"<<endl;
	return 0;

}


/**********************************************************/
void Initial_Data(void)
{
int cellspec=0;
char NACAn[5];
cout << endl << "---------------------------------------------------------------------------";
cout << endl << "                COMPUTATIONAL STRUCTURED O-MESH GENERATOR" << endl;
cout <<         "                       CREATED BY: J.O.MacLeod " << endl;
cout << "  " <<endl;
cout << "---------------------------------------------------------------------------" << endl;
do {
	cout << "Input the geometry you wish to model" << endl;
	cout << "1. Circle       2. Ellipse    3. NACA 4 or 5 series    4. Input file" << endl;
	cout << ": ";
	cin >> geometrycheck;
	correct=true;
	if (geometrycheck > 4 )
		 {
		 cout << "Invalid value, please choose again" << endl;
		 correct = false;
		 }
} while (correct == false);

if(geometrycheck==3) {
	do{
	cout<< "Input Aerofoil Number" << endl;
	cout << ": ";
	cin >> NACAn;
  correct=true;
	NACAseries=(unsigned)strlen(NACAn);
	//cout << NACAseries;
  NACA << NACAn;
	NACA >> NACAnumber;
	if (NACAseries <4 || NACAseries > 5) {
		cout << "Invalid value, please choose again" << endl;
		correct = false;
		}
	//Assign values for the calculation of the aerofoil.
	//Find the 4 digit series parmameters
	if (NACAseries==4){
		//Find the thickness of the aerofoil in % of chord.
		string thickness=NACAnumber.substr(2,2);
		std::stringstream(thickness)>>tc;
		thickness=NACAnumber.substr(0,1);
		//Find the % camber
		std::stringstream(thickness)>>p;
		p=p/100;
		thickness=NACAnumber.substr(1,1);
		//Find the location of max camber (in 1/10ths chord)
		std::stringstream(thickness)>>m;
		m=m/10;
		//cout << p << " " << m << endl;
	}
	//Find the 5 digit series values
	if(NACAseries==5) {
		string thickness=NACAnumber.substr(3,2);
		std::stringstream(thickness)>>tc;

		thickness=NACAnumber.substr(0,1);
		std::stringstream(thickness)>>p;
		p=p*3/20;

		thickness=NACAnumber.substr(1,1);
		std::stringstream(thickness)>>m;
		m=m/20;

		thickness=NACAnumber.substr(2,1);
		std::stringstream(thickness)>>reflex;

		if (reflex==0) {
			if (m==0.05) {
				r=0.0580;
				k1=361.4;
			}
			if (m==0.10) {
				r=0.1260;
				k1=51.640;
			}
			if (m==0.15) {
				r=0.2025;
				k1=15.957;
			}
			if (m==0.20) {
				r=0.2900;
				k1=6.643;
			}
			if (m==0.25) {
				r=0.3910;
				k1=3.230;
			}
		}

		if (reflex==1) {
			if (m==0.10) {
				r=0.1300;
				k1=51.990;
			}
			if (m==0.15) {
				r=0.2170;
				k1=15.793;
			}
			if (m==0.20) {
				r=0.3180;
				k1=6.520;
			}
			if (m==0.25) {
				r=0.4410;
				k1=3.191;
			}
		}
		k1=k1*p*20/3;
		//cout << p << " " << m << " " << k1 << " " << r << endl;
	}

	}while (correct==false);
}

if (geometrycheck==2){
	do {
		cout << "Input the thickness to chord ratio, t/c" << endl;
		cout << ": ";
		cin >> tc;
		correct = true;
		 if (tc >1)
				{
				cout << "Invalid value, please choose again" << endl;
				correct = false;
				}
	  if (tc <=0)
	     {
	     cout << "Invalid value, please choose again" << endl;
	     correct = false;
	     }
	} while (correct == false);
}

if (geometrycheck==4) {
	std::string oss;
	cout << "Input file name" << endl;
	getline(cin, oss);
	std::ifstream infile;
	infile.open(oss.c_str());
	cout << "File opened successfully. Actually using it to be added..." << endl;
}

do {
cout << "Input the distance, Nfar, for the mesh to end (must be an integer). Minimum of 5" << endl;
	cout << ": ";
	cin >> nfar;
	correct = true;
	 if (nfar < 5 )
			{
			cout << "Invalid value, please choose again" << endl;
			correct = false;
			}
} while (correct == false);

do {
	cout << "Input the number of points in the i direction (minimum of 10)" << endl;
	cout << ": ";
	cin >> imax;
	correct = true;
	 if (imax < 10 )
			{
			cout << "Invalid value, please choose again" << endl;
			correct = false;
			}
} while (correct == false);

do {
	cout << "Input the number of points in the trailing edge truncation" << endl;
	cout << ": ";
	cin >> it;
	correct = true;
	 if (it < 1 )
			{
			cout << "Invalid value, please choose again" << endl;
			correct = false;
			}
} while (correct == false);

do {
	cout << "Input the number of points in the j direction (minimum of 10)" << endl;
	cout << ": ";
	cin >> jmax;
	correct = true;
	 if (jmax < 10 )
			{
			cout << "Invalid value, please choose again" << endl;
			correct = false;
			}
} while (correct == false);



do {
	cout << "How do you wish to specify the first cell?" << endl;
	cout << "1. Cell Height         2. Cell Aspect Ratio" << endl;
	cout << ": ";
	cin >> cellspec;
	correct = true;
	if (cellspec !=1 && cellspec !=2)
			{
			cout << "Invalid value, please choose again" << endl;
			correct = false;
			}
  else if (cellspec ==1)
  {
		do {
    	cout << "Input starting cell height"<< endl;
      cout << ": " ;
      cin >> cellheight;

			correct = true;
      if (cellheight<=0)
      {
        cout << "Invalid value, please choose again" << endl;
  			correct = false;
      }
		} while (correct==false);

  }
  else if (cellspec ==2)
  {
		do {
				cout << "Input starting cell aspect ratio"<< endl;
	      cout << ": " ;
	      cin >> cellratio;
	      correct = true;
	      if (cellratio<=0)
	      {
	        cout << "Invalid value, please choose again" << endl;
	        correct = false;
	      }
			} while (correct==false);
		}
// Keep performing the loop until a correct value is given
} while (correct == false);





}

/* Define the structure & format of the functions */
struct Point {
  double x,z;
  Point(double x, double z) : x(x), z(z) {}
};


/* Find the inner points on the geometry defined */
struct Point calcInner(double theta) {

	double innerx=0.5 * cos(theta);
	double innerz= tc/2.0 * sin(theta);

	return Point(innerx,innerz);
}

/* Find the outer points of the mesh, on the boundary */
struct Point calcOuter(double theta) {
  double outerx;
  double outerz;
	// Step through to find the correct face, then calculate the lengths
	if (theta<=M_PI/4) {
			outerx=nfar;
			outerz=nfar*tan(theta);
	}

	else if (theta>M_PI/4 && theta<=M_PI*3/4){
			outerx=-nfar*tan(theta-M_PI/2);
			outerz=nfar;
	}

	else if (theta>M_PI*3/4 && theta<= M_PI*5/4) {
		  outerx=-nfar;
			outerz=-nfar*tan(theta);
	}
	else if (theta>M_PI*5/4 && theta<=M_PI*7/4) {
			outerx=nfar*tan(theta-M_PI/2);
			outerz=-nfar;
	}
	else if (theta>M_PI*7/4 && theta <=2*M_PI) {
			outerx=nfar;
			outerz=nfar*tan(theta);
	}
	else if (theta>2*M_PI) {
		cout << "ERROR: Code attempted to step beyond imax limit" << endl;
		cout << theta << endl;
		cout << "_______________________________________________" << endl;
	}
	  return Point(outerx,outerz);
}

struct Point innernormals(int i)
	{
		double innernormalx, innernormalz, mag;
// first and last cells will always be [1,0]T. No i-1 or i+1 point, so simple define.
		if (i==1) {
			innernormalx=1.0;
			innernormalz=0.0;
		}
		else if (i==imax) {
	 	  innernormalx=1.0;
		  innernormalz=0.0;
	  }
		// find the normals through i, from central difference, and vector transformation.
		else {
			innernormalx= (z[ind(i+1,1)]-z[ind(i-1,1)]);
			innernormalz= -(x[ind(i+1,1)]-x[ind(i-1,1)]);
			mag = sqrt(pow(innernormalz,2)+pow(innernormalx,2));
			innernormalx = innernormalx/mag;
			innernormalz = innernormalz/mag;


		}
	return Point(innernormalx, innernormalz);
	}

struct Point outernormals(double i)
	{
	double outernormalx, outernormalz, mag;
// first and last cells will always be [1,0]T. No i-1 or i+1 point, so simple define.
	if (i==1) {
		outernormalx= (z[ind(i+1,jmax)]-z[ind(imax+it-1,jmax)]);
		outernormalz= -(x[ind(i+1,jmax)]-x[ind(imax+it-1,jmax)]);
		mag = sqrt(pow(outernormalz,2.0)+pow(outernormalx,2.0));
		outernormalx = outernormalx/mag;
		outernormalz = outernormalz/mag;
	}
	else if (i==imax+it) {
		outernormalx= (z[ind(2,jmax)]-z[ind(i-1,jmax)]);
		outernormalz= -(x[ind(2,jmax)]-x[ind(i-1,jmax)]);
		mag = sqrt(pow(outernormalz,2.0)+pow(outernormalx,2.0));
		outernormalx = outernormalx/mag;
		outernormalz = outernormalz/mag;
	}
	else {
		outernormalx= (z[ind(i+1,jmax)]-z[ind(i-1,jmax)]);
		outernormalz= -(x[ind(i+1,jmax)]-x[ind(i-1,jmax)]);
		mag = sqrt(pow(outernormalz,2.0)+pow(outernormalx,2.0));
		outernormalx = outernormalx/mag;
		outernormalz = outernormalz/mag;

	}

	return Point(outernormalx, outernormalz);
}

void find_Boundary(void)	{
		int i=1;
		double theta, dtheta= (2*M_PI)/(imax-1), k2;
		double trailinglengthx, trailinglengthz, outerlengthx, outerlengthz;

		/*  FIND THE INNER AND OUTER POINTS,
			    ALONG THETA (I) DIRECTION */
	  if(geometrycheck==1 || geometrycheck==2){
			for(theta=0; theta <=2*M_PI; theta+=dtheta) {
					Point inner = calcInner(theta);
					Point outer = calcOuter(theta);
					x[ind(i,1)]=inner.x;
					z[ind(i,1)]=inner.z;
					x[ind(i,jmax)]=outer.x;
					z[ind(i,jmax)]=outer.z;
					i++;
			}

		}

		if(geometrycheck==3) {
			i=1;
			double a0=0.2969, a1=-0.1260, a2=-0.3516, a3=0.2843, a4=-0.1015;

			for (theta=0; theta<M_PI; theta+=dtheta) {
			Point inner = calcInner(theta);
			Point outer = calcOuter(theta);
			x[ind(i,1)]=inner.x;
			z[ind(i,1)]=(tc/20)*(a0*pow(inner.x+0.5,0.5)+a1*(inner.x+0.5)
			+a2*pow(inner.x+0.5,2)+a3*pow(inner.x+0.5,3)+a4*pow(inner.x+0.5,4));

			x[ind(i,jmax)]=outer.x;
			z[ind(i,jmax)]=outer.z;
			i++;
			}

			for (theta=theta; theta<2*M_PI ; theta+=dtheta ) {
				Point inner = calcInner(theta);
				Point outer = calcOuter(theta);
				x[ind(i,1)]=inner.x;
				z[ind(i,1)]=(-tc/20)*(a0*pow(inner.x+0.5,0.5)+a1*(inner.x+0.5)
				+a2*pow(inner.x+0.5,2)+a3*pow(inner.x+0.5,3)+a4*pow(inner.x+0.5,4));

				x[ind(i,jmax)]=outer.x;
				z[ind(i,jmax)]=outer.z;
				i++;
			}

			Point inner = calcInner(2*M_PI);
			Point outer = calcOuter(2*M_PI);
			x[ind(imax,1)]=inner.x;
			z[ind(imax,1)]=(-tc/20)*(a0*pow(inner.x+0.5,0.5)+a1*(inner.x+0.5)
			+a2*pow(inner.x+0.5,2)+a3*pow(inner.x+0.5,3)+a4*pow(inner.x+0.5,4));

			x[ind(imax,jmax)]=outer.x;
			z[ind(imax,jmax)]=outer.z;

			if (NACAseries==4) {
				//Camber line addition for the the rear.
				 i=1;
				do {
					z[ind(i,1)]+=p/pow(1-m,2)*(1-2*m+2*m*(x[ind(i,1)]+0.5)-pow(x[ind(i,1)]+0.5,2));
					z[ind(imax-i,1)]+=p/pow(1-m,2)*(1-2*m+2*m*(0.5+x[ind(imax-i,1)])-pow(0.5+x[ind(imax-i,1)],2));
					//cout << x[ind(i,1)] << " " << z[ind(i,1)] << endl;
					i++;
				} while (x[ind(i,1)]+0.5 >= m);
				//Camber line addition for the front.
				do {
					z[ind(i,1)]+=p/pow(m,2)*(2*m*(x[ind(i,1)]+0.5)-pow(x[ind(i,1)]+0.5,2));
					z[ind(imax-i,1)]+=p/pow(m,2)*(2*m*(x[ind(imax-i,1)]+0.5)-pow(x[ind(imax-i,1)]+0.5,2));
					i++;
				} while (i < imax-i);
			}


			if (NACAseries==5) {
				if (reflex==0) {
					i=1;
					do {
						z[ind(i,1)]+=k1*pow(r,3)/6*(1-(x[ind(i,1)]+0.5));
						z[ind(imax-i,1)]+=k1*pow(r,3)/6*(1-(x[ind(imax-i,1)]+0.5));
						//cout << x[ind(i,1)] << " " << z[ind(i,1)] << endl;
						i++;
					} while (x[ind(i,1)]+0.5 >= r);
					//Camber line addition for the front.
					do {
						z[ind(i,1)]+=k1/6*(pow(x[ind(i,1)]+0.5,3)-3*r*pow(x[ind(i,1)]+0.5,2)
						+pow(r,2)*(3-r)*(x[ind(i,1)]+0.5));

						z[ind(imax-i,1)]+=k1/6*(pow(x[ind(imax-i,1)]+0.5,3)-3*r*pow(x[ind(imax-i,1)]+0.5,2)
						+pow(r,2)*(3-r)*(x[ind(imax-i,1)]+0.5));

						i++;
					} while (i < imax-i);
				}
				if (reflex==1) {
					k2= 3/(1-r)*pow(r-m,2)-pow(r,3);
					i=1;
					do {
						z[ind(i,1)]+=k1/6*(k2*pow(x[ind(i,1)]+0.5-r,3)-k2*pow(1-r,3)*
						(x[ind(i,1)]+0.5)-pow(r,3)*(x[ind(i,1)]+0.5-1));

						z[ind(imax-i,1)]+=k1/6*(k2*pow(x[ind(i,1)]+0.5-r,3)-k2*pow(1-r,3)*
						(x[ind(i,1)]+0.5)-pow(r,3)*(x[ind(i,1)]+0.5-1));
						//cout << x[ind(i,1)] << " " << z[ind(i,1)] << endl;
						i++;
					} while (x[ind(i,1)]+0.5 >= r);
					//Camber line addition for the front.
					do {
						z[ind(i,1)]+=k1/6*(pow(x[ind(i,1)]+0.5-r,3)-k2*pow(1-r,3)*
						(x[ind(i,1)]+0.5)-pow(r,3)*(x[ind(i,1)]+0.5-1));
						z[ind(imax-i,1)]+=k1/6*(pow(x[ind(imax-i,1)]+0.5-r,3)-k2*pow(1-r,3)*
						(x[ind(imax-i,1)]+0.5)-pow(r,3)*(x[ind(imax-i,1)]+0.5-1));
						i++;
					} while (i < imax-i);
				}
			}
		}

		trailinglengthx= (x[ind(1,1)]-x[ind(imax,1)])/it;
		trailinglengthz = (z[ind(1,1)]-z[ind(imax,1)])/it;
		outerlengthx=(x[ind(1,jmax)]-x[ind(imax,jmax)])/it;
		outerlengthz=(z[ind(1,jmax)]-z[ind(imax,jmax)])/it;
		for (i=1; i<=it; i++) {
		// for imax, also on the trailing edge, so superimpose the first cell
			x[ind(imax+i,1)]=x[ind(imax,1)]+i*trailinglengthx;
			z[ind(imax+i,1)]=z[ind(imax,1)]+i*trailinglengthz;
			x[ind(imax+i,jmax)]=x[ind(imax,jmax)]+i*outerlengthx;
			z[ind(imax+i,jmax)]=z[ind(1,jmax)]+i*outerlengthz;
		}
}

void find_cellheight(void) {
	int i, test=0;
	double rms=2, *error, errorxj, errorzj, lengthj;
	error = (double *) calloc(imax+it, sizeof(double));

	while (fabs(rms-1)>0.0001) {
		rms=0;
		if (test>500) {
			break;
		}
 	 /* Peform the TFI code calculating just the 2nd ring to find alpha array */
 	 TFI(3);

 	 for (i=1; i<=imax+it; i++) {
 		 //cell components in the j direction
 		 errorxj=x[ind(i,2)]-x[ind(i,1)];
 		 errorzj=z[ind(i,2)]-z[ind(i,1)];
 		 lengthj=sqrt(pow(errorxj,2)+pow(errorzj,2));
 		 //Compute the ratio of cell height compared to the desired.
 		 error[i]=lengthj/cellheight;
 		 //adjust each st accordingly: Increase in st -> decrease in cell height.
 		 st[i]=st[i] + (error[i]-1);
 		 rms+= pow(error[i],2);
 	 }

 	 rms=sqrt(rms/imax);
 	 cout << "Error:   " << fabs(rms-1) <<  endl;
	 test++;
  }

}

void find_cellratio(void) {
	int i, test=0;
	double rms=2, *error, errorxj, errorzj, errorxi, errorzi, lengthj, lengthi;
	error = (double *) calloc(imax+it, sizeof(double));

	while (fabs(rms-1)>0.0001) {
		rms=0;
		if (test>500) {
			break;
		}
		/* Peform the TFI code calculating just the 2nd ring to find alpha array */
		TFI(3);
		for (i=1; i<imax+it; i++) {
		//Cell components in the j direction
		errorxj=x[ind(i,2)]-x[ind(i,1)];
		errorzj=z[ind(i,2)]-z[ind(i,1)];
		lengthj=sqrt(pow(errorxj,2)+pow(errorzj,2));
		// Cell components in the i direction (body face)
		errorxi=x[ind(i+1,1)]-x[ind(i,1)];
		errorzi=z[ind(i+1,1)]-z[ind(i,1)];
		lengthi=sqrt(pow(errorxi,2)+pow(errorzi,2));
		//Compute the ratio of current aspect ratio compared to desired.
		error[i]=(lengthj/lengthi)/cellratio;
		//adjust each st accordingly: Increase in st -> decrease in cell ratio.
		st[i]=st[i] + (error[i]-1);
		rms+= pow(error[i],2);
		}
		rms=sqrt(rms/(imax+it-1));
		cout << "Error:   " << fabs(rms-1) << endl;
	}

}

// Function to perform improved transfinite interpolation [C.B.Allen, 2008].
int TFI(int jvalue) {
	int i,j;
	double zeta, psi1, psi2, psi3, psi4, st2=1/2;
	//double psi5, psi6;
	double fnOuterx, fnOuterz;
	//double fnInnerx, fnInnerz;

	for (j=2; j<jvalue; j++){
		for(i=1; i<imax+it; i++) {
				Point innernorm = innernormals(i);
				fnOuterx=x[ind(i,1)]+fabs(x[ind(i,jmax)]-x[ind(i,1)])*innernorm.x;
				fnOuterz=z[ind(i,1)]+fabs(z[ind(i,jmax)]-z[ind(i,1)])*innernorm.z;
				/* OPTIONS FOR INCLUDING ORTHOGONALITY AT THE BOUNDARY */
				//Point outernorm = outernormals(i);
				//fnInnerx=x[ind(i,jmax)]-fabs(x[ind(i,jmax)]-x[ind(i,1)])*outernorm.x;
				//fnInnerz=z[ind(i,jmax)]-fabs(z[ind(i,jmax)]-z[ind(i,1)])*outernorm.z;

				// Find the values relevant for transfinite interpolation
				zeta=(j-1)*pow((jmax-1),-1);
				psi1=-(tanh(st[i]*(zeta-1)))/tanh(st[i]);
				psi2=1-psi1;
				psi3=pow(psi2,st2);
				psi4=1-psi3;
				//psi5=1-pow(psi2,2);
				//psi6=1-psi5;


				// Update the values of x and z, assuming they have the same form.
				x[ind(i,j)]=(psi1*x[ind(i,1)]+psi2*(psi4*fnOuterx+psi3*x[ind(i,jmax)]));
				z[ind(i,j)]=(psi1*z[ind(i,1)]+psi2*(psi4*fnOuterz+psi3*z[ind(i,jmax)]));

			}

		}
	//map i=1 to i=imax points. j=1 and j=jmax have already been assigned.
	for (j=2; j<jmax; j++) {
		x[ind(imax+it,j)]=x[ind(1,j)];
		z[ind(imax+it,j)]=z[ind(1,j)];
	}
	return 0;
}

/* Function to smooth out the TFI edges. Laplace smoothing risks smoothing into the surface. */
void Mesh_smoothing(void) {
	int i=1, j=1, narraylength;
	double riplusx, riplusz, test=0;

	double *tempx, *tempz;

	narraylength=(jmax)*(imax+it);
	tempx = (double *) calloc(narraylength, sizeof(double));
	tempz = (double *) calloc(narraylength, sizeof(double));



		//create a loop of smoothing stepping proportional to alpha (bigger alpha=less smooth mesh)
	while(test<0) {

		for (i=2; i<imax+it; i++ ) {
			for (j=3; j<jmax; j++) {
				//FIND THE VECTOR ALONG I OF EACH CELL
				riplusx=x[ind(i+1,j)]-x[ind(i,j)];
				riplusz=z[ind(i+1,j)]-z[ind(i,j)];

				//CHECK FOR NEGATIVE VOLUME CELLS. CONTINUE SMOOTHING IF SO.
				if (i<1.0*imax/4) {
					if (riplusx >=  0) {
						test = 0;
					}
					else if (riplusz <= 0) {
						test = 0;
					}
					else {
						test=test;
					}
				}
				//PERFORM LAPLACIAN SMOOTHING
		  	tempx[ind(i,j)]=(x[ind(i+1,j)]+x[ind(i-1,j)]+x[ind(i,j+1)]+x[ind(i,j-1)])/4;
				tempz[ind(i,j)]=(z[ind(i+1,j)]+z[ind(i-1,j)]+z[ind(i,j+1)]+z[ind(i,j-1)])/4;
			}
		}
		for (j=2; j<jmax; j++) {
			//PEFRORM LAPLACE SMOOTHING AT THE TRAILING EDGE, WHERE THERE IS NO I-1 POINT.
			x[ind(1,j)]=(x[ind(2,j)]+x[ind(imax-1,j)]+x[ind(1,j+1)]+x[ind(1,j-1)])/4;
			z[ind(1,j)]=(z[ind(2,j)]+z[ind(imax-1,j)]+z[ind(1,j+1)]+z[ind(1,j-1)])/4;
			//REPLACE IMAX POINT AS THE SAME AS I=1 POINT.
			x[ind(imax,j)]=x[ind(1,j)];
			z[ind(imax,j)]=z[ind(1,j)];
		}
		for (i=2; i<imax+it; i++) {
			for (j=3; j<jmax; j++) {
				x[ind(i,j)]=tempx[ind(i,j)];
				z[ind(i,j)]=tempz[ind(i,j)];
			}
		}

	  test++;
	}

}

//Index function to find the appropriate x and z value along the 1 arrays.
int ind(int ivalue, int jvalue)
{
   return ivalue + jvalue*(imax+it);
}
//WRITE TO FILE
void tecplotwrite(void)
		{
		int i,j,kmax=1;
		double y=0;
		std::ofstream myfile;
		myfile.open ("CpAjm13368_output.dat");
		myfile << "VARIABLES = " << "X " << "Y "<< "Z" << endl;
		myfile << "ZONE I = " << imax+it << ", J = " << jmax << ", K = "<< kmax << ", F = POINT" << endl;
		// Set the number of decimal places to 4
		myfile << std::fixed;
		myfile << std::setprecision(8);

    // Write values to text file.
		    for  (j=1; j<=jmax; j++)
		    {
					for  (i=1; i<=imax+it; i++)
					{
		        myfile << x[ind(i,j)] << " " << y << " " << z[ind(i,j)] << endl;
		      }
		    }
}
