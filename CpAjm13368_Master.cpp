/***************************************************************************/
/*  TWO-DIMENSIONAL STRUCTURED O-MESH GENERATOR AROUND UNIT CHORD ELLIPSE  */
/*  INPUTS NFAR AS A MULTIPLE OF CHORD LENGTH (TYPICALLY <50C)    				 */
/*  CREATED BY J.M.MACLEOD                                        				 */
/***************************************************************************/

#include <iostream>
#include <ios>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <fstream>

using std::cout;
using std::cin;
using std::endl;

/*** Global Variables ***/
/* Scalars */
bool correct;
int nfar, imax, jmax;
double tc, cellheight, cellratio;

/* Solution arrays */
double *x, *z;

/*** Function Prototypes ***/
/* Standard */

void Initial_Data(void);         // READ IN INPUT DATA
void find_Boundary(void);   		 // FIND THE BOUNDARY VALUES
void TFI(void);									 // USE TRANSFINITE INTERPOLATION TO GENERATE INNER POINTS
void calcnormals(void);					 // CALCULATE NORMALS AT INNER AND OUTER BOUNDARY
int ind(int, int);							 // FUNCTION STEP ALONG THE 1D ARRAYS OF X AND Z.
void techplotwrite(void);				 // WRITE TO TXT FILE TO BE DISPLAYED IN TECPLOT

/*********** END OF FUNCTION PROTOTYPES *************/
/********************************/
/*           MAIN CODE          */
/********************************/

int main(int argc, char *argv[])
{
	int narraylength;

	Initial_Data();


	narraylength=(jmax+1)*(imax+1);
	x = (double *) calloc(narraylength, sizeof(double));
	z = (double *) calloc(narraylength, sizeof(double));

	find_Boundary();

	TFI();

	techplotwrite();

	cout << endl << "Mesh complete, and available, in CpAjm13368_output.dat" << endl;
	cout << "  " <<endl;
	cout << "---------------------------------------------------------------------------"<<endl;
	return 0;

}


/**********************************************************/
void Initial_Data(void)
{
int cellspec, cellcheck, cellcont;
cout << "Computational structured o-mesh generator" << endl;
cout << "for elliptical body with a unit chord. Created by J.O.MacLeod " << endl;
cout << "  " <<endl;
cout << "---------------------------------------------------------------------------"<<endl;
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


cellcheck =nfar/(jmax);
do {
	cout << "How do you wish to specify the first cell?" << endl;
	cout << "1. Cell Height         2. Cell Aspect Ratio" << endl;
	cout << ": ";
	cin >> cellspec;
	correct = true;

  if (cellspec ==1)
  {
    cout << "Input starting cell height"<< endl;
      cout << ": " ;
      cin >> cellheight;
      correct = true;
      if (cellheight<=0)
      {
        cout << "Invalid value, please choose again" << endl;
  			correct = false;
      }
			else if (cellheight>cellcheck) {
				cout << "WARNING: Your starting cell height is greater than Nfar/jmax." << endl;
				cout << "This means points will be clustered at the far field boundary. Do you still want to continue?" << endl;
				cout << "1. Yes					2. No" << endl;
				cout << ":	" ;
				cin >> cellcont;
				if (cellcont == 1) {
					cout << "Continuing..." << endl;
				}
				else if (cellcont == 2) {
					cout << "Generation terminated. Select a new value" << endl;
					correct = false;
				}
				else if (cellcont!=1 && cellcont!=2){
					cout << "Invalid value, please choose again" << endl;
	  			correct = false;
				}
			}
  }
  if (cellspec ==2)
  {
    cout << "Input starting cell aspect ratio"<< endl;
      cout << ": " ;
      cin >> cellratio;
      correct = true;
      if (cellratio<=0)
      {
        cout << "Invalid value, please choose again" << endl;
        correct = false;
      }
  }
  if (cellspec !=1 && cellspec !=2)
    	{
			cout << "Invalid value, please choose again" << endl;
			correct = false;
			}
	} while (correct == false);

}

/* Define the structure & format of the functions */
struct Point {
  double x,z;
  Point(double x, double z) : x(x), z(z) {}
};


/* Find the inner points on the geometry defined */
struct Point calcInner(double theta) {

// Step through using polar coordinate conversion.
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
		cout << "Error: Code attempted to step beyond imax limit" << endl;
		cout << theta << endl;
		cout << "_______________________________________________" << endl;
	}
	  return Point(outerx,outerz);
}

struct Point innernormals(int i)
{
	double innernormalx, innernormalz, mag;

	if (i==1) {
		innernormalx=1.0;
		innernormalz=0.0;
	}
	else if (i==imax) {
 	  innernormalx=1.0;
	  innernormalz=0.0;
  }
	else {
		innernormalx= (z[ind(i+1,1)]-z[ind(i-1,1)]);
		innernormalz= -(x[ind(i+1,1)]-x[ind(i-1,1)]);
		mag = sqrt(pow(innernormalz,2)+pow(innernormalx,2));
		innernormalx = innernormalx/mag;
		innernormalz = innernormalz/mag;


	}
	/* DEBUGGING CODE
	cout<< std::fixed;
	cout << std::setprecision(4);
	cout <<endl;
	cout << innernormalx << "		" << innernormalz << "		" << i<< endl; */
	return Point(innernormalx, innernormalz);
}

struct Point outernormals(double i)
	{
	double outernormalx, outernormalz, mag;

	if (i==1) {
		outernormalx=1.0;
		outernormalz=0.0;
	}
	else if (i==imax) {
		outernormalx=1.0;
		outernormalz=0.0;
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

void find_Boundary(void)
		{
		int i=1;
		double theta, dtheta= (2*M_PI)/(imax-1);

		/*  FIND THE INNER AND OUTER POINTS,
			    ALONG THETA (I) DIRECTION */
		for(theta=0; theta <2*M_PI; theta+=dtheta) {
				Point inner = calcInner(theta);
				Point outer = calcOuter(theta);
				x[ind(i,1)]=inner.x;
				z[ind(i,1)]=inner.z;
				x[ind(i,jmax)]=outer.x;
				z[ind(i,jmax)]=outer.z;
				i++;
		}
x[ind(imax,1)]=x[ind(1,1)];
z[ind(imax,1)]=z[ind(1,1)];
x[ind(imax,jmax)]=x[ind(1,jmax)];
z[ind(imax,jmax)]=z[ind(1,jmax)];
}

void TFI(void) {
	int i,j;
	double zeta, psi01, psi11, psi02, alpha=2.5;

		for(i=1; i<imax; i++) {
			Point innernorm = innernormals(i);
			Point outernorm = outernormals(i);

			for (j=2; j<jmax; j++){
				// Find the values relevant for transfinite interpolation
				zeta=(j-1)*pow((jmax-1),-1);
				psi02=(exp(alpha*zeta)-1-alpha*zeta)/(exp(alpha)-1-alpha);
				psi01=1-psi02;
				psi11=zeta -psi02;

				/* DEBUGGING CODE
				cout<< std::fixed;
				cout << std::setprecision(4);
				cout <<endl;
				cout << j<< "		" << zeta << "		" << psi01 << "		" << psi11 << "		" << psi02 << endl; */
				// Update the values of x and z, assuming they have the same form.
				x[ind(i,j)]=(psi01*x[ind(i,1)]+psi11*innernorm.x+psi02*x[ind(i,jmax)]);
				z[ind(i,j)]=(psi01*z[ind(i,1)]+psi11*innernorm.z+psi02*z[ind(i,jmax)]);

			}

		}
	for (j=2; j<jmax; j++) {
		x[ind(imax,j)]=x[ind(1,j)];
		z[ind(imax,j)]=z[ind(1,j)];
	}
}

int ind(int ivalue, int jvalue)
{
   return ivalue + jvalue*(imax);
}

void techplotwrite(void)
		{
		int i,j,kmax=1;
		double y=0;
		std::ofstream myfile;
		myfile.open ("CpAjm13368_output.dat");
		myfile << "VARIABLES = " << "X " << "Y "<< "Z" << endl;
		myfile << "ZONE I = " << imax << ", J = " << jmax << ", K = "<< kmax << ", F = POINT" << endl;
		myfile << std::fixed;
		myfile << std::setprecision(4);


		    for  (j=1; j<=jmax; j++)
		    {
					for  (i=1; i<=imax; i++)
					{
		        myfile << x[ind(i,j)] << " " << y << " " << z[ind(i,j)] << endl;
		      }
		    }
}
