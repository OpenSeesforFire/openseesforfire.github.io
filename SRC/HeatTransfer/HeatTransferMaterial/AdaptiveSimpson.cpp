#include <math.h>  // include file for fabs and sin
#include <stdio.h> // include file for printf
#include <fstream>
#include <stdlib.h>
#include <iostream>

using namespace std;
 
// this function does numerical quadrature with adaptive Simpson method 
//
// Recursive auxiliary function for adaptiveSimpsons() function below
//                                                                                                 
double 
adaptiveSimpsonsAux(double (*f)(double), double a, double b, double epsilon,                 
                         double S, double fa, double fb, double fc, int bottom) 
	{                 
	double c = (a + b)/2, h = b - a;                                                                  
	double d = (a + c)/2, e = (c + b)/2;                                                              
	double fd = f(d), fe = f(e);                                                                      
	double Sleft = (h/12)*(fa + 4*fd + fc);                                                           
	double Sright = (h/12)*(fc + 4*fe + fb);                                                          
	double S2 = Sleft + Sright;                                                                       
	if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)                                                    
		return S2 + (S2 - S)/15;                                                                        
	return adaptiveSimpsonsAux(f, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +                    
		adaptiveSimpsonsAux(f, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);                     
	}         
 

//
// Adaptive Simpson's Rule
//
double 
adaptiveSimpsons(double (*f)(double),   // ptr to function
                           double a, double b,  // interval [a,b]
                           double epsilon,  // error tolerance
                           int maxRecursionDepth) // recursion cap
	{           
	double c = (a + b)/2, h = b - a;                                                                  
	double fa = f(a), fb = f(b), fc = f(c);                                                           
	double S = (h/6)*(fa + 4*fc + fb);                                                                
	return adaptiveSimpsonsAux(f, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);                   
	}                                                                                                   
 
double testfunc(double x)
	{
	return (100/(x*x)) * sin(10/x);
	}


//generate the specific heat for steel
double 
getSpecificHeatSteel(double trial_temp)
{
    double cp;
    if (trial_temp < 20.0){
			  cp = 439.7841332;
		} else if ((20.0 <= trial_temp) && (trial_temp < 600.0)) {
			cp = 425.0 + 0.773 * trial_temp - 1.69 * 1e-3 * pow(trial_temp,2) 
				 + 2.22 * 1e-6 * pow(trial_temp,3);
		} else if ((600.0 <= trial_temp) && (trial_temp < 735.0)) {
			cp = 666.0 + 13002.0 / (738.0 - trial_temp);
		} else if ((735.0 <= trial_temp) && (trial_temp < 900.0)) {
			cp = 545.0 + 17820.0 / (trial_temp - 731.0);
		} else if ((900.0 <= trial_temp) && (trial_temp <= 1200.0)) {
			cp = 650.0;
		} else {
			std::cout << "CarbonSteelEC3::getSpecificHeat() - trial temperature " 
				   << trial_temp << " out of bounds, cap values used\n";
			cp = 650.0;
			}
    return cp*7850.0;
}

//generate specific heat for concrete
double moist = 0.03;

double 
getSpecificHeatConcrete(double trial_temp)
{
  double cp; 
  if (trial_temp <= 100.0) {
		cp = 900.0;
		} else if ((100.0 < trial_temp) && (trial_temp <= 200.0)) {
			if (moist == 0.0) {
				cp = trial_temp + 800.0;
				} else if (moist == 0.015) {
					if ((100.0 < trial_temp) && (trial_temp <= 115.0)) {
						cp = 1470.0;
						} else if ((115.0 < trial_temp) && (trial_temp <= 200.0)) {
							cp = -5.529411765 * trial_temp + 2105.882353;
						}
				} else if (moist == 0.03) {
					if ((100.0 < trial_temp) && (trial_temp <= 115.0)) {
						cp = 2020.0;
						} else if ((115.0 < trial_temp) && (trial_temp <= 200.0)) {
							cp = -12.0 * trial_temp + 3400.0;
						} 
					} else {
						cout << "ConcreteEC2::getSpecificHeat() - specific heat values "
							<< "are not available for moisture level " << moist << " .\n";
						exit(-1);
					}
		} else if ((200.0 < trial_temp) && (trial_temp <= 400.0)) {
			cp = 0.5 * trial_temp + 900.0;
			} else {
				cp = 1100.0;
			}

			return cp*2300.0;
}

 
int main(){
/* double I = adaptiveSimpsons(sin, 0, 1, 0.000000001, 10);*/ // compute integral of sin(x)
                                                          // from 0 to 1 and store it in 
                                                          // the new variable I
    double T =10.0;
	double I;

    std::ofstream output;
	//output.open("concreteEnthalpy");

	//Note: For concrete, at T = 1140, the integration using adaptive Simpson quadrature
	//gives inaccurate results, as S(a,b) = S(a,(a+b)/2) + S((a+b)/2,b), this should be
	//pure coincidence
	output.open("concreteEnthalpy");
	for(int i=0;i<120;i++){
		I = adaptiveSimpsons(getSpecificHeatConcrete, 20.0, T, 1e-6, 800);
		output << T << "  " << I << "   " << std::endl;
		T = T + 10.0;
		}
	output.close();

	double EHP = adaptiveSimpsons(getSpecificHeatConcrete, 20.0, 1141, 1e-6, 800);
	printf("EHP = %lf\n",EHP); // print the result
 return 0;
}


//int main(){
// double I = adaptiveSimpsons(sin, 0, 1, 0.000000001, 10); // compute integral of sin(x)
//                                                          // from 0 to 1 and store it in 
//                                                          // the new variable I
// printf("I = %lf\n",I); // print the result
// return 0;
//}