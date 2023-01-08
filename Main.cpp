#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <math.h>

//Algorithm Describtion
//*
//The algorithm is an implementation of the Fourier series. Given 
//a user-defined function s(x), interval x0 & xf, and a degree n, 
//a Fourier series of degree n, which approximates s(x) is generated 
//and printed to the console.
//*
//by Hicham Ben Abdallah
//version1.0 12/31/2022

const double x0 = -2.0, xf = 1.0;
const unsigned short c_fourierDegree = 25; //mathematically referred to as n, default degree of 25

const unsigned long int c_integrationIterations = 10000; //default number of subintervals is 10,000
const int c_precision = 5; //default precision of 5 decimal places

const double P = xf - x0;
double a[c_fourierDegree + 1], b[c_fourierDegree + 1];

double FofX(double x) //s(x)
{
    //return ((unsigned long)(x) % 2)?(0.0):(1.0); //square wave x0 = 0, xf = 2
    //return x * x * x; //x cubed
    //return x; //sawtooth
    return 1 / (1 + exp(-4 * x)); //sigmoid //default function

    //return /*define your own function*/
}

double cosTerm(double x, unsigned short n)
{
    return FofX(x) * cos(2 * M_PI * n * x / P);
}

double sinTerm(double x, unsigned short n)
{
    return FofX(x) * sin(2 * M_PI * n * x / P);
}

double integrateFuncOfXandN(double(*f)(double, unsigned short), double a, double b, unsigned short n, unsigned long it = c_integrationIterations) //trapezoid rule implementation
{
    double deltaX = (b - a) / it;
    double sum = 0.0;
    for (unsigned long i = 0; i < it; i++) sum += ((*f)(a + i * deltaX, n) + (*f)(a + i * deltaX + deltaX, n));
    return sum * deltaX / 2;
}

int main()
{
    for (unsigned short n = 0; n <= c_fourierDegree; n++)
    {
        a[n] = (2.0 / P) * integrateFuncOfXandN(&cosTerm, x0, xf, n);
        b[n] = (2.0 / P) * integrateFuncOfXandN(&sinTerm, x0, xf, n);
    }

    std::cout.precision(c_precision);
    std::cout << "Fourier series for s(x) of degree " << c_fourierDegree << std::fixed << "\n\n";
    std::cout << a[0] / 2;
    for (unsigned short n = 1; n <= c_fourierDegree; n++)
    {
        std::cout << (a[n] > 0.0 ? '+' : ' ') << a[n] << "*cos(" << 2 * M_PI * n / P << "*x)"; //cosine terms
        std::cout << (b[n] > 0.0 ? '+' : ' ') << b[n] << "*sin(" << 2 * M_PI * n / P << "*x)"; //sine terms
    }
    std::cout << "\n";
}