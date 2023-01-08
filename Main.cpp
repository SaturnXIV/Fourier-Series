#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <math.h>

const double x0 = -2.0, xf = 1.0;
const unsigned short c_fourierDegree = 25; //mathematically referred to as n, default degree of 25

const unsigned long int c_integrationIterations = 10000; //default number of subinterval is 10,000
const int c_precision = 5; //default precision of 5 decimal places

const double P = xf - x0;
double a[c_fourierDegree + 1], b[c_fourierDegree + 1];

double linInterpolFofX(unsigned int subIntervals, double* slope, double* knownY, double x)
{
    return knownY[(unsigned int)x] + slope[(unsigned int)x] * (x - (double)(unsigned int)x);
}
double* linInterpolFourierSeries(unsigned int subIntervals, double* slope, double* knownY, unsigned short fourierDeg)
{
    const double P = subIntervals; const unsigned long it = 1000000;
    double* fourierCoefficients = new double[2 * fourierDeg]; //even->a terms & odd-> b terms
    for (unsigned short n = 0; n <= fourierDeg; n++)
    {
        double coef;
        double deltaX = subIntervals / it;
        double sum = 0.0;

        for (unsigned long i = 0; i < it; i++) sum += (linInterpolFofX(subIntervals, slope, knownY, i * deltaX) * cos(2 * M_PI * n * i * deltaX / P) + linInterpolFofX(subIntervals, slope, knownY, i * deltaX + deltaX) * cos(2 * M_PI * n * i * (i * deltaX + deltaX) / P));
        coef = sum * deltaX / 2;
        fourierCoefficients[2*n] = (2.0 / P) * coef;

        sum = 0.0;
        for (unsigned long i = 0; i < it; i++) sum += (linInterpolFofX(subIntervals, slope, knownY, i * deltaX) * sin(2 * M_PI * n * i * deltaX / P) + linInterpolFofX(subIntervals, slope, knownY, i * deltaX + deltaX) * sin(2 * M_PI * n * i * (i * deltaX + deltaX) / P));
        coef = sum * deltaX / 2;
        fourierCoefficients[2*n + 1] = (2.0 / P) * coef;
    }
}

double FofX(double x) //s(x)
{
    //return ((unsigned long)(x) % 2)?(0.0):(1.0); //square wave x0 = 0, xf = 2
    //return x * x * x; //x cubed
    //return x; //sawtooth
    return 1 / (1 + exp(-4 * x)); //sigmoid

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

double integrateFuncOfX(double(*f)(double), double a, double b, unsigned long it = c_integrationIterations) //trapezoid rule implementation
{
    double deltaX = (b - a) / it;
    double sum = 0.0;
    for (unsigned long i = 0; i < it; i++) sum += ((*f)(a + i * deltaX) + (*f)(a + i * deltaX + deltaX));
    return sum * deltaX / 2;
}

double integrateFuncOfXandN(double(*f)(double, unsigned short), double a, double b, unsigned short n, unsigned long it = c_integrationIterations) //trapezoid rule implementation
{
    double coef;
    double deltaX = (b - a) / it;
    double sum = 0.0;
    for (unsigned long i = 0; i < it; i++) sum += ((*f)(a + i * deltaX, n) + (*f)(a + i * deltaX + deltaX, n));
    coef = sum * deltaX / 2;
    return coef;
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
        std::cout << (a[n] > 0.0 ? '+' : ' ') << a[n] << "*cos(" << 2 * M_PI * n / P << "*x)";
        std::cout << (b[n] > 0.0 ? '+' : ' ') << b[n] << "*sin(" << 2 * M_PI * n / P << "*x)";
    }
    std::cout << "\n";
}