#ifndef BISECTION_METHOD_H_INCLUDED
#define BISECTION_METHOD_H_INCLUDED
#include <cmath>
#include <iomanip>

double legendre(double x, size_t N)
{
    double P;
    if (N == 0) P = 1;
    else if (N == 1) P = x;
    else if (N > 1)
        P = ((2*N-1)*legendre(x, N-1)*x)/N - ((N-1)*legendre(x, N-2))/N;
    return P;
}

double f(double x)
{
    return sqrt(1-sin(x)*sin(x)/2);
}

double fMeler (double x)
{
    return 1/(1+x*x);
}

/*Метод бисекции (половинного деления)

    Делим отрезок пополам, пока его длина не будет меньше 2*epsilon.
    Как только мы достигли этого условия, берем X, который находится на середине получившегося отрезка.*/
    double bisection_method (double a, double b, double epsilon, size_t N)
    {
        double c, X/*, delta*/;
        size_t counter = 0;

        while ((b - a) > 2*epsilon)
        {
            c = (a+b)/2;
            if (legendre(a, N)*legendre(c, N) <= 0)
                b = c;
            else
                a = c;
            counter++;
        }
        X = (a + b)/2;
        //delta = (b - a)/2;
        return X;
        /*std::cout << "x = " << std::setprecision(15) << X << std::endl;
        std::cout << "Длина последнего отрезка: " << delta << std::endl;
        std::cout << "Абсолютная величина невязки для приближенного решения: " << fabs(legendre(X) - 0) << std::endl;
        std::cout << "Количество шагов для достижения точности эпсилон: " << counter << std::endl << std::endl;*/
    }


#endif // BISECTION_METHOD_H_INCLUDED
