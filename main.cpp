#include <iostream>
#include "bisection_method.h"
#define pi 3.1415926535

using namespace std;

int main()
{
    setlocale(0, "Russian");
    cout << "Задание 5.2. КФ Гаусса, ее узлы и коэффициенты. Вычисление интегралов при помощи КФ Гаусса."
         << endl << "КФ Мелера, ее узлы и коэффициенты. Вычисление интегралов при помощи КФ Мелера" << endl << endl;

    double a, b;
    cout << "Введите нижний предел интегрирования: а = ";
    cin >> a;
    cout << "Введите верхний предел интегрирования: b = ";
    cin >> b;

    for(size_t N = 1; N <= 8; N++)
    {
        cout << "--------------------------------------------------" << endl << "N = " << N << endl;
        cout << "Узлы КФ Гаусса: ";
        double h, x1[10000], x2[10000], y1[10000], y2[10000], x[100], coeff[100];
        size_t counter = 0;
        h = (b-a)/10000;
        x1[counter] = a;
        x2[counter] = x1[counter] + h;
        while (x2[counter] <= b)
        {
            y1[counter] = legendre(x1[counter], N);
            y2[counter] = legendre(x2[counter], N);
            if (y1[counter]*y2[counter] <= 0)
            {
                counter++;
                x1[counter] = x2[counter-1];
                x2[counter] = x2[counter-1] + h;
            }
            else
            {
                x1[counter] = x2[counter];
                x2[counter] = x2[counter] + h;
            }
        }
        //cout << "Число отрезков перемены знака: " << counter << endl << endl;
        size_t counter_2 = 0, i = 1;
        do
        {
            x[i] = bisection_method(x1[counter_2], x2[counter_2], 0.000000000001, N);
            cout << "x[" << i << "] = " << x[i] << "; ";
            counter_2 ++;
            i++;
        } while (counter_2 != counter);

        cout << endl << "Проверка (x[k] + x[N+1-k] = 0): ";
        for (size_t k = 1; k <= N; k++)
        {
            cout << x[k] + x[N+1-k] << "; ";
        }
        cout << endl << endl << "Коэффициенты КФ Гаусса: " << endl;

        for(size_t k = 1; k <= N; k++)
        {
            coeff[k] = (2*(1-x[k]*x[k]))/(N*N*legendre(x[k],N-1)*legendre(x[k],N-1));
            if (coeff[k] <= 0) { cout << "Коэффициент отрицателен!"; return 0;}
            cout << "A[" << k << "] = " << coeff[k] << "; ";
        }
        cout << endl;
        double sum = 0;
        for (size_t k = 1; k <= N; k++)
        {
            cout << "Проверка (C[N+1-" << k << "]-C[" << k << "]=0): " << coeff[N+1-k] - coeff[k] << endl;
            sum += coeff[k];
        }
        cout << "Проверка (sum=2): sum = " << sum << endl;

    }

    cout << "------------------------------------------" << endl << "Рассматриваем функцию f(x)=sqrt(1-sin(x)*sin(x)/2)" << endl;
    cout << "Введите нижний предел интегрирования: а = ";
    cin >> a;
    cout << "Введите верхний предел интегрирования: b = ";
    cin >> b;
    for (size_t N = 5; N <= 8; N++)
    {
        cout << "--------------------------------------------------" << endl << "N = " << N << endl;
        cout << "Узлы КФ Гаусса: ";
        double h, x1[10000], x2[10000], y1[10000], y2[10000], x[100], coeff[100];
        size_t counter = 0;
        h = (b-a)/10000;
        x1[counter] = a;
        x2[counter] = x1[counter] + h;
        while (x2[counter] <= b)
        {
            y1[counter] = legendre(x1[counter], N);
            y2[counter] = legendre(x2[counter], N);
            if (y1[counter]*y2[counter] <= 0)
            {
                counter++;
                x1[counter] = x2[counter-1];
                x2[counter] = x2[counter-1] + h;
            }
            else
            {
                x1[counter] = x2[counter];
                x2[counter] = x2[counter] + h;
            }
        }
        //cout << "Число отрезков перемены знака: " << counter << endl << endl;
        size_t counter_2 = 0, i = 1;
        do
        {
            x[i] = bisection_method(x1[counter_2], x2[counter_2], 0.000000000001, N);
            cout << "x[" << i << "] = " << x[i] << "; ";
            counter_2 ++;
            i++;
        } while (counter_2 != counter);

        /*cout << endl << "Проверка (x[k] + x[N+1-k] = 0): ";
        for (size_t k = 1; k <= N; k++)
        {
            cout << x[k] + x[N+1-k] << "; ";
        }*/
        cout << endl << endl << "Коэффициенты КФ Гаусса: " << endl;

        for(size_t k = 1; k <= counter; k++)
        {
            coeff[k] = (2*(1-x[k]*x[k]))/(N*N*legendre(x[k],N-1)*legendre(x[k],N-1));
            if (coeff[k] <= 0) { cout << "Коэффициент отрицателен!"; return 0;}
            cout << "A[" << k << "] = " << coeff[k] << "; ";
        }
        cout << endl;
        double sum = 0;
        for (size_t k = 1; k <= N; k++)
        {
            sum += coeff[k]*f((b-a)/2*x[k]+(b+a)/2);
        }
        cout << "Значение интеграла: " << (b-a)/2*sum << endl;
    }

    cout << "---------------------------------------" << endl;
    cout << "Рассматриваем КФ Мелера с f(x)=1/(1+x^2)" << endl;
    double N0, x[100]/*, coeff[100]*/;
    for(size_t s = 1; s <= 3; s++)
    {
        cout << "-----------------------------------" << endl << "Введите число узлов: N" << s << " = ";
        cin >> N0;
        cout << "Узлы КФ Мелера: " << endl;
        for (size_t k = 1; k <= N0; k++)
        {
            x[k] = cos((2*k-1)*pi/(2*N0));
            cout << "x[" << k << "] = " << x[k] << "; ";
        }
        //cout << "Коэффициенты КФ Мелера: " << endl;
        /*for(size_t k = 1; k <= N0; k++)
        {
            coeff[k] = pi/N0;
            cout << "A[" << k << "] = " << coeff[k] << "; ";
        }*/
        double sum = 0;
        for(size_t k = 1; k <= N0; k++)
        {
            sum += fMeler(x[k]);
        }
        cout << endl << "Значение интеграла: " << setprecision(12) << (pi*sum)/N0 << endl;
    }



    return 0;
}
