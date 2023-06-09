#include <iostream>
#include <boost/math/quadrature/gauss.hpp>
#include "bisection_method.h"

using namespace std;
double pi = acos(-1.0);

int main()
{
    setlocale(0, "Russian");
    cout << "������� 5.2. �� ������, �� ���� � ������������. ���������� ���������� ��� ������ �� ������."
         << endl << "�� ������, �� ���� � ������������. ���������� ���������� ��� ������ �� ������" << endl << endl;

    double a, b, a1, b1;
    /*cout << "������� ������ ������ ��������������: � = ";
    cin >> a;
    cout << "������� ������� ������ ��������������: b = ";
    cin >> b;*/
    cout << "������ ������ ��������������: � = -1" << endl << "������� ������ ��������������: b = 1" << endl;
    a1 = -1;
    b1 = 1;

    for(size_t N = 1; N <= 8; N++)
    {
        cout << "--------------------------------------------------" << endl << "N = " << N << endl;
        cout << "���� �� ������: ";
        double h, x1[10000], x2[10000], y1[10000], y2[10000], x[100], coeff[100];
        size_t counter = 0;
        h = (b1-a1)/10000;
        x1[counter] = a1;
        x2[counter] = x1[counter] + h;
        while (x2[counter] <= b1)
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
        //cout << "����� �������� �������� �����: " << counter << endl << endl;
        size_t counter_2 = 0, i = 1;
        do
        {
            x[i] = bisection_method(x1[counter_2], x2[counter_2], 0.000000000001, N);
            cout << "x[" << i << "] = " << x[i] << "; ";
            counter_2 ++;
            i++;
        } while (counter_2 != counter);

        cout << endl << "�������� (x[k] + x[N+1-k] = 0): ";
        for (size_t k = 1; k <= N; k++)
        {
            cout << x[k] + x[N+1-k] << "; ";
        }
        cout << endl << endl << "������������ �� ������: " << endl;

        for(size_t k = 1; k <= N; k++)
        {
            coeff[k] = (2*(1-x[k]*x[k]))/(N*N*legendre(x[k],N-1)*legendre(x[k],N-1));
            if (coeff[k] <= 0) { cout << "����������� �����������!"; return 0;}
            cout << "A[" << k << "] = " << coeff[k] << "; ";
        }
        cout << endl;
        double sum = 0;
        for (size_t k = 1; k <= N; k++)
        {
            cout << "�������� (C[N+1-" << k << "]-C[" << k << "]=0): " << coeff[N+1-k] - coeff[k] << endl;
            sum += coeff[k];
        }
        cout << "�������� (sum=2): sum = " << sum << endl;

    }

    cout << "------------------------------------------" << endl << "������������� ������� f(x)=sqrt(1-sin(x)*sin(x)/2)" << endl;
    cout << "������� ������ ������ �������������� (��:0): � = ";
    cin >> a;
    cout << "������� ������� ������ �������������� (��:1.5708): b = ";
    cin >> b;
    using boost::math::quadrature::gauss;
    auto exact = [](double x) { return sqrt(1-sin(x)*sin(x)/2); };
    cout << "������ �������� ���������: " << gauss<double, 7>::integrate(exact, a, b) << endl;

    for (size_t N = 5; N <= 8; N++)
    {
        cout << "--------------------------------------------------" << endl << "N = " << N << endl;
        cout << "���� ��: ";
        double h, x1[10000], x2[10000], y1[10000], y2[10000], x[100], coeff[100], coeffA[100];
        size_t counter = 0;
        h = (b1-a1)/10000;
        x1[counter] = a1;
        x2[counter] = x1[counter] + h;
        while (x2[counter] <= b1)
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
        //cout << "����� �������� �������� �����: " << counter << endl << endl;
        size_t counter_2 = 0, i = 1;
        do
        {
            x[i] = bisection_method(x1[counter_2], x2[counter_2], 0.000000000001, N);
            cout << "x[" << i << "] = " << x[i] << "; ";
            counter_2 ++;
            i++;
        } while (counter_2 != counter);

        //������������ �� ������
        for(size_t k = 1; k <= N; k++)
        {
            coeff[k] = (2*(1-x[k]*x[k]))/(N*N*legendre(x[k],N-1)*legendre(x[k],N-1));
            if (coeff[k] <= 0) { cout << "����������� �����������!"; return 0;}
        }

        cout << endl << endl << "������������ �� �������� �� ������: " << endl;

        for(size_t k = 1; k <= N; k++)
        {
            coeffA[k] = (b - a) * coeff[k] / 2;
            if (coeffA[k] <= 0) { cout << "����������� �����������!"; return 0;}
            cout << "A[" << k << "] = " << coeffA[k] << "; ";
        }
        cout << endl;
        double sum = 0;
        for (size_t k = 1; k <= N; k++)
        {
            sum += coeffA[k]*f((b-a)*x[k]/2+(b+a)/2);
        }
        cout << "�������� ���������: " << sum << endl;
    }

    cout << endl << endl << "---------------------------------------" << endl;
    cout << "������������� �� ������ � f(x)=1/(1+x^2)" << endl;

    double N0, x[100];
    for(size_t s = 1; s <= 3; s++)
    {
        cout << "-----------------------------------" << endl << "������� ����� �����: N" << s << " = ";
        cin >> N0;
        cout << "���� �� ������: " << endl;
        for (size_t k = 1; k <= N0; k++)
        {
            x[k] = cos((2*k-1)*pi/(2*N0));
            cout << "x[" << k << "] = " << x[k] << "; ";
        }

        double sum = 0;
        for(size_t k = 1; k <= N0; k++)
        {
            sum += fMeler(x[k]);
        }
        cout << endl << "�������� ���������: " << setprecision(12) << pi * sum / N0 << endl;
    }

    return 0;
}
