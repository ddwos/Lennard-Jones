#include <iostream>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

class Point
{
    // Class to define points
public:
    double x, y, z;
    Point(double a, double b, double c)
    {
        x = a;
        y = b;
        z = c;
    }
};

double distance(Point p1, Point p2)
{
    // distance between the points
    return sqrt(pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2) +
                pow((p1.z - p2.z), 2));
}

double potential(double r, int a = 6)
{
    // Lennard-Jones potential between two molecules at a distance r
    return 4 * (pow(1 / r, 2 * a) - pow(1 / r, a));
}

Point draw(double R)
{
    // drawing single point in a sphere of radius R

    double r = (double)(rand()) / ((double)(RAND_MAX / R));
    double phi = (double)(rand()) / ((double)(RAND_MAX / (2 * M_PI)));
    double theta = (double)(rand()) / ((double)(RAND_MAX / (M_PI)));
    Point p(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r *
            cos(theta));
    return p;
}

vector<Point> draw_points(int N, double R)
{
    // drawing N points in a sphere of radius R

    vector<Point> v;
    for (int i = 0; i < N; i++)
        v.push_back(draw(R));
    return v;
}

double energy(vector<Point> points, int a)
{
    // total energy for certain molecules distribution
    double result = 0;
    double r;
    int n = points.size();
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            r = distance(points.at(i), points.at(j));
            result += potential(r, a);
        }
    }
    return result;
}

double potential_derivative(double r, int a = 6)
{
    // derivative of Lennard-Jones potential between two molecules at a
    // distance r
    return 4 * (-2 * a * pow(r, -2 - 2 * a) + a * pow(r, -2 - a));
}

vector<double> grad_energy(vector<Point> points, int a = 6)
{
    // potential energy gradient, consecutive coorrdinates are partial
    // derivatives of energy in respect to x0, y0, z0, x1, y1, z1, x2 etc...
    vector<double> result;
    double r;
    double sumx, sumy, sumz;
    int n = points.size();
    for (int i = 0; i < n; i++)
    {
        sumx = 0;
        sumy = 0;
        sumz = 0;
        for (int j = 0; j < n; j++)
        {
            if (j == i)
            {
                continue;
            }

            r = distance(points.at(i), points.at(j));

            sumx += (points.at(i).x - points.at(j).x) *
                    potential_derivative(r, a);
            sumy += (points.at(i).y - points.at(j).y) *
                    potential_derivative(r, a);
            sumz += (points.at(i).z - points.at(j).z) *
                    potential_derivative(r, a);
        }
        result.push_back(sumx);
        result.push_back(sumy);
        result.push_back(sumz);
    }
    return result;
}

vector<double> change_coordinates(vector<Point> points)
{
    // function changes vector which elements are points to vector which
    // elements are coordinates
    vector<double> result;
    for (int i = 0; i < points.size(); i++)
    {
        result.push_back(points.at(i).x);
        result.push_back(points.at(i).y);
        result.push_back(points.at(i).z);
    }
    return result;
}

vector<Point> change_point(vector<double> points)
{
    // function changes vector which elements are coordinates to vector which
    // elements are points
    int n = points.size();
    Point p(0, 0, 0);
    vector<Point> result;
    for (int i = 0; i < n / 3; i++)
    {
        result.push_back(p);
    }
    int j = 0;
    for (int i = 0; i < n / 3; i++)
    {
        result.at(i).x += points.at(j);
        j++;
        result.at(i).y += points.at(j);
        j++;
        result.at(i).z += points.at(j);
        j++;
    }
    return result;
}

vector<double> operator*(double x, vector<double> &A)
{
    // vector multiplication by scalar
    vector<double> B;
    B = A;
    for (int i = 0; i < A.size(); i++)
    {
        B.at(i) = B.at(i) * x;
    }
    return B;
}

vector<double> operator+(vector<double> &A, vector<double> &B)
{
    // vector addition
    vector<double> C;
    C = B;
    for (int i = 0; i < B.size(); i++)
    {
        C.at(i) += A.at(i);
    }
    return C;
}

double norm(vector<double> vec)
{
    // norm of a vector
    int n = vec.size();
    double result = 0;
    for (int i = 0; i < n; i++)
        result += vec.at(i) * vec.at(i);
    return sqrt(result);
}

double findk(vector<double> p, vector<double> d, int w = 6)
{
    // finding the optimal parameter k for the method of conjugate gradients
    double rho = 0.001;
    vector<double> vec, pom;
    int j = 0;
    double a = pow(rho, j);
    pom = a * d;
    vec = p + pom;
    double delta = 20;
    while (energy(change_point(vec), w) > energy(change_point(p), w) - a * a *
            delta * norm(d) * norm(d))
    {
        j++;
        a = pow(rho, j);
        pom = a * d;
        vec = p + pom;
        if (j > 15)
        {
            break;
        }
    }
    return a;
}

int main()
{
    srand(time(NULL));

// method of conjugate gradients
beginning:
    int a = 1;          // potential parameter
    int N = 7;          // number of points
    int limit = 100000; // iteration limit
    vector<Point> points = draw_points(N, 1 * 1.1);
    vector<double> points1, points2, d, gradient1, gradient2, pom1, pom2;
a:
    cout << "Resetting iterations \n";
    int iterations = 0;
    points1 = change_coordinates(points);
    double B;
    gradient1 = grad_energy(change_point(points1), a);
    gradient2 = grad_energy(change_point(points2), a);
    d = (-1) * gradient1;
    double epsilon = 0.0001;
    double k = 0.01;

    while (norm(gradient1) >= epsilon)
    {
        points2 = points1;
        k = findk(points1, d, a);
        pom1 = k * d;
        points1 = points1 + pom1;
        gradient1 = grad_energy(change_point(points1), a);
        gradient2 = grad_energy(change_point(points2), a);
        B = (norm(gradient1) / norm(gradient2)) * (norm(gradient1) /
                norm(gradient2));
        pom2 = B * d;
        d = (-1) * gradient2;
        d = pom2 + d;
        iterations++;
        if (iterations % 20000 == 0)
            cout << "Iteration number: " << iterations << endl;
        if (iterations >= limit)
        {
            cout << "ITERATION NUMBER EXCEEDED"
                 << "\n";
            break;
        }
    }

    if (energy(change_point(points1), a) > -0.9)
    {
        // protection for values close to zero
        cout << "ERROR_0" << endl;
        goto beginning;
    }

    for (int i = 0; i < points1.size(); i++)
    {
        if (abs(points1.at(i)) > 10)
        {
            // protection for too large coordinates
            cout << "ERROR_1" << endl;
            goto beginning;
        }
    }

    // calculations for successive values of parameter a
    if (a < 6)
    {
        a++;
        cout << "increasing a -> " << a << endl;
        points = change_point(points1);
        goto a;
    }

    cout << "END" << endl;

    // writing out the coordinates of the points
    for (int i = 0; i < change_point(points1).size(); i++)
    {
        cout << change_point(points1).at(i).x << " " <<
             change_point(points1).at(i).y << " " <<
             change_point(points1).at(i).z << "\n";
    }

    // writing out the total energy
    cout << "energy: " << energy(change_point(points1), a) << "\n";
    return 0;
}
