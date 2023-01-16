#include <iostream>
#include <cmath>
#include "vector.h"
#include "helpers.h"

class polygon {
private:
    int n;
    vector *points = new vector[n];
    vector *intersect_start = new vector[n];
    vector *intersect_end = new vector[n];

    bool getPoints();

    vector calculate_centroid();

    double calculate_area();

    double calculate_Ix();

    double calculate_Iy();

    double find_min_y();

    double find_max_y();

    vector centroid;
    double Ix;
    double Iy;
    double area;
    double min_y;
    double max_y;

public:
    double E;
    double SIGMA_Y;
    double EPSILON_Y;
    double Y_STEP;
    double PHI_MAX;
    double PHI_STEP;

    polygon();

    void setPoint(vector, int);

    vector getCentroid();

    double calculate_moment(double, double &);
};

polygon::polygon() {
    int a;
    double b;

    cout << "Enter E Modulus: ";
    cin >> b;
    this->E = b;

    cout << "Enter Yield Stress: ";
    cin >> b;
    this->SIGMA_Y = b;

    this->EPSILON_Y = this->SIGMA_Y / this->E;

    cout << "Enter Y Increment Step: ";
    cin >> b;
    this->Y_STEP = b / 1000;

    cout << "Increase Curvature Upto: ";
    cin >> b;
    this->PHI_MAX = b;

    cout << "Enter Curvature Increment Step: ";
    cin >> b;
    this->PHI_STEP = b;

    cout << "Enter Number of Polygon Vertexes: ";
    cin >> a;
    this->n = a;

    cout << "Scale: mm" << endl << endl;
    if (getPoints()) {
        centroid = calculate_centroid();
        area = calculate_area();
        Ix = calculate_Ix();
        Iy = calculate_Iy();
        min_y = find_min_y();
        max_y = find_max_y();

        for (int i = 0; i < n - 1; i++) {
            intersect_start[i] = points[i];
            intersect_end[i] = points[i + 1];
        }
        intersect_start[n - 1] = points[n - 1];
        intersect_end[n - 1] = points[0];
    }
}

void polygon::setPoint(vector v, int i) {
    this->points[i] = v;
}

bool polygon::getPoints() {
    vector v;
    for (int i = 0; i < this->n; i++) {
        cout << "X" << i + 1 << ": ";
        cin >> v.x;
        cout << "Y" << i + 1 << ": ";
        cin >> v.y;
        cout << endl;

        v.x /= 1000;
        v.y /= 1000;

        setPoint(v, i);
    }
    return true;
}

double polygon::calculate_area() {
    double a = 0;

    for (int i = 0; i < this->n; ++i)
        a += points[i].x * points[i + 1].y - points[i + 1].x * points[i].y;

    return a / 2;
}

vector polygon::calculate_centroid() {
    vector c;
    double signedArea = 0;
    double a = 0;

    for (int i = 0; i < this->n; ++i) {
        a = points[i].x * points[i + 1].y - points[i + 1].x * points[i].y;
        signedArea += a / 2;
        c.x += (points[i].x + points[i + 1].x) * a;
        c.y += (points[i].y + points[i + 1].y) * a;
    }

    c.x /= (6 * signedArea);
    c.y /= (6 * signedArea);

    return c;
}

double polygon::calculate_Ix() {
    double Ix = 0;
    double a = 0;
    double area = this->area;

    for (int i = 0; i < this->n; ++i) {
        a = points[i].x * points[i + 1].y - points[i + 1].x * points[i].y;
        Ix += (pow(points[i].y, 2) + points[i].y * points[i + 1].y + pow(points[i + 1].y, 2)) * a / 12;
    }

    return Ix - area * pow(this->centroid.y, 2);
}

double polygon::calculate_Iy() {
    double Iy = 0;
    double a = 0;
    double area = this->area;

    for (int i = 0; i < this->n; ++i) {
        a = points[i].x * points[i + 1].y - points[i + 1].x * points[i].y;
        Iy += (pow(points[i].x, 2) + points[i].x * points[i + 1].x + pow(points[i + 1].x, 2)) * a / 12;
    }

    return Iy - area * pow(this->centroid.x, 2);
}

double polygon::find_max_y() {
    double max_y = points[0].y;
    for (int i = 1; i < this->n; i++)
        if (points[i].y > max_y)
            max_y = points[i].y;
    return max_y;
}

double polygon::find_min_y() {
    double min_y = points[0].y;
    for (int i = 1; i < this->n; i++)
        if (points[i].y < min_y)
            min_y = points[i].y;
    return min_y;
}

vector polygon::getCentroid() {
    return this->centroid;
}

double polygon::calculate_moment(double phi, double &PNA) {
    double T = 0; // Tension Force
    double C = 0; // Compression Force
    int intersection_num = 0;
    double temp_x[this->n]; // temp
    double sigma;
    double a; // dA
    double x; // intersection length in each y
    double PNA_y = PNA;
    int k; // Compare C & T
    double M = 0; // result
    double yp = EPSILON_Y / phi; // elastic core height
    double check = 0;

    while (true) {
        // Compression Force
        for (double y = this->min_y; y <= PNA_y; y += Y_STEP) {
            for (int i = 0; i < this->n; i++) {
                if ((this->intersect_start[i].y < y && this->intersect_end[i].y > y) ||
                    (this->intersect_start[i].y > y && this->intersect_end[i].y < y)) {
                    // X = ( Y - y1 ) * ( x2 - x1 ) / ( y2 - y1 ) + x1
                    temp_x[intersection_num] = (y - intersect_start[i].y) *
                                               (intersect_end[i].x - intersect_start[i].x) /
                                               (intersect_end[i].y - intersect_start[i].y) +
                                               intersect_start[i].x;
                    intersection_num++;
                }
            }

            sort(temp_x, intersection_num);

            for (int i = 0; i < intersection_num - 1; i += 2)
                x += temp_x[i + 1] - temp_x[i];

            if ((PNA_y - y) > yp) sigma = SIGMA_Y;
            else sigma = phi * (PNA_y - y) * E;

            a = x * Y_STEP;
            C += a * sigma;
            M += a * sigma * fabs(y - PNA_y);

            x = 0;
            intersection_num = 0;
        }

        // Tension Force
        for (double y = PNA_y; y <= this->max_y; y += Y_STEP) {
            for (int i = 0; i < this->n; i++) {
                if ((this->intersect_start[i].y < y && this->intersect_end[i].y > y) ||
                    (this->intersect_start[i].y > y && this->intersect_end[i].y < y)) {
                    // X = ( Y - y1 ) * ( x2 - x1 ) / ( y2 - y1 ) + x1
                    temp_x[intersection_num] = (y - intersect_start[i].y) *
                                               (intersect_end[i].x - intersect_start[i].x) /
                                               (intersect_end[i].y - intersect_start[i].y) +
                                               intersect_start[i].x;
                    intersection_num++;
                }
            }

            sort(temp_x, intersection_num);

            for (int i = 0; i < intersection_num - 1; i += 2)
                x += temp_x[i + 1] - temp_x[i];

            if ((y - PNA_y) > yp) sigma = SIGMA_Y;
            else sigma = phi * (y - PNA_y) * E;

            a = x * Y_STEP;
            T += a * sigma;
            M += a * sigma * fabs(y - PNA_y);

            x = 0;
            intersection_num = 0;
        }

        // Comparison
        if (T > C) k = 1;
        else k = -1;

        PNA_y += k * Y_STEP;

        if ((T - C) * check < 0) break;
        check = T - C;
        if (check == 0) break;

        M = 0;
        C = 0;
        T = 0;
    }

    PNA = PNA_y;
    return M;
}
