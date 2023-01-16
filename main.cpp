#include <iostream>
#include <fstream>
#include "polygon.h"

using namespace std;

int main() {
    polygon p;
    vector c = p.getCentroid();
    double PNA = c.y;
    ofstream FILE_Moment;
    ofstream FILE_PNA;

    FILE_Moment.open("export_moment.txt");
    FILE_PNA.open("export_pna.txt");

    FILE_Moment.precision(3);

    for (double i = 0; i < p.PHI_MAX; i += p.PHI_STEP) {
        FILE_Moment << fixed << p.calculate_moment(i, PNA) / 1000 << endl;
        FILE_PNA << fixed << PNA * 1000 << endl;
    }

    cout.precision(3);
    cout << fixed << "Mp: " << p.calculate_moment(1e300, PNA) / 1000 << " KN.m" << endl;
    cout << "PNA: " << PNA * 1000 << " mm";

    return 0;
}
