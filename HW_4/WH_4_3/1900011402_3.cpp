#include <fstream>
using namespace std;

int main()
{
    const double SIGMA = 10;
    const double RHO = 28;
    const double BETA = 8.0 / 3.0;
    double t_ini = 0;
    double t_final = 10.0;
    int stepnum = 100000;
    double y1_ini = 12;
    double y2_ini = 4;
    double y3_ini = 0;

    //准备把数据存入文件
    ofstream fout;
    fout.open("Lorenz_attractor_data.txt");
    if (!fout)
        return -1;
    fout.clear();

    double steplength = (t_final - t_ini) / stepnum;
    double y1 = y1_ini, y2 = y2_ini, y3 = y2_ini;
    fout << 0 << ' ' << y1 << ' ' << y2 << ' ' << y3 << endl;
    for (int i = 1; i <= stepnum; i++)
    {
        double temp1 = y1, temp2 = y2, temp3 = y3;
        y1 += steplength * (temp2 * temp3 - BETA * temp1);
        y2 += steplength * SIGMA * (temp3 - temp2);
        y3 += steplength * ((RHO - temp1) * temp2 - temp3);
        fout << i * steplength << ' ' << y1 << ' ' << y2 << ' ' << y3 << endl;
    }

    fout.close();
    return 0;
}