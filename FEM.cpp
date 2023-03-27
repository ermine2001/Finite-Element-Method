#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

using namespace std;

struct node
{
    int nr;
    double x, y = 0.0;  // coordinates
    int BC = 0;         // boundary condition
};

/* Storage of the node numbers of each element */
struct element
{
    int ID[4] = { 0, 0, 0, 0 };
};

/* creating a structure to keep most important info about mesh; object of this structure is a singleton */
struct grid
{
    int nN, nE; // number of nodes and number of elements in mesh

} siatka; 

/* structure with information about material properties, determining the Fourier-Kirchoff equation */
struct global_data
{
    int SimulationTime;
    int SimulationStepTime;
    int Conductivity;
    int Alfa;
    int Tot;
    int InitialTemp;
    int Density;
    int SpecificHeat;

} glob;

/* It's worth noticing that every quadrangle uses the same shape functions, and also same coefficients(ksi and eta) and weights in Gauss-Legendre quadrature. 
In order to avoid writing that info many times we use one universal structure: */

struct UniversalElement
{
    double ksi2[4][4];
    double eta2[4][4];
    double ksi3[9][4];
    double eta3[9][4];
    double tabWagi1[2] = { 1,1 };
    double tabWagi2[3] = { (5.0 / 9.0), (8.0 / 9.0), (5.0 / 9.0) };
    double tabWagi3[4] = { ((18.0 - sqrt(30.0)) / 36.0), ((18.0 + sqrt(30.0)) / 36.0), ((18.0 - sqrt(30.0)) / 36.0),((18.0 + sqrt(30.0)) / 36.0) };
    double** funkcje_ksztaltu=new double *[2];

} el_uni;


/* function to return the proper shape function depending on its number */
double fk4(int nr, double ksi, double eta) {
    if (nr == 0) return (0.25 * (1 - ksi) * (1 - eta));
    if (nr == 1) return (0.25 * (1 + ksi) * (1 - eta));
    if (nr == 2) return (0.25 * (1 + ksi) * (1 + eta));
    if (nr == 3) return (0.25 * (1 - ksi) * (1 + eta));   
}

/* fk4_reverse was created to solve the problem of reversing order of points while computing values for HBC matrix and P wector on 3rd and 4th
side of the element - basically that situation happens when we must choose the opposite direction while going through the points */
double fk4_reverse(int nr, double ksi, double eta) {
    if (nr == 0) return (0.25 * (1 - ksi) * (1 + eta));
    if (nr == 1) return (0.25 * (1 + ksi) * (1 + eta));
    if (nr == 2) return (0.25 * (1 + ksi) * (1 - eta));
    if (nr == 3) return (0.25 * (1 - ksi) * (1 - eta));
}

// function to compute value of shape functions with given ksi and eta factors, with 2,3, or 4 point igtegration scheme 
void oblicz_N(int ilosc, double *N[]) { 

    int il_pc = sqrt(ilosc);
    double* ksi = new double[il_pc];
    double* eta = new double[il_pc];
    for (int i = 0; i < il_pc; i++) {
        ksi[i] = 0.0;
    }

    if (il_pc == 2) {
        double ksi2[2] = {{-(1.0 / sqrt(3.0))}, {1.0 / sqrt(3.0)}};
        double eta2[2] = { {-(1.0 / sqrt(3.0))}, {1.0 / sqrt(3.0)} };
        ksi = ksi2;
        eta = eta2;
    }
    if (il_pc == 3) {
        double ksi3[3] = { -(double)sqrt(3.0 / 5.0) , 0.0, (double)sqrt(3.0 / 5.0) };
        double eta3[3] = { -(double)sqrt(3.0 / 5.0) , 0.0, (double)sqrt(3.0 / 5.0) };
        ksi = ksi3;
        eta = eta3;
    }
    if (il_pc == 4) {
        double ksi4[4] = { (-(double)sqrt(3.0 / 7.0 + (2.0 / 7.0 * sqrt(6.0 / 5.0)))), ((double)sqrt(3.0 / 7.0 - (2.0 / 7.0 * sqrt(6.0 / 5.0)))), ((double)sqrt(3.0 / 7.0 + (2.0 / 7.0 * sqrt(6.0 / 5.0)))), (-(double)sqrt(3.0 / 7.0 - (2.0 / 7.0 * sqrt(6.0 / 5.0)))) };
        double eta4[4] = { (-(double)sqrt(3.0 / 7.0 + (2.0 / 7.0 * sqrt(6.0 / 5.0)))), ((double)sqrt(3.0 / 7.0 - (2.0 / 7.0 * sqrt(6.0 / 5.0)))), ((double)sqrt(3.0 / 7.0 + (2.0 / 7.0 * sqrt(6.0 / 5.0)))), (-(double)sqrt(3.0 / 7.0 - (2.0 / 7.0 * sqrt(6.0 / 5.0)))) };
        ksi = ksi4;
        eta = eta4;
    }
    int pom2 = 0;

    for (int i = 0; i < il_pc; i++) {
        for (int j = 0; j < il_pc; j++) {

          N[pom2][0] = 0.25 * (1 - ksi[i]) * (1 - eta[j]);
          N[pom2][1] = 0.25 * (1 + ksi[i]) * (1 - eta[j]);
          N[pom2][2] = 0.25 * (1 + ksi[i]) * (1 + eta[j]);
          N[pom2][3] = 0.25 * (1 - ksi[i]) * (1 + eta[j]);
          pom2++;
        }
    }
}

// Function providing Gauss elimination algorithm, neccessary to compute the system of equations 
bool gauss(int n, double** AB, double* X)
{
    const double eps = 1e-12; //The constant approximation of zero
    int i, j, k;
    double m, s;
    double pom[16][17] = {};

    // elimination of coefficients:
    for (i = 0; i < n - 1; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            if (fabs(AB[i][i]) < eps) { cout << "ERR NO 1"; return false; }
            m = -AB[j][i] / AB[i][i];
            for (k = i + 1; k <= n; k++)
                AB[j][k] += m * AB[i][k];
        }
    }

    // solving for unknowns.
    for (i = n - 1; i >= 0; i--)
    {
        s = AB[i][n];
        for (j = n - 1; j >= i + 1; j--)
            s -=AB[i][j] * X[j];
        if (fabs(AB[i][i]) < eps) { cout << "ERR NO 2"; return false; }
        X[i] = s / AB[i][i];
    }
    double min = X[0];
    double max = X[0];
   
    for (int i = 0; i <= n - 1; i++) {
        if (X[i] >= max) max = X[i];
        if (X[i] <= min) min = X[i];
    }
    cout << "   Temp Min: " << min;
    cout << "   Temp Max: " << max;
   
    return 0;
}


int main() {
    
   

    double w = 0.25;
    double eta_pom[4] = { -(double)(1.0 / sqrt(3.0)), -(double)(1.0 / sqrt(3.0)), (double)(1.0 / sqrt(3.0)), (double)(1.0 / sqrt(3.0)) };
    double ksi_pom[4] = { -(double)(1.0 / sqrt(3.0)), (double)(1.0 / sqrt(3.0)), -(double)(1.0 / sqrt(3.0)),(double)(1.0 / sqrt(3.0)) };
    double eta_pom3[9] = { -(double(sqrt(3.0 / 5.0))), -(double(sqrt(3.0 / 5.0))), -(double(sqrt(3.0 / 5.0))), 0, 0, 0, (double(sqrt(3.0 / 5.0))) , (double(sqrt(3.0 / 5.0))) , (double(sqrt(3.0 / 5.0))) };
    double ksi_pom3[9] = { -(double(sqrt(3.0 / 5.0))), 0, (double(sqrt(3.0 / 5.0))), -(double(sqrt(3.0 / 5.0))) , 0, (double(sqrt(3.0 / 5.0))) , -(double(sqrt(3.0 / 5.0))) ,0, (double(sqrt(3.0 / 5.0))) };

    double ksi = 0.0;
    double eta = 0.0;


    double deriv_ksi[4][4]{ w * (1 - eta), -w * (1 - eta), w * (1 + eta), -w * (1 + eta) };    // Derivatives with respect to ksi
    double deriv_eta[4][4] = { -w * (1 - ksi), -w * (1 + ksi), w * (1 + ksi), w * (1 - ksi) }; // Derivatives with respect to eta

    cout << "2-point scheme. \n Ksi derivatives: " << endl;
    for (int i = 0; i < 4; i++)
    {


        eta = eta_pom[i];
        
        deriv_ksi[i][0] = -w * (1 - eta); 
        el_uni.ksi2[i][0] = -w * (1 - eta);
        cout << deriv_ksi[i][0] << " ";

        deriv_ksi[i][1] = w * (1 - eta);
        el_uni.ksi2[i][1] = w * (1 - eta);
        cout << deriv_ksi[i][1] << " ";

        deriv_ksi[i][2] = w * (1 + eta);
        el_uni.ksi2[i][2] = w * (1 + eta);
        cout << deriv_ksi[i][2] << " ";

        deriv_ksi[i][3] = -w * (1 + eta);
        el_uni.ksi2[i][3] = -w * (1 + eta);
        cout << deriv_ksi[i][3] << " ";
        cout << endl;
    }
    cout << "Eta derivatives :" << endl;
    for (int i = 0; i < 4; i++)
    {
        ksi = ksi_pom[i];
        
        deriv_eta[i][0] = -w * (1 - ksi);
        el_uni.eta2[i][0] = -w * (1 - ksi);
        cout << deriv_eta[i][0] << " ";

        deriv_eta[i][1] = -w * (1 + ksi);
        el_uni.eta2[i][1] = -w * (1 + ksi);
        cout << deriv_eta[i][1] << " ";

        deriv_eta[i][2] = w * (1 + ksi);
        el_uni.eta2[i][2] = w * (1 + ksi);
        cout << deriv_eta[i][2] << " ";

        deriv_eta[i][3] = w * (1 - ksi);
        el_uni.eta2[i][3] = w * (1 - ksi);
        cout << deriv_eta[i][3] << " ";
        cout << endl;

    }

    cout << " \n 3-point scheme. Array with ksi derivatives: " << endl;
    for (int i = 0; i < 9; i++)
    {
        eta = eta_pom3[i];
        el_uni.ksi3[i][0] = -w * (1 - eta);
        cout << el_uni.ksi3[i][0] << " ";

        el_uni.ksi3[i][1] = w * (1 - eta);
        cout << el_uni.ksi3[i][1] << " ";


        el_uni.ksi3[i][2] = w * (1 + eta);
        cout << el_uni.ksi3[i][2] << " ";


        el_uni.ksi3[i][3] = -w * (1 + eta);
        cout << el_uni.ksi3[i][3] << " ";
        cout << endl;
    }

    //eta
    cout << " Eta derivatives array: " << endl;

    for (int i = 0; i < 9; i++)
    {
        ksi = ksi_pom3[i];
        el_uni.eta3[i][0] = -w * (1 - ksi);
        cout << el_uni.eta3[i][0] << " ";
        el_uni.eta3[i][1] = -w * (1 + ksi);
        cout << el_uni.eta3[i][1] << " ";
        el_uni.eta3[i][2] = w * (1 + ksi);
        cout << el_uni.eta3[i][2] << " ";
        el_uni.eta3[i][3] = w * (1 - ksi);
        cout << el_uni.eta3[i][3] << " ";
        cout << endl;
    }
    cout << endl;
    cout << " helper tables for ksi and eta " << endl;
    cout << "The 2-point integration scheme: " << endl;
    cout << " ksi: " << endl;


    for (int i = 0; i < 4; i++)
    {
        cout << i + 1 << " " << ksi_pom[i] << endl;
    }
    cout << " eta: " << endl;
    for (int i = 0; i < 4; i++)
    {
        cout << i + 1 << " " << eta_pom[i] << endl;
    }
    cout << endl << "The 3-point integration scheme " << endl;
    cout << " ksi: " << endl;
    for (int i = 0; i < 9; i++)
    {
        cout << i + 1 << " " << ksi_pom3[i] << endl;
    }
    cout << " eta: " << endl;
    for (int i = 0; i < 9; i++)
    {
        cout << i + 1 << " " << eta_pom3[i] << endl;
    }

    /* Reading data from a file */
    cout.precision(9);
    string line;
    string pom;
    ifstream file;
    file.open("grid1.txt", ios::in);

    if (file.good() == false) cout << " FILE ERROR " << endl;


    string pom22[30]; // strings from getline
    int line_number;
    int pozycja_start;
    int pozycja_end;
    int ilosc_cyfr;
    int liczba_pom;
    int pom1[11]; //strings converted to numbers
    string test;
    string a, nr, x, y, BC;
    int node_nr = 1;

    // global_data and grid: lines 1-10
    for (int i = 1; i < 11; i++)
    {
        getline(file, line);
        test = line;
        pozycja_start = test.find_first_of("0123456789");
        pozycja_end = test.find_last_of("0123456789");
        ilosc_cyfr = pozycja_end - pozycja_start + 1;
        pom22[i - 1] = test;
        test = test.substr(pozycja_start, ilosc_cyfr + 1);
        liczba_pom = stoi(test);
        pom1[i - 1] = liczba_pom;
    }
    siatka.nN = pom1[8];
    siatka.nE = pom1[9]; 

    cout << " DEBUG " << "\nsiatka.nN: " << siatka.nN << "\nsiatka.nE:  " << siatka.nE << endl;

    do
    {
        string b;
        getline(file, line);
        b = line;

    } while (0);

    //storage of nodes and elements in vectors

    vector <node> pomND(siatka.nN + 2);
    vector <element> pomEL(siatka.nE + 2);

    cout << "Grid: " << endl;

    for (int i = 0; i < (siatka.nN); i++) // from line 12: nodes

    {
        getline(file, line);
        a = line;

        int start_nr, start_x, start_y;
        int end_nr, end_x, end_y;
        int il_nr, il_x, il_y;


        start_nr = a.find_first_of("123456789");
        end_nr = (a.find_first_of(",")) - 1;
        il_nr = end_nr - start_nr + 1;

        start_x = (a.find_first_of(".")) - 1; ///
        end_x = (a.find_last_of(",")) - 1;
        il_x = end_x - start_x + 1;

        start_y = (a.find_last_of(",")) + 2;
        end_y = (a.find_last_of("0123456789"));
        il_y = end_y - start_x + 1;
        x = a.substr(start_x, il_x);
        y = a.substr(start_y, il_y);
        nr = a.substr(start_nr, il_nr);
        pomND[node_nr - 1].x = stod(x);
        pomND[node_nr - 1].y = stod(y);
        pomND[node_nr - 1].nr = stoi(nr);
        cout << "Node no. " << pomND[node_nr - 1].nr << " coordinates:   x: \t" << pomND[node_nr - 1].x << "     y:\t" << pomND[node_nr - 1].y << " flag: " << pomND[node_nr - 1].BC << endl;
        node_nr++;

    }

    do
    {
        string b;
        getline(file, line);
        b = line;

    } while (0);
    cout << "New lines: " << endl;
    // reading elements
    string word;
    int i = 0;
    for (i = 0; i < siatka.nE; i++)
    {

        file >> word;
        file >> word;

        pomEL[i].ID[0] = stoi(word);
        cout << pomEL[i].ID[0] << " ";

        file >> word;
        pomEL[i].ID[1] = stoi(word);
        cout << pomEL[i].ID[1] << " ";

        file >> word;
        pomEL[i].ID[2] = stoi(word);
        cout << pomEL[i].ID[2] << " ";

        file >> word;
        pomEL[i].ID[3] = stoi(word);
        cout << pomEL[i].ID[3] << " ";

        cout << endl;
    }
    for (int h = 0; h < siatka.nE; h++)
    {


        cout << pomEL[h].ID[0] << " ";
        cout << pomEL[h].ID[1] << " ";
        cout << pomEL[h].ID[2] << " ";
        cout << pomEL[h].ID[3] << " ";
        cout << endl;
    }
    file >> word; // saving *BC
    string nowy;
    nowy = word;
    cout << nowy << endl;
    for (int h = 0; h < siatka.nN; h++)

        while (!file.eof()) {

            file >> word;
            int nowy2;

            for (int i = 0; i < siatka.nN; i++)
            {
                nowy2 = stoi(word);

                if (nowy2 == pomND[i].nr)
                {

                    pomND[i].BC += 1;
                    cout << "nr " << i + 1 << ": " << pomND[i].BC << endl;
                    file >> word;

                }
                else
                {
                    int roznica;
                    roznica = nowy2 - pomND[i].nr;
                    for (int h = 0; h < (roznica - 1); h++)
                    {
                        i++;
                    }
                }
            }

        }
    for (int h = 0; h < siatka.nN; h++)
    {
        cout << " flag of node no. " << pomND[h].nr << ":  " << pomND[h].BC << endl;
    }

    file.close();

   /* end of file */

    cout << "\n lines 1-8:  " << endl;
    cout << "glob.SimulationTime: " << (glob.SimulationTime = pom1[0]) << endl;
    cout << "glob.SimulationStepTime: " << (glob.SimulationStepTime = pom1[1]) << endl;
    cout << "glob.Conductivity: " << (glob.Conductivity = pom1[2]) << endl;
    cout << "glob.alfa: " << (glob.Alfa = pom1[3]) << endl;
    cout << "glob.Tot: " << (glob.Tot = pom1[4]) << endl;
    cout << "glob.InitialTemp: " << (glob.InitialTemp = pom1[5]) << endl;
    cout << "glob.Density: " << (glob.Density = pom1[6]) << endl;
    cout << "glob.SpecificHeat: " << (glob.SpecificHeat = pom1[7]) << endl;
    cout << "\n lines 9 i 10:  " << endl;
    cout << "siatka.nN: " << (siatka.nN = pom1[8]) << endl;
    cout << "siatka.nE: " << (siatka.nE = pom1[9]) << endl;

    cout << "---------Number of nodes: " << siatka.nN << endl;
    cout << "--------Number of elements: " << siatka.nE << endl;


    //  *Matrix H* , *Matrix C*

    int pc_H = 4;// point scheme
     
    el_uni.funkcje_ksztaltu = new double* [pc_H];
    for (int i = 0; i < pc_H; i++)
    {
        el_uni.funkcje_ksztaltu[i] = new double[4];
    }

    for (int i = 0; i < pc_H; i++) {
        for (int j = 0; j < 4; j++) {
            el_uni.funkcje_ksztaltu[i][j] = 0.0;
        
        }
    }

    oblicz_N(pc_H, el_uni.funkcje_ksztaltu);
    cout << " ~ SHAPE FUNCTIONS VALUES ~\n";
    for (int i = 0; i < pc_H; i++) {
        for (int j = 0; j < 4; j++) {
            cout << el_uni.funkcje_ksztaltu[i][j] << " ";

        } 
        cout << endl;
    }
    double* detJ = new double[pc_H]; // Determinants of matrices formed by 4 points, index 0: 1 point, etc.
    for (int i = 0; i < pc_H; i++) {
        detJ[i] = 0.0;
    }
    double** m_j = new double* [pc_H]; // array of Jacobian matrices
    // generating dimensions 
    for (int i = 0; i < pc_H; i++)
        m_j[i] = new double[4];
    for (int i = 0; i < pc_H; i++) {
        for (int j = 0; j < 4; j++) {
            m_j[i][j] = 0.0; 
        }
    }
    double** DNDX = new double* [pc_H];
    double** DNDY = new double* [pc_H];
    for (int i = 0; i < pc_H; i++)
    {
        DNDX[i] = new double[4];
        DNDY[i] = new double[4];
    }
    for (int i = 0; i < pc_H; i++) {
        for (int j = 0; j < 4; j++) {
            DNDX[i][j] = 0.0;
            DNDY[i][j] = 0.0;
        }
    }
   double  eta_pom2p[4] = { -(double)(1.0 / sqrt(3.0)), -(double)(1.0 / sqrt(3.0)), (double)(1.0 / sqrt(3.0)), (double)(1.0 / sqrt(3.0)) };
   double  ksi_pom2p[4] = { -(double)(1.0 / sqrt(3.0)), (double)(1.0 / sqrt(3.0)), -(double)(1.0 / sqrt(3.0)),(double)(1.0 / sqrt(3.0)) };
   double  eta_pom3p[9] = { -(double(sqrt(3.0 / 5.0))), -(double(sqrt(3.0 / 5.0))), -(double(sqrt(3.0 / 5.0))), 0, 0, 0, (double(sqrt(3.0 / 5.0))) , (double(sqrt(3.0 / 5.0))) , (double(sqrt(3.0 / 5.0))) };
   double  ksi_pom3p[9] = { -(double(sqrt(3.0 / 5.0))), 0, (double(sqrt(3.0 / 5.0))), -(double(sqrt(3.0 / 5.0))) , 0, (double(sqrt(3.0 / 5.0))) , -(double(sqrt(3.0 / 5.0))) ,0, (double(sqrt(3.0 / 5.0))) };
  // 
   double** C_wszystkiepc = new double* [pc_H * 4];



    double** HX_wszystkie_pc = new double* [pc_H * 4];
    double** HY_wszystkie_pc = new double* [pc_H * 4];// For storing all X and Y matrices of points (FOR ONE ELEMENT!)
    double** H_X_Y_pc = new double* [pc_H * 4]; // The sum of hx and hy, multiplied by detJ and cond
    double** C_pc = new double* [pc_H * 4]; //  detJ and cond - multiplying

    for (int i = 0; i < (pc_H * 4); i++)
    {
        HX_wszystkie_pc[i] = new double[4];
        HY_wszystkie_pc[i] = new double[4];
        H_X_Y_pc[i] = new double[4];
        C_wszystkiepc[i] = new double[4];
        C_pc[i] = new double[4];

    }
    for (int i = 0; i < (pc_H * 4); i++) {
        for (int j = 0; j < 4; j++) {
            HX_wszystkie_pc[i][j] = 0.0;
            HY_wszystkie_pc[i][j] = 0.0;
            H_X_Y_pc[i][j] = 0.0;
            C_wszystkiepc[i][j] = 0.0;
            C_pc[i][j] = 0.0;

        }
    }
    double tabx[4] = {}; // x from the grid; a new iteration of the loop over elements = new x and y.
    double taby[4] = {}; // y - from grid;
    double row_map[4][4] = { {}, {}, {}, {} };
    double col_map[4][4] = { {}, {}, {}, {} };
    double matrixH_X[4][4] = { {}, {}, {}, {} };
    double matrixH_Y[4][4] = { {}, {}, {}, {} };
    double matrix_C[4][4] = { {}, {}, {}, {} };

    double H_element[4][4] = {};
    double C_element[4][4] = {};
    double C_wlasciwe_elementu[4][4] = {};
   // tablica_H_G - the global H matrix, which is created by aggregation, with a size of NxN(siatka.nN x siatka.nN).
    double** tablica_H_G = new double* [siatka.nN];
    double** tablica_HiHBC_G = new double* [siatka.nN];
    double** tablica_HiC50= new double* [siatka.nN]; // H with boundary condition + C/time
    double** tablica_C50 = new double* [siatka.nN];

    double** tablica_C = new double* [siatka.nN];
    for (int i = 0; i < siatka.nN; i++)
    {
        tablica_H_G[i] = new double[siatka.nN];
        tablica_HiHBC_G[i] = new double[siatka.nN];
        tablica_C[i] = new double[siatka.nN];
        tablica_C50[i] = new double[siatka.nN];
        tablica_HiC50[i] = new double[siatka.nN];
    }
    for (int i = 0; i < siatka.nN; i++) {
        for (int j = 0; j < siatka.nN; j++) {
            tablica_H_G[i][j] = 0.0;
            tablica_HiHBC_G[i][j] = 0.0;
            tablica_C[i][j] = 0.0;
            tablica_C50[i][j] = 0.0;
            tablica_HiC50[i][j] = 0.0;
        }
    }
            double** H_local = new double* [4 * siatka.nE];
            double** HBC_local = new double* [4 * siatka.nE];
            double** C_local = new double* [4 * siatka.nE];

            for (int i = 0; i < 4 * siatka.nE; i++) {

                H_local[i] = new double[4];
                HBC_local[i] = new double[4];
                C_local[i] = new double[4];
            }
            for (int i = 0; i < 4 * siatka.nE; i++) {
                for (int j = 0; j < 4; j++) {
                    H_local[i][j] = 0.0;
                    H_local[i][j] = 0.0;
                    C_local[i][j] = 0.0;
                }
            }

            /*  LOOP OVER ALL ELEMENTS OF THE GRID */
            for (int E = 0; E < siatka.nE; E++) { 

                 //Loading nodes belonging to elements:

                for (int i = 0; i < 4; i++) {

                    int tmp;
                    tmp = pomEL[E].ID[i];
                    tabx[i] = pomND[tmp - 1].x;
                    cout << "tabx: " << tabx[i] << " ";
                    taby[i] = pomND[tmp - 1].y;
                    cout << " taby: " << taby[i] << endl;
                    int kas = 0;
                }
                /*  Loop over integration points (pc_H) */
                for (int j = 0; j < pc_H; j++) {

                    //Calculating the Jacobian matrices for the integration points.
                    for (int i = 0; i < 4; i++) {
                        if (pc_H == 4) {
                            m_j[j][0] = m_j[j][0] + el_uni.ksi2[j][i] * tabx[i];   //dx/dksi
                            m_j[j][1] = m_j[j][1] + el_uni.ksi2[j][i] * taby[i];   //dy/dksi
                            m_j[j][2] = m_j[j][2] + el_uni.eta2[j][i] * tabx[i];   //dx/deta
                            m_j[j][3] = m_j[j][3] + el_uni.eta2[j][i] * taby[i];   //dy/deta
                        }
                        else if (pc_H == 9) {
                            m_j[j][0] = m_j[j][0] + el_uni.ksi3[j][i] * tabx[i];   //dx/dksi
                            m_j[j][1] = m_j[j][1] + el_uni.ksi3[j][i] * taby[i];   //dy/dksi
                            m_j[j][2] = m_j[j][2] + el_uni.eta3[j][i] * tabx[i];   //dx/deta
                            m_j[j][3] = m_j[j][3] + el_uni.eta3[j][i] * taby[i];   //dy/deta
                        }

                    }


                    cout << "\n Jacobians matrices for 4 integration points: " << endl;


                    for (int i = 0; i < 4; i++) {
                        cout << m_j[j][i] << "  ";
                    }

                    detJ[j] = m_j[j][0] * m_j[j][3] - m_j[j][1] * m_j[j][2];
                    cout << "determinant det J of " << j + 1 << " integration point: " << detJ[j] << endl;

                    double pom = 0.0;
                    pom = m_j[j][0];
                    m_j[j][0] = m_j[j][3];
                    m_j[j][3] = pom;
                    m_j[j][1] *= (double)(-1.0);
                    m_j[j][2] *= (double)(-1.0);


                    cout << "\n   Jacobian matrix of " << j + 1 << " integration point; with swapping along the diagonal and reversing signs, etc: " << endl;

                    for (int i = 0; i < 4; i++) {
                        cout << m_j[j][i] << "  ";
                    }


                    m_j[j][0] *= (double)(1.0) / detJ[j];
                    m_j[j][3] *= (double)(1.0) / detJ[j];
                    m_j[j][1] *= (double)(1.0) / detJ[j];
                    m_j[j][2] *= (double)(1.0) / detJ[j];


                    cout << "\n Jacobian matrices * detJ: \n: " << endl;

                    for (int i = 0; i < 4; i++) {
                        cout << m_j[j][i] << " ";
                    }


                    for (int i = 0; i < 4; i++) {

                        if (pc_H == 4) {
                            DNDX[j][i] = m_j[j][0] * el_uni.ksi2[j][i] + m_j[j][1] * el_uni.eta2[j][i];
                            DNDY[j][i] = m_j[j][2] * el_uni.ksi2[j][i] + m_j[j][3] * el_uni.eta2[j][i];

                            
                        }
                        if (pc_H == 9) {
                            DNDX[j][i] = m_j[j][0] * el_uni.ksi3[j][i] + m_j[j][1] * el_uni.eta3[j][i];
                            DNDY[j][i] = m_j[j][2] * el_uni.ksi3[j][i] + m_j[j][3] * el_uni.eta3[j][i];
                            
                        }

                    }

                    for (int i = 0; i < 4; i++) {
                        matrixH_X[i][0] = DNDX[j][i] * DNDX[j][0];
                        matrixH_X[i][1] = DNDX[j][i] * DNDX[j][1];
                        matrixH_X[i][2] = DNDX[j][i] * DNDX[j][2];
                        matrixH_X[i][3] = DNDX[j][i] * DNDX[j][3];

                        matrix_C[i][0] = el_uni.funkcje_ksztaltu[j][i] * el_uni.funkcje_ksztaltu[j][0];
                        matrix_C[i][1] = el_uni.funkcje_ksztaltu[j][i] * el_uni.funkcje_ksztaltu[j][3];
                        matrix_C[i][2] = el_uni.funkcje_ksztaltu[j][i] * el_uni.funkcje_ksztaltu[j][2];
                        matrix_C[i][3] = el_uni.funkcje_ksztaltu[j][i] * el_uni.funkcje_ksztaltu[j][1];


                        HX_wszystkie_pc[i + 4 * j][0] = matrixH_X[i][0];
                        HX_wszystkie_pc[i + 4 * j][1] = matrixH_X[i][1];
                        HX_wszystkie_pc[i + 4 * j][2] = matrixH_X[i][2];
                        HX_wszystkie_pc[i + 4 * j][3] = matrixH_X[i][3];

                        C_wszystkiepc[i + 4 * j][0] = matrix_C[i][0];
                        C_wszystkiepc[i + 4 * j][1] = matrix_C[i][1];
                        C_wszystkiepc[i + 4 * j][2] = matrix_C[i][2];
                        C_wszystkiepc[i + 4 * j][3] = matrix_C[i][3];

                    }
                   
                    //The next iteration of the integration point loop creates a matrix for the next point. To avoid overwriting the values that we still need, they are copied to other arrays

                    for (int i = 0; i < 4; i++) {
                        matrixH_Y[i][0] = DNDY[j][i] * DNDY[j][0];
                        matrixH_Y[i][1] = DNDY[j][i] * DNDY[j][1];
                        matrixH_Y[i][2] = DNDY[j][i] * DNDY[j][2];
                        matrixH_Y[i][3] = DNDY[j][i] * DNDY[j][3];

                        HY_wszystkie_pc[i + 4 * j][0] = matrixH_Y[i][0];
                        HY_wszystkie_pc[i + 4 * j][1] = matrixH_Y[i][1];
                        HY_wszystkie_pc[i + 4 * j][2] = matrixH_Y[i][2];
                        HY_wszystkie_pc[i + 4 * j][3] = matrixH_Y[i][3];
                    }

                    for (int i = 0; i < 4; i++) {
                       
                        H_X_Y_pc[i + 4 * j][0] = (HX_wszystkie_pc[i + 4 * j][0] + HY_wszystkie_pc[i + 4 * j][0]) * glob.Conductivity * detJ[j];
                        H_X_Y_pc[i + 4 * j][1] = (HX_wszystkie_pc[i + 4 * j][1] + HY_wszystkie_pc[i + 4 * j][1]) * glob.Conductivity * detJ[j];
                        H_X_Y_pc[i + 4 * j][2] = (HX_wszystkie_pc[i + 4 * j][2] + HY_wszystkie_pc[i + 4 * j][2]) * glob.Conductivity * detJ[j];
                        H_X_Y_pc[i + 4 * j][3] = (HX_wszystkie_pc[i + 4 * j][3] + HY_wszystkie_pc[i + 4 * j][3]) * glob.Conductivity * detJ[j];
                        C_pc[i + 4 * j][0] = C_wszystkiepc[i + 4 * j][0] * glob.SpecificHeat * glob.Density * detJ[j];
                        C_pc[i + 4 * j][1] = C_wszystkiepc[i + 4 * j][1] * glob.SpecificHeat * glob.Density * detJ[j];
                        C_pc[i + 4 * j][2] = C_wszystkiepc[i + 4 * j][2] * glob.SpecificHeat * glob.Density * detJ[j];
                        C_pc[i + 4 * j][3] = C_wszystkiepc[i + 4 * j][3] * glob.SpecificHeat * glob.Density * detJ[j];

                    }



                    for (int i = 0; i < 4; i++) {
                        H_element[i][0] += H_X_Y_pc[i + 4 * j][0];
                        H_element[i][1] += H_X_Y_pc[i + 4 * j][1];
                        H_element[i][2] += H_X_Y_pc[i + 4 * j][2];
                        H_element[i][3] += H_X_Y_pc[i + 4 * j][3];
                        //matrix C
                        C_element[i][0] += C_pc[i + 4 * j][0];
                        C_element[i][1] += C_pc[i + 4 * j][1];
                        C_element[i][2] += C_pc[i + 4 * j][2];
                        C_element[i][3] += C_pc[i + 4 * j][3];

                    }
                    // Swap row 2 with row 4
                    for (int i = 0; i < 4; i++) {
                        for (int k = 0; k < 4; k++){
                            if (i == 1) { C_wlasciwe_elementu[i][k] = C_element[3][k]; }
                             if (i == 3) { C_wlasciwe_elementu[i][k] = C_element[1][k]; }
                           if(i!=1&&i!=3) { C_wlasciwe_elementu[i][k] = C_element[i][k]; }

                        }
                    }



                } 

                //End of loop over integration points

                cout << " H_X_Y_PC: \n";
                for (int i = 0; i < pc_H; i++) {
                    cout << "pcH = " << i << endl;
                    for (int j = 0; j < 4; j++) {
                        cout << H_X_Y_pc[j + 4 * i][0] << " ";
                        cout << H_X_Y_pc[j + 4 * i][1] << " ";
                        cout << H_X_Y_pc[j + 4 * i][2] << " ";
                        cout << H_X_Y_pc[j + 4 * i][3] << " ";
                        cout << endl;

                    }
                    cout << endl;
                }

               //  Creating the final local H matrices for each element
                for (int j = 0; j < 4; j++) {   //rows 
                    for (int i = 0; i < 4; i++) //columns
                    {

                        H_local[j + 4 * E][i] = H_element[j][i]; // Transformation of 4x4 to 36x4(referring to rows)
                        C_local[j + 4 * E][i] = C_wlasciwe_elementu[j][i];

                    }
                }
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        H_element[i][j] = 0.0;
                        C_element[i][j] = 0.0;
                        C_wlasciwe_elementu[i][j] = 0.0;
                    }
                }

                cout << endl;
                // ----------------------* MATRIX H AGGREGATION *. --------------------------------

                for (int j = 0; j < 4; j++) {
                    for (int i = 0; i < 4; i++) {
                        row_map[j][i] = pomEL[E].ID[j];
                        col_map[j][i] = pomEL[E].ID[i];

                    }
                }
                cout << "Mapping of rows and columns of H matrix \n";
                for (int j = 0; j < 4; j++) {
                    for (int i = 0; i < 4; i++) {
                        cout << row_map[j][i] << ",";
                        cout << col_map[j][i] << ",";
                        cout << "  ";

                    }cout << endl;
                }
                //system("PAUSE");


                for (int h = 0; h < 4; h++)
                {
                    for (int f = 0; f < 4; f++) {
                        int w = row_map[h][f];
                        int k = col_map[h][f];
                        tablica_H_G[w - 1][k - 1] += H_local[h + 4 * E][f];

                        // * C matrix aggregation*

                        tablica_C[w - 1][k - 1] += C_local[h + 4 * E][f];

                    }
                }
                //----end of matrix H aggregation---------------------------------------------

                 
                for (int j = 0; j < pc_H; j++)
                {
                    for (int i = 0; i < 4; i++) {
                        m_j[j][i] = 0.0;
                    }
                }




            }//the end of loop over elements E


            cout << "H local for all elements\n";
            for (int e = 0; e < 9; e++) {
                cout << "H element no: " << e + 1 << endl;
                for (int j = 0; j < 4; j++) {
                    for (int i = 0; i < 4; i++) {
                        cout << H_local[j + 4 * e][i] << " ";
                    }
                    cout << endl;
                }
                cout << endl;
            }
            cout << " C_local for all elements \n";
            for (int e = 0; e < 9; e++) {
                cout << "C element no: " << e + 1 << endl;
                for (int j = 0; j < 4; j++) {
                    for (int i = 0; i < 4; i++) {
                        cout << C_local[j + 4 * e][i] << " ";
                    }
                    cout << endl;
                }
                cout << endl;
            }


            // aggregation result: GLOBAL MATRIX H      

            //wypisywanie
            cout << "\n global matrix H, without boundary condition (tablica_H_G): \n";
            for (int i = 0; i < siatka.nN; i++) {
                for (int j = 0; j < siatka.nN; j++) {
                    cout << tablica_H_G[i][j] << " ";
                }
                cout << endl;
            }
            system("PAUSE");

            cout << "\n global matrix C: \n";
            for (int i = 0; i < siatka.nN; i++) {
                for (int j = 0; j < siatka.nN; j++) {
                    cout << tablica_C[i][j] << " ";
                }
                cout << endl;
            }
            system("PAUSE");
            // Matrix C with SimulationStepTime
            for (int i = 0; i < siatka.nN; i++) {
                for (int j = 0; j < siatka.nN; j++) {
                    tablica_C50[i][j] = tablica_C[i][j] / glob.SimulationStepTime;
                }
            }
        


            for (int i = 0; i < 4; i++)
            {
                delete[] HX_wszystkie_pc[i];
                delete[] HY_wszystkie_pc[i];
                delete[] H_X_Y_pc[i];
                delete[] m_j[i];
            }

            delete[] HX_wszystkie_pc;
            delete[] HY_wszystkie_pc;
            delete[] H_X_Y_pc;
            delete[] m_j;


            delete[] detJ;
            //   -------------------MATRIX HBC and VECTOR P ---------------
                   // 1. Loop over elements
                   // 2. In each element: check 4 edges, if there is a flag 1 on both nodes creating the edge - create HBC matrices in each integration point;
                   // 3. Sum HBC from each integration point for the edge
                   // 4. Sum all edges for the element - HBC matrix for the element [4E][4]
                   // 5. Filling the HBC_local matrix [4 number of elements][4]

double ksi_eta2[2] = { -(1.0 / sqrt(3.0)), 1.0 / sqrt(3.0) };
double ksi_eta3[3] = { -(double)sqrt(3.0 / 5.0) , 0.0, (double)sqrt(3.0 / 5.0) };
double ksi_eta4[4] = { (-(double)sqrt(3.0 / 7.0 + (2.0 / 7.0 * sqrt(6.0 / 5.0)))), ((double)sqrt(3.0 / 7.0 - (2.0 / 7.0 * sqrt(6.0 / 5.0)))), ((double)sqrt(3.0 / 7.0 + (2.0 / 7.0 * sqrt(6.0 / 5.0)))), (-(double)sqrt(3.0 / 7.0 - (2.0 / 7.0 * sqrt(6.0 / 5.0)))) };
ksi = 0.0;
eta = 0.0;
double x_[4] = { 0.0, 0.025, 0.025 , 0.0 }; 
double y_[4] = { 0.0, 0.0, 0.025, 0.025 }; 
int pcHBC;

cout << "\n Please provide the number of integration points for HBC matrix and P vector (choose 2, 3, or 4):\n";
cin >> pcHBC;
//helper arrays for HBC 
double** WektorP_local= new double* [siatka.nE];
for (int i = 0; i < siatka.nE; i++) {
    WektorP_local[i] = new double[4];

}
double** tabD = new double* [pcHBC];
double** tabP = new double* [pcHBC];
double** tabG = new double* [pcHBC];
double** tabL = new double* [pcHBC];
//helper arrays for P vector 
double** wektorD = new double* [pcHBC];
double** wektorPR = new double* [pcHBC]; // right wall
double** wektorG = new double* [pcHBC];
double** wektorL = new double* [pcHBC];
double wektorD_all[4] = {};
double wektorPR_all[4] = {};
double wektorG_all[4] = {};
double wektorL_all[4] = {};
double wektorP_elementu[4] = {};
for (int i = 0; i < pcHBC; i++) {  // generating
    tabD[i] = new double[4];
    tabP[i] = new double[4];
    tabG[i] = new double[4];
    tabL[i] = new double[4];
    wektorD[i] = new double[4];
    wektorPR[i] = new double[4];
    wektorG[i] = new double[4];
    wektorL[i] = new double[4];
}
for (int i = 0; i < pcHBC; i++) { 
    for (int j = 0; j < 4; j++) {
        tabD[i][j] = 0.0;
        tabP[i][j] = 0.0;
        tabG[i][j] = 0.0;
        tabL[i][j] = 0.0;
    }
}
double tabD_all[4][4] = {};
double  tabG_all[4][4] = {};
double tabP_all[4][4] = {};
double  tabL_all[4][4] = {};
double det_J = 0.0; // the same value for HBC matrix and P vector
double HBC_el[4][4] = {};



//  -------- hbc - loop over elements ---------

for (int e = 0; e < siatka.nE; e++) {
    cout << "~~~~~~~~~~~~ ELEMENT NO. " << e + 1 << " ~~~~~~~~~~~~ \n";

    x_[0] = pomND[pomEL[e].ID[0] - 1].x;   y_[0] = pomND[pomEL[e].ID[0] - 1].y;
    x_[1] = pomND[pomEL[e].ID[1] - 1].x;   y_[1] = pomND[pomEL[e].ID[1] - 1].y;
    x_[2] = pomND[pomEL[e].ID[2] - 1].x;   y_[2] = pomND[pomEL[e].ID[2] - 1].y;
    x_[3] = pomND[pomEL[e].ID[3] - 1].x;   y_[3] = pomND[pomEL[e].ID[3] - 1].y;
    for (int i = 0; i < 4; i++) {
        cout << " node no " << i + 1 << ", x:" << x_[i] << ", y: " << y_[i] << endl;
    }
    for (int i = 0; i < pcHBC; i++) { // loop over integration points

        if ((pomND[pomEL[e].ID[0] - 1].BC == 1) && (pomND[pomEL[e].ID[1] - 1].BC == 1)) {
            // lower edge
            det_J = sqrt(pow(pomND[pomEL[e].ID[0] - 1].x - pomND[pomEL[e].ID[1] - 1].x, 2) + pow(pomND[pomEL[e].ID[0] - 1].y - pomND[pomEL[e].ID[1] - 1].y, 2)) / 2.0;
            if (pcHBC == 2) { ksi = ksi_eta2[i]; }
            if (pcHBC == 3) { ksi = ksi_eta3[i]; }
            if (pcHBC == 4) { ksi = ksi_eta4[i]; }
            eta = -1.0;

            for (int j = 0; j < 4; j++) {
                tabD[i][j] = fk4(j, ksi, eta);
                wektorD[i][j] = fk4(j, ksi, eta);

            }

            cout << "tabD[i][0] : " << tabD[i][0];

            // Hbc - Summing up the whole edge (from all integration points)
            for (int k = 0; k < 4; k++) { 
                for (int j = 0; j < 4; j++) {
                    if (pcHBC == 2) {
                        tabD_all[k][j] += tabD[i][k] * tabD[i][j] * el_uni.tabWagi1[i] * glob.Alfa * det_J;
                    }
                    if (pcHBC == 3) {
                        tabD_all[k][j] += tabD[i][k] * tabD[i][j] * el_uni.tabWagi2[i] * glob.Alfa * det_J;
                    }
                    if (pcHBC == 4) {
                        tabD_all[k][j] += tabD[i][k] * tabD[i][j] * el_uni.tabWagi3[i] * glob.Alfa * det_J;
                    }
                }
            }

            // P vector - summing up the lower edge (from all integration points)
            for (int k = 0; k < 4; k++) {
                if (pcHBC == 2) {
                    wektorD_all[k] += wektorD[i][k] * glob.Tot * glob.Alfa * det_J * el_uni.tabWagi1[i];
                }
                if (pcHBC == 3) {
                    wektorD_all[k] += wektorD[i][k] * glob.Tot * glob.Alfa * det_J * el_uni.tabWagi2[i];
                }
                if (pcHBC == 4) {
                    wektorD_all[k] += wektorD[i][k] * glob.Tot * glob.Alfa * det_J * el_uni.tabWagi3[i];
                }

            }



                    }
                    // right edge
                    if ((pomND[pomEL[e].ID[1] - 1].BC) * (pomND[pomEL[e].ID[2] - 1].BC)) {
                        det_J = sqrt(pow(pomND[pomEL[e].ID[1] - 1].x - pomND[pomEL[e].ID[2] - 1].x, 2) + pow(pomND[pomEL[e].ID[1] - 1].y - pomND[pomEL[e].ID[2] - 1].y, 2)) / 2.0;

                        ksi = 1.0;
                        if (pcHBC == 2) { eta = ksi_eta2[i]; }
                        if (pcHBC == 3) { eta = ksi_eta3[i]; }
                        if (pcHBC == 4) { eta = ksi_eta4[i]; }

                        for (int j = 0; j < 4; j++) {
                            tabP[i][j] = fk4(j, ksi, eta);
                            wektorPR[i][j] = fk4(j, ksi, eta);

                        }
                        
                        for (int k = 0; k < 4; k++) {
                            for (int j = 0; j < 4; j++) {
                                if (pcHBC == 2) {
                                    tabP_all[k][j] += tabP[i][k] * tabP[i][j] * el_uni.tabWagi1[i] * glob.Alfa * det_J;
                                }
                                if (pcHBC == 3) {
                                    tabP_all[k][j] += tabP[i][k] * tabP[i][j] * el_uni.tabWagi2[i] * glob.Alfa * det_J;
                                }
                                if (pcHBC == 4) {
                                    tabP_all[k][j] += tabP[i][k] * tabP[i][j] * el_uni.tabWagi3[i] * glob.Alfa * det_J;
                                }
                            }
                        }
                        // P vector - summing up the right edge (from all integration points)
                        for (int k = 0; k < 4; k++) {
                            if (pcHBC == 2) {
                                wektorPR_all[k] += wektorPR[i][k] * glob.Tot * glob.Alfa * det_J * el_uni.tabWagi1[i];
                            }
                            if (pcHBC == 3) {
                                wektorPR_all[k] += wektorPR[i][k] * glob.Tot * det_J * glob.Alfa * el_uni.tabWagi2[i];
                            }
                            if (pcHBC == 4) {
                                wektorPR_all[k] += wektorPR[i][k] * glob.Tot * det_J * glob.Alfa * el_uni.tabWagi3[i];
                            }

                        }
                    }
                    // upper edge - reversed order of shape functions
                    if (pomND[pomEL[e].ID[2] - 1].BC == 1 && pomND[pomEL[e].ID[3] - 1].BC == 1) {
                        det_J = sqrt(pow(pomND[pomEL[e].ID[2] - 1].x - pomND[pomEL[e].ID[3] - 1].x, 2) + pow(pomND[pomEL[e].ID[2] - 1].y - pomND[pomEL[e].ID[3] - 1].y, 2)) / 2.0;
                        if (pcHBC == 2) { ksi = ksi_eta2[i]; }
                        if (pcHBC == 3) { ksi = ksi_eta3[i]; }
                        if (pcHBC == 4) { ksi = ksi_eta4[i]; }

                        eta = 1.0;
                       
                        tabG[i][0] = fk4_reverse(3, ksi, eta);
                        tabG[i][1] = fk4_reverse(2, ksi, eta);
                        tabG[i][2] = fk4_reverse(1, ksi, eta);
                        tabG[i][3] = fk4_reverse(0, ksi, eta);
                        
                        wektorG[i][0] = fk4_reverse(3, ksi, eta);
                        wektorG[i][1] = fk4_reverse(2, ksi, eta);
                        wektorG[i][2] = fk4_reverse(1, ksi, eta);
                        wektorG[i][3] = fk4_reverse(0, ksi, eta);
                        for (int k = 0; k < 4; k++) {
                            for (int j = 0; j < 4; j++) {
                                if (pcHBC == 2) {
                                    tabG_all[k][j] += tabG[i][k] * tabG[i][j] * el_uni.tabWagi1[i] * glob.Alfa * det_J;
                                }
                                if (pcHBC == 3) {
                                    tabG_all[k][j] += tabG[i][k] * tabG[i][j] * el_uni.tabWagi2[i] * glob.Alfa * det_J;  
                                }
                                if (pcHBC == 4) {
                                    tabG_all[k][j] += tabG[i][k] * tabG[i][j] * el_uni.tabWagi3[i] * glob.Alfa * det_J;
                                  }
                            }
                        }
                       // P vector - summing up the upper edge (from all integration points)
                        for (int k = 0; k < 4; k++) {
                            if (pcHBC == 2) {
                                wektorG_all[k] += wektorG[i][k] * glob.Tot * glob.Alfa * det_J * el_uni.tabWagi1[i];
                            }
                            if (pcHBC == 3) {
                                wektorG_all[k] += wektorG[i][k] * glob.Tot * glob.Alfa * det_J * el_uni.tabWagi2[i];
                            }
                            if (pcHBC == 4) {
                                wektorG_all[k] += wektorG[i][k] * glob.Tot * glob.Alfa * det_J * el_uni.tabWagi3[i];
                            }

                        }
                    }
                    //left edge - reversed order of shape functions
                    if (pomND[pomEL[e].ID[3] - 1].BC == 1 && pomND[pomEL[e].ID[0] - 1].BC == 1) {
                        det_J = sqrt(pow(pomND[pomEL[e].ID[3] - 1].x - pomND[pomEL[e].ID[0] - 1].x, 2) + pow(pomND[pomEL[e].ID[3] - 1].y - pomND[pomEL[e].ID[0] - 1].y, 2)) / 2.0;

                        ksi = -1.0;
                        if (pcHBC == 2) { eta = ksi_eta2[i]; }
                        if (pcHBC == 3) { eta = ksi_eta3[i]; }
                        if (pcHBC == 4) { eta = ksi_eta4[i]; }
                        
                        tabL[i][0] = fk4_reverse(3, ksi, eta);
                        tabL[i][1] = fk4_reverse(2, ksi, eta);
                        tabL[i][2] = fk4_reverse(1, ksi, eta);
                        tabL[i][3] = fk4_reverse(0, ksi, eta);
                        wektorL[i][0] = fk4_reverse(3, ksi, eta);
                        wektorL[i][1] = fk4_reverse(2, ksi, eta);
                        wektorL[i][2] = fk4_reverse(1, ksi, eta);
                        wektorL[i][3] = fk4_reverse(0, ksi, eta);
                        for (int k = 0; k < 4; k++) {
                            for (int j = 0; j < 4; j++) {
                                if (pcHBC == 2) {
                                    tabL_all[k][j] += tabL[i][k] * tabL[i][j] * el_uni.tabWagi1[i] * glob.Alfa * det_J;
                                   }
                                if (pcHBC == 3) {
                                    tabL_all[k][j] += tabL[i][k] * tabL[i][j] * el_uni.tabWagi2[i] * glob.Alfa * det_J;  
                                }
                                if (pcHBC == 4) {
                                    tabL_all[k][j] += tabL[i][k] * tabL[i][j] * el_uni.tabWagi3[i] * glob.Alfa * det_J;
                                }
                            }
                        }
                        // P vector - summing up the left edge (from all integration points)
                        for (int k = 0; k < 4; k++) {
                            if (pcHBC == 2) {
                                wektorL_all[k] += wektorL[i][k] * glob.Tot * glob.Alfa * det_J * el_uni.tabWagi1[i];
                            }
                            if (pcHBC == 3) {
                                wektorL_all[k] += wektorL[i][k] * glob.Tot * glob.Alfa * det_J * el_uni.tabWagi2[i];
                            }
                            if (pcHBC == 4) {
                                wektorL_all[k] += wektorL[i][k] * glob.Tot * glob.Alfa * det_J * el_uni.tabWagi3[i];
                            }

                        }
                    }
                    det_J = 0.0;
                }
                // the end of loop over integration points

                cout << "End of loop over integration points. Collecting 4 edge matrices of the element into one. \n ";
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        HBC_el[i][j] = tabD_all[i][j] + tabP_all[i][j] + tabG_all[i][j] + tabL_all[i][j];
                        cout << HBC_el[i][j] << " ";
                    }cout << endl;
                }
                cout << "Collecting 4 edge vectors of the element into one array: element vector[4]: \n ";
                for (int i = 0; i < 4; i++) {
                    wektorP_elementu[i] = wektorD_all[i] + wektorPR_all[i] + wektorG_all[i] + wektorL_all[i];
                }
                cout << "\n  HBC of the lower edge: \n";
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        cout << tabD_all[i][j] << " ";
                    }
                    cout << endl;
                }
                cout << "\n  HBC of the right edge: \n";
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        cout << tabP_all[i][j] << " ";
                    }
                    cout << endl;
                }
                cout << "\n  HBC of the upper edge: \n";
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        cout << tabG_all[i][j] << " ";
                    }
                    cout << endl;
                }
                cout << "\n  HBC of the left edge: \n";
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        cout << tabL_all[i][j] << " ";
                    }
                    cout << endl;
                }
                // Entering HBC of each element into one array. Each element has one 4x4 matrix regardless of the number of integration points 
                for (int j = 0; j < 4; j++) {//wierwsze 
                    for (int i = 0; i < 4; i++) //kolumny
                    {
                        HBC_local[j + 4 * e][i] = HBC_el[j][i];
                    }
                }
                // Entering the vector of each element [4] into one array (wektorP_local[liczba el][4]
                for (int i = 0; i < 4; i++) {
                    WektorP_local[e][i] = wektorP_elementu[i];

                }
                // Clearing arrays after processing one element
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        tabD_all[i][j] = 0.0;
                        tabP_all[i][j] = 0.0;
                        tabG_all[i][j] = 0.0;
                        tabL_all[i][j] = 0.0;
                        HBC_el[i][j] = 0.0;
                    }
                    cout << endl;
                }
                // clearing vector P helper arrays 
                for (int i = 0; i < 4; i++) {
                    wektorD_all[i] = 0.0;
                    wektorPR_all[i] = 0.0;
                    wektorG_all[i] = 0.0;
                    wektorL_all[i] = 0.0;
                    wektorP_elementu[i] = 0.0;

                }



            } //~~~~~~~~~~~~~~`the end of loop over elements in HBC 

            cout << " --- HBC ---- " << endl;
            for (int e = 0; e < siatka.nE; e++) {
                cout << "HBC element no: " << e + 1 << endl;
                for (int j = 0; j < 4; j++) {
                    for (int i = 0; i < 4; i++) {
                        cout << HBC_local[j + 4 * e][i] << " ";
                    }
                    cout << endl;
                }
                cout << endl;
            }

            cout <<"\n----, WektorP_local: ----- \n";
             for (int e = 0; e < siatka.nE; e++) {
                 cout << "p vector, el no: " << e + 1 << endl;
                
                     for (int i = 0; i < 4; i++) {
                         cout << WektorP_local[e][i] << " ";
                     }
                     cout << endl << endl;
               
             }
             system("PAUSE");

            // adding HBC local to H local     
            for (int el = 0; el < siatka.nE; el++) {
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        H_local[i + 4 * el][j] = H_local[i + 4 * el][j] + HBC_local[i + 4 * el][j];
                    }
                }
            }
            //---------------------------------aggregation of H+HBC. Loop over elements ----------------------------
            for (int E = 0; E < siatka.nE; E++) {
                for (int j = 0; j < 4; j++) {
                    for (int i = 0; i < 4; i++) {
                        row_map[j][i] = pomEL[E].ID[j];
                        col_map[j][i] = pomEL[E].ID[i];

                    }
                }
                for (int h = 0; h < 4; h++)
                {
                    for (int f = 0; f < 4; f++) {
                        int w = row_map[h][f];
                        int k = col_map[h][f];
                        tablica_HiHBC_G[w - 1][k - 1] += H_local[h + 4 * E][f];


                    }
                }
            }
         // ----------------the end of aggregation of H + HBC -----------------------------


            for (int e = 0; e < siatka.nE; e++) {
                cout << "H+HBC element no: " << e + 1 << endl;
                for (int j = 0; j < 4; j++) {
                    for (int i = 0; i < 4; i++) {
                        cout << H_local[j + 4 * e][i] << " ";
                    }
                    cout << endl;
                }
                cout << endl;
            }

            // printing
            cout << "\n global matrix H+HBC (tablica_HiHBC_G)\n";
            for (int i = 0; i < siatka.nN; i++) {
                for (int j = 0; j < siatka.nN; j++) {
                    cout << tablica_HiHBC_G[i][j] << " ";
                }
                cout << endl;
            }
            system("PAUSE");

            // summing H+hbc and C/50
            cout << "H+ C/simulation step time: \n";
            for (int i = 0; i < siatka.nN; i++) {
                for (int j = 0; j < siatka.nN; j++) {
                    tablica_HiC50[i][j] = tablica_C50[i][j] + tablica_HiHBC_G[i][j];
                    cout << tablica_HiC50[i][j] << " ";
                }
                cout << endl;
            }
            system("PAUSE");
            // --------------------- aggregation of vector P  ---------------------
            double map_vecP[4];
            double P_global[16] = {}; // vector P aggregated
            for (int E = 0; E < siatka.nE; E++) {
                    for (int i = 0; i < 4; i++) {
                      
                        map_vecP[i] = pomEL[E].ID[i];

                    }
                
                    for (int f = 0; f < 4; f++) {
                        int P_pom = map_vecP[f];
                        P_global[P_pom-1] += WektorP_local[E][f];


                    }
                
            }
            cout << "\n Aggregated P vector (P_global): \n";
            for (int i = 0; i < siatka.nN; i++) {
                cout <<i << ", "<< P_global[i] << endl;
            }
            system("PAUSE");

              
         
            // arrays of the extended matrix for Gaussian elimination.

            double** matrix_gauss = new double* [siatka.nN];
            for (int i = 0; i < siatka.nN; i++) {
                matrix_gauss[i] = new double [siatka.nN + 1];

            }
            // GAUSS for steady state 
            cout << "MATRIX GAUSS \n";
            // Entering into this array columns 1-16 with H matrix with boundary conditions and column 17 with P vector (aggregated)
            for (int i = 0; i < siatka.nN; i++) {
                for (int j = 0; j < siatka.nN + 1; j++) {
                    if (j != 16)  matrix_gauss[i][j] = tablica_HiHBC_G[i][j];
                    else matrix_gauss[i][j] = P_global[i];
                    cout << matrix_gauss[i][j] << " ";
                }
                cout << endl;
            }
            system("PAUSE");
            

            double *temperatures = new double [siatka.nN];
            for (int i = 0; i < siatka.nN; i++) {
                temperatures[i] = 0.0;
            }

            // Performing Gauss elimination for steady state
            gauss(16, matrix_gauss, temperatures);

            cout << "\n temperatures for steady state \n";
            for (int i = 0; i < siatka.nN; i++) {
                cout << temperatures[i] << endl;
            }
            system("PAUSE");


            // -------------------- simulation. Unsteady state solution ---------------------------------
            double* P_sym = new double[siatka.nN]; // the entire simulation P vector (second part of the equation after the minus sign)
            double *t0= new double[siatka.nN]; // Initial temperature vector 
            for (int i = 0; i < 16; i++) {
                t0[i] = glob.InitialTemp; // temp = 100,  will change in each iteration
                P_sym[i] = 0.0;
            }

            double t_pom[16] = {};
            
            system("PAUSE");
            // creating an extended matrix H+C/dt + vector P in column 17 
            for (int i = 0; i < siatka.nN; i++) {
                for (int j = 0; j < siatka.nN + 1; j++) {
                    if (j != 16)  matrix_gauss[i][j] = tablica_HiC50[i][j]; 
                    else matrix_gauss[i][j] = t_pom[i];
                }
            } 
            for (int i = 0; i < 16; i++) {
                temperatures[i] = 0.0;
            }
           
           
            system("PAUSE");
            for (int i = 0; i < 16; i++) {
                temperatures[i] = 0.0;
            }
           


            //---------------------The simulation of changes over time -  loop.----------------------------------
            cout << "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
            cout << "\n  Computed temperatures at mesh nodes. Unsteady state. \n";
            for (int i = glob.SimulationStepTime; i <= glob.SimulationTime; i = i + glob.SimulationStepTime) {
                double tabC_pom[16][16];
                // prawa strona - duzy wektor 
                cout << " \n T = " << i << ":  ";
                    for (int k = 0; k < siatka.nN; k++) {
                        for (int j = 0; j < siatka.nN; j++) {
                            tabC_pom[k][j] = tablica_C50[k][j]; // helper array for C
                            tablica_C50[k][j] *= t0[j];

                            t_pom[k] += tablica_C50[k][j];
                            tablica_C50[k][j] = tabC_pom[k][j];

                        }
                        t_pom[k] += P_global[k];
                        
                    }
                    // left side of equation (before minus sign)
                    
                    for (int k = 0; k < siatka.nN; k++) {
                        for (int j = 0; j < siatka.nN + 1; j++) {
                            matrix_gauss[k][j] = 0.0;
                            if (j != 16)  matrix_gauss[k][j] = tablica_HiC50[k][j]; 
                            else matrix_gauss[k][j] = t_pom[k];
                           
                        }
                    }
                
                for (int k = 0; k < 16; k++) {
                    temperatures[k] = 0.0;
                }
                gauss(16, matrix_gauss, temperatures);

                for (int k = 0; k < 16; k++) {
                    t0[k] = temperatures[k];
                    t_pom[k] = 0.0;
                }
            }
            

            return 0;

        }

                     

