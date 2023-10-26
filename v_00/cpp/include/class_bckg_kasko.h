#ifndef CLASS_BCKG_KASKO_H
#define CLASS_BCKG_KASKO_H

#include "class_bckg.h"

class BCKG_KASKO : public BACKGROUND
{

private:
    double  a1, a2;
    double  b0, phi;
    double  f, k, n;
    double  rho_b, p_b, gdfd;
    double  Te, Ti, Ti2Te;
    double  Ve, Vi, Vc, Vb;

public:

    /// structure, storing background parameters
    struct BCKGR_PARAM
    {
        double  a1, a2;
        double  b0, phi;
        double  f, k, n;
        double  rho_b, p_b, gdfd;
        double  Te, Ti, Ti2Te;
        double  Ve, Vi, Vc, Vb;
    }   BckgParams;

    BCKGR_PARAM * pBckg;                /// pointer to background parameters structure

    /// destructor
    ~BCKG_KASKO(void) {}

    /// constructor
    BCKG_KASKO(void) : BACKGROUND()
    {
        pBckg = &BckgParams;            /// pointer to structure, containing parameters
        read_bckg_config(  InOutFnames[MODEL].c_str()  );
    }


    /// constructor
    BCKG_KASKO(INI_BCKG_SETT* pIniSett) : BACKGROUND(pIniSett)
    {
        pBckg = &BckgParams;            /// pointer to structure, containing parameters
        read_bckg_config(  InOutFnames[MODEL].c_str()  );
    }


    /// read background parameters
    void read_bckg_config(const char* infile)
    {
        std::cout << "configuration input file: " << infile << "\n";

        char* tmp = NULL;
        int   out = 0;

        FILE *fid;
        if( (fid = fopen(infile, "r")) == NULL )
        {
            out = printf("\nError: file %s was not opened\a\a\a\n", infile);
            exit(0);
        }

        char tline[256];

        /// skip header
        for (int j = 1; j<15; j++){ tmp = fgets(tline, 256, fid); }

        /// read <a1>, <a2>
        out = fscanf(fid, "%lf%lf\n", &a1, &a2);        //std::cout << "a1 = " << a1 << ",\ta2 = " << a2 << "\n";

        /// skip comments
        for (int j = 1; j<4; j++) { tmp = fgets(tline, 256, fid); }

        /// read <b0>, <phi>
        out = fscanf(fid, "%lf%lf\n", &b0, &phi);       //std::cout << "b0 = " << b0 << ",\tphi = " << phi << "\n";

        /// skip comments
        for (int j = 1; j<4; j++) { tmp = fgets(tline, 256, fid); }

        /// read <f>, <k>, <n>
        out = fscanf(fid, "%lf%lf%lf\n", &f, &k, &n);   //std::cout << "f = " << f << ",\tk = " << k << ,\tn = " << n << "\n";

        /// skip comments
        for (int j = 1; j<4; j++) { tmp = fgets(tline, 256, fid); }

        /// read <rho_b>, <p_b>
        out = fscanf(fid, "%lf%lf\n", &rho_b, &p_b);    //std::cout << "rho_b = " << rho_b << ",\tp_b = " << p_b << "\n";

        /// skip comments
        for (int j = 1; j<4; j++) { tmp = fgets(tline, 256, fid); }

        /// read <gdfd>
        out = fscanf(fid, "%lf\n", &gdfd);              //std::cout << "gdfd = " << gdfd << "\n";

        /// skip comments
        for (int j = 1; j<4; j++) { tmp = fgets(tline, 256, fid); }

        /// read <Ti>
        out = fscanf(fid, "%lf\n", &Ti);                //std::cout << "Ti = " << Ti << "\n";
        tmp = fgets(tline, 256, fid);

        /// skip comments
        for (int j = 1; j<4; j++) { tmp = fgets(tline, 256, fid); }

        /// read <Vi>
        out = fscanf(fid, "%lf\n", &Vi);                //std::cout << "Vi = " << Vi << "\n";
        fclose(fid);

        /// check temperature
        if (  ( Ti < 0.0 ) || ( Ti > 0.5 )  )
        {
            out = printf("\nError! Illegal value of Ti = %lf. Legal range is [0, 0.5].\n\n",Ti);
            exit(0);
        }

        /// evaluate velocities
        double* Vel = new double[4];
        Te = 0.5 - Ti;
        Ti2Te = eval_velocity_harris_family(Ti,Te,Vi, Vel);
        Ve = Vel[0];    Vi = Vel[1];    Vc = Vel[2];    Vb = Vel[3];
        delete Vel;

        /// assign values to the structure fields
        double pi = com::Pi<double>();      /// pi constant

        phi = phi/180.0*pi;                 /// cast phi to radians

        (pBckg->a1) = a1;       (pBckg->a2)  = a2;
        (pBckg->b0) = b0;       (pBckg->phi) = phi;
        (pBckg->f ) = f;        (pBckg->k  ) = k;
        (pBckg->p_b) = p_b;     (pBckg->n  ) = n;
        (pBckg->gdfd) = gdfd;
        (pBckg->rho_b) = rho_b;
        (pBckg->Te) = Te;       (pBckg->Ti) = Ti;       (pBckg->Ti2Te) = Ti2Te;
        (pBckg->Ve) = Ve;       (pBckg->Vi) = Vi;
        (pBckg->Vc) = Vc;       (pBckg->Vb) = Vb;

        /// self-control output
        std::cout.setf( std::ios::showpos ); std::cout.setf( std::ios::fixed, std::ios::floatfield );
        std::cout << "\nDimensionless Background Parameters:\n\n";
        std::cout << "a1\t\t"    << (pBckg->a1)    << "\n";
        std::cout << "a2\t\t"    << (pBckg->a2)    << "\n";
        std::cout << "b0\t\t"    << (pBckg->b0)    << "\n";
        std::cout << "phi\t\t"   << (pBckg->phi)   << " rad\n";
        std::cout << "f\t\t"     << (pBckg->f)     << "\n";
        std::cout << "k\t\t"     << (pBckg->k)     << "\n";
        std::cout << "n\t\t"     << (pBckg->n)     << "\n";
        std::cout << "rho_b\t\t" << (pBckg->rho_b) << "\n";
        std::cout << "p_b\t\t"   << (pBckg->p_b)   << "\n";
        std::cout << "gdfd\t\t"  << (pBckg->gdfd)  << "\n\n";
        std::cout << "Ti\t\t"    << (pBckg->Ti)    << "\n";
        std::cout << "Te\t\t"    << (pBckg->Te)    << "\n";
        std::cout << "Ti2Te\t\t" << (pBckg->Ti2Te) << "\n\n";
        std::cout << "Ve\t\t"    << (pBckg->Ve)    << "\n";
        std::cout << "Vi\t\t"    << (pBckg->Vi)    << "\n";
        std::cout << "Vc\t\t"    << (pBckg->Vc)    << "\n";
        std::cout << "Vb\t\t"    << (pBckg->Vb)    << "\n\n";
    } // end of the function <bckg_read_config>


    /// calculate background plasma quantities
    void calc_bckg_qnt(size_t dim0, size_t dim1, size_t dim2, double* X1, double* X2, double*** QNT)
    {
        double x, z;                        //std::cout << dim0 << ", " << dim1 << ", " << dim2 << "\n";

        for (size_t i1 = 0; i1 < dim1; i1++)
        {
            x = X1[i1] / (pNC->L);

            for (size_t i2 = 0; i2 < dim2; i2++)
            {
                z = X2[i2] / (pNC->L);       //std::cout << "x = " << x << ", z = " << z <<"\n";

                /// 1. Psi:
                if ( dim0 > 0){
                    QNT[0][i1][i2] = (log((cosh(sin(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*sin(phi-k*atan((a2+z)/(a1+x)))*pow(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z,k*(-1.0/2.0)))*sqrt(f*f+1.0)+f*cos(cos(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*cos(phi-k*atan((a2+z)/(a1+x)))*pow(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z,k*(-1.0/2.0))))*1.0/sqrt((n*n)*pow(x*x+z*z,n-1.0)+(b0*b0)*(k*k)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),-k-1.0)+b0*k*n*cos(phi-atan(z/x)*(n-1.0)-atan((a2+z)/(a1+x))*(k+1.0))*pow(x*x+z*z,n*(1.0/2.0)-1.0/2.0)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0/2.0)*2.0)));
                    if ( dim0 == 1 ) {   QNT[0][i1][i2] = (pNC->Psi0)*QNT[0][i1][i2];   }
                }

                /// 2. Rho
                if ( dim0 > 1){
                    QNT[1][i1][i2] = (pNC->Rho0)*(rho_b + exp(-2.0*QNT[0][i1][i2]));
                    if ( dim0 == 2 ) {   QNT[0][i1][i2] = (pNC->Psi0)*QNT[0][i1][i2];   }
                }

                /// 3. Pg
                if ( dim0 > 2){
                    QNT[2][i1][i2] = (pNC->P0)*(p_b + 0.5*exp(-2.0*QNT[0][i1][i2]));
                    QNT[0][i1][i2] = (pNC->Psi0)*QNT[0][i1][i2];
                }

                /// 4. Bx
                if ( dim0 > 3){
                    QNT[3][i1][i2] = (pNC->B0)*(+(((f*sin(cos(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*cos(phi-k*atan((a2+z)/(a1+x)))*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)))*(n*z*cos(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0)-1.0)-n*x*sin(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0)-1.0)+b0*k*cos(phi-k*atan((a2+z)/(a1+x)))*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0)*(a2*2.0+z*2.0)*(1.0/2.0)-(b0*k*sin(phi-k*atan((a2+z)/(a1+x)))*(a1+x)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)))/(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z))-sinh(sin(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*sin(phi-k*atan((a2+z)/(a1+x)))*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)))*sqrt(f*f+1.0)*(n*x*cos(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0)-1.0)+n*z*sin(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0)-1.0)+b0*k*sin(phi-k*atan((a2+z)/(a1+x)))*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0)*(a2*2.0+z*2.0)*(1.0/2.0)+(b0*k*cos(phi-k*atan((a2+z)/(a1+x)))*(a1+x)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)))/(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z)))*1.0/sqrt((n*n)*pow(x*x+z*z,n-1.0)+(b0*b0)*(k*k)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),-k-1.0)+b0*k*n*cos(phi-atan(z/x)*(n-1.0)-atan((a2+z)/(a1+x))*(k+1.0))*pow(x*x+z*z,n*(1.0/2.0)-1.0/2.0)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0/2.0)*2.0)+(cosh(sin(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*sin(phi-k*atan((a2+z)/(a1+x)))*pow(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z,k*(-1.0/2.0)))*sqrt(f*f+1.0)+f*cos(cos(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*cos(phi-k*atan((a2+z)/(a1+x)))*pow(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z,k*(-1.0/2.0))))*1.0/pow((n*n)*pow(x*x+z*z,n-1.0)+(b0*b0)*(k*k)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),-k-1.0)+b0*k*n*cos(phi-atan(z/x)*(n-1.0)-atan((a2+z)/(a1+x))*(k+1.0))*pow(x*x+z*z,n*(1.0/2.0)-1.0/2.0)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0/2.0)*2.0,3.0/2.0)*((n*n)*z*pow(x*x+z*z,n-2.0)*(n*2.0-2.0)-(b0*b0)*(k*k)*(k*2.0+2.0)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),-k-2.0)*(a2*2.0+z*2.0)*(1.0/2.0)+b0*k*n*sin(phi-atan(z/x)*(n-1.0)-atan((a2+z)/(a1+x))*(k+1.0))*pow(x*x+z*z,n*(1.0/2.0)-1.0/2.0)*((x*(n-1.0))/(x*x+z*z)+((a1+x)*(k+1.0))/(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z))*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0/2.0)*2.0-b0*k*n*cos(phi-atan(z/x)*(n-1.0)-atan((a2+z)/(a1+x))*(k+1.0))*pow(x*x+z*z,n*(1.0/2.0)-1.0/2.0)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-3.0/2.0)*(k+1.0)*(a2*2.0+z*2.0)+b0*k*n*z*cos(phi-atan(z/x)*(n-1.0)-atan((a2+z)/(a1+x))*(k+1.0))*pow(x*x+z*z,n*(1.0/2.0)-3.0/2.0)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0/2.0)*(n-1.0)*2.0)*(1.0/2.0))*sqrt((n*n)*pow(x*x+z*z,n-1.0)+(b0*b0)*(k*k)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),-k-1.0)+b0*k*n*cos(phi-atan(z/x)*(n-1.0)-atan((a2+z)/(a1+x))*(k+1.0))*pow(x*x+z*z,n*(1.0/2.0)-1.0/2.0)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0/2.0)*2.0))/(cosh(sin(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*sin(phi-k*atan((a2+z)/(a1+x)))*pow(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z,k*(-1.0/2.0)))*sqrt(f*f+1.0)+f*cos(cos(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*cos(phi-k*atan((a2+z)/(a1+x)))*pow(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z,k*(-1.0/2.0)))));
                }

                /// 5. By
                if ( dim0 > 4){
                    QNT[4][i1][i2] = (pNC->B0)*gdfd;
                }

                /// 6. Bz
                if ( dim0 > 5){
                    QNT[5][i1][i2] = (pNC->B0)*(-(((f*sin(cos(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*cos(phi-k*atan((a2+z)/(a1+x)))*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)))*(n*x*cos(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0)-1.0)+n*z*sin(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0)-1.0)+b0*k*cos(phi-k*atan((a2+z)/(a1+x)))*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0)*(a1*2.0+x*2.0)*(1.0/2.0)+(b0*k*sin(phi-k*atan((a2+z)/(a1+x)))*(a2+z)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)))/(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z))+sinh(sin(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*sin(phi-k*atan((a2+z)/(a1+x)))*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)))*sqrt(f*f+1.0)*(n*z*cos(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0)-1.0)-n*x*sin(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0)-1.0)-b0*k*sin(phi-k*atan((a2+z)/(a1+x)))*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0)*(a1*2.0+x*2.0)*(1.0/2.0)+(b0*k*cos(phi-k*atan((a2+z)/(a1+x)))*(a2+z)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)))/(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z)))*1.0/sqrt((n*n)*pow(x*x+z*z,n-1.0)+(b0*b0)*(k*k)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),-k-1.0)+b0*k*n*cos(phi-atan(z/x)*(n-1.0)-atan((a2+z)/(a1+x))*(k+1.0))*pow(x*x+z*z,n*(1.0/2.0)-1.0/2.0)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0/2.0)*2.0)-(cosh(sin(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*sin(phi-k*atan((a2+z)/(a1+x)))*pow(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z,k*(-1.0/2.0)))*sqrt(f*f+1.0)+f*cos(cos(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*cos(phi-k*atan((a2+z)/(a1+x)))*pow(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z,k*(-1.0/2.0))))*1.0/pow((n*n)*pow(x*x+z*z,n-1.0)+(b0*b0)*(k*k)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),-k-1.0)+b0*k*n*cos(phi-atan(z/x)*(n-1.0)-atan((a2+z)/(a1+x))*(k+1.0))*pow(x*x+z*z,n*(1.0/2.0)-1.0/2.0)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0/2.0)*2.0,3.0/2.0)*(-(n*n)*x*pow(x*x+z*z,n-2.0)*(n*2.0-2.0)+(b0*b0)*(k*k)*(k*2.0+2.0)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),-k-2.0)*(a1*2.0+x*2.0)*(1.0/2.0)+b0*k*n*sin(phi-atan(z/x)*(n-1.0)-atan((a2+z)/(a1+x))*(k+1.0))*pow(x*x+z*z,n*(1.0/2.0)-1.0/2.0)*((z*(n-1.0))/(x*x+z*z)+((a2+z)*(k+1.0))/(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z))*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0/2.0)*2.0+b0*k*n*cos(phi-atan(z/x)*(n-1.0)-atan((a2+z)/(a1+x))*(k+1.0))*pow(x*x+z*z,n*(1.0/2.0)-1.0/2.0)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-3.0/2.0)*(k+1.0)*(a1*2.0+x*2.0)-b0*k*n*x*cos(phi-atan(z/x)*(n-1.0)-atan((a2+z)/(a1+x))*(k+1.0))*pow(x*x+z*z,n*(1.0/2.0)-3.0/2.0)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0/2.0)*(n-1.0)*2.0)*(1.0/2.0))*sqrt((n*n)*pow(x*x+z*z,n-1.0)+(b0*b0)*(k*k)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),-k-1.0)+b0*k*n*cos(phi-atan(z/x)*(n-1.0)-atan((a2+z)/(a1+x))*(k+1.0))*pow(x*x+z*z,n*(1.0/2.0)-1.0/2.0)*pow(pow(a1+x,2.0)+pow(a2+z,2.0),k*(-1.0/2.0)-1.0/2.0)*2.0))/(cosh(sin(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*sin(phi-k*atan((a2+z)/(a1+x)))*pow(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z,k*(-1.0/2.0)))*sqrt(f*f+1.0)+f*cos(cos(n*atan(z/x))*pow(x*x+z*z,n*(1.0/2.0))-b0*cos(phi-k*atan((a2+z)/(a1+x)))*pow(a1*x*2.0+a2*z*2.0+a1*a1+a2*a2+x*x+z*z,k*(-1.0/2.0)))));
                }

                /// 7. Vx
                if ( dim0 > 6){   QNT[6][i1][i2] = 0.0;   }

                /// 8. Vy
                if ( dim0 > 7){   QNT[7][i1][i2] = (pNC->Va)*Vb;   }

                /// 9. Vz
                if ( dim0 > 8){   QNT[8][i1][i2] = 0.0;   }
            }
        }
    } // end of the function <calc_bckg_qnt>


    /// control output in text file
    void textfile_write( const std::string& file_name)
    {
        char buffer [256];  int counter;

        counter = sprintf (buffer, "\nDimensionless Background Parameters:\n\n");
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "a1\t:\t%+g\n"       ,pBckg->a1);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "a2\t:\t%+g\n"       ,pBckg->a2);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "b0\t:\t%+g\n"       ,pBckg->b0);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "phi\t:\t%+g rad\n"  ,pBckg->phi);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "f\t:\t%+g\n"        ,pBckg->f);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "k\t:\t%+g\n"        ,pBckg->k);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "n\t:\t%+g\n"        ,pBckg->n);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "rho_b\t:\t%+g\n"    ,pBckg->rho_b);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "p_b\t:\t%+g\n"      ,pBckg->p_b);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "gdfd\t:\t%+g\n\n"   ,pBckg->gdfd);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "Ti\t:\t%+g\n"       ,pBckg->Ti);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "Te\t:\t%+g\n"       ,pBckg->Te);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "Ti2Te\t:\t%+g\n\n"  ,pBckg->Ti2Te);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "Ve\t:\t%+f\t\telectron velocity\n"  ,pBckg->Ve);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "Vi\t:\t%+f\t\tion velocity\n"       ,pBckg->Vi);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "Vc\t:\t%+f\t\tcurrent velocity\n"   ,pBckg->Vc);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "Vb\t:\t%+f\t\tbulk velocity\n\n"    ,pBckg->Vb);
        TextOut.insert( TextOut.end(), buffer );

        int out = std::remove(file_name.c_str());     // delete old file if any

        std::ofstream file;
        file.open(file_name.c_str());
            for (size_t j=0; j<TextOut.size(); j++)  { file << TextOut[j].data(); }
        file.close();
    } // end of the function <textfile_write>
};

#endif /// CLASS_BCKG_KASKO_H
