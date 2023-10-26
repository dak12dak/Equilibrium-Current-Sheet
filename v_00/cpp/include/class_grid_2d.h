#ifndef CLASS_GRID_2D_H
#define CLASS_GRID_2D_H

#include "common.h"

class GRID_2D
{

private:
    const size_t spacedim;      /// space dimension cannot be changed!
    size_t dimx, dimz;          /// number of nodes
    double xmin, xmax, dx;      /// ---- in normalized units ----
    double zmin, zmax, dz;      /// ---- in normalized units ----

public:

    /// structure, storing grid parameters
    struct GRID_PARAM
    {
        size_t spacedim;        /// space dimension
        size_t dimx, dimz;      /// number of nodes
        double xmin, xmax, dx;  /// ---- in user-defined units ----
        double zmin, zmax, dz;  /// ---- in user-defined units ----
        double L;               /// scaling factor
    }   GridParams;

    GRID_PARAM * pGrid;         /// pointer to grid parameters structure

    std::vector<double> vX;     /// in user-defined units
    std::vector<double> vZ;     /// in user-defined units

    /// destructor
    ~GRID_2D(void) {}

    /// constructor
    GRID_2D(std::string infile, double scale) : spacedim(2)
    {
        pGrid = &GridParams;                /// pointer to grid parameters structure

        (pGrid->spacedim) = spacedim;       /// space dimension

        (pGrid->L) = scale;                 /// space scale

        read_grid_param(infile.c_str());    /// read input file


        double x[(const size_t)dimx];

        com::make_uniform_vector<double>( dimx, pGrid->xmin, pGrid->xmax, pGrid->dx, &x[0] );

        vX.reserve(dimx);       memcpy(vX.data(), (const void*)x, sizeof x);


        double z[(const size_t)dimz];

        com::make_uniform_vector<double>( dimz, pGrid->zmin, pGrid->zmax, pGrid->dz, &z[0] );

        vZ.reserve(dimz);       memcpy(vZ.data(), (const void*)z, sizeof z);
    } // end of the constructor


    /// read grid and normalization parameters
    void read_grid_param(const char* infile)
    {
        char  tline[256];        //std::cout <<"grid parameters input file: " << infile << "\n";
        char* tmp = NULL;
        int   out = 0;      int c = 0;
        double  a = 0.0; double b = 0.0;

        FILE *fid;
        if( (fid = fopen(infile, "r")) == NULL )
        {
            out = printf("Error: file %s was not opened\a\a\a\n", infile);
            exit(0);
        }

        out = fscanf(fid, "%d\n", &c);    tmp = fgets(tline, 256, fid);     //std::cout << "dim = " << c << "\n";
        if ( (size_t)c != (pGrid->spacedim) )
        {
            out = printf("Error: inconsistent space dimension: %d\a\a\a\n", c);
            fclose(fid); exit(0);
        }

        tmp = fgets(tline, 256, fid);                                       //std::cout << tline;
        out = fscanf(fid, "%lf\n", &a);   tmp = fgets(tline, 256, fid);     //std::cout << "xmin = " << a << "\n";
        out = fscanf(fid, "%lf\n", &b);   tmp = fgets(tline, 256, fid);     //std::cout << "xmax = " << b << "\n";
        out = fscanf(fid, "%d \n", &c);   tmp = fgets(tline, 256, fid);     //std::cout << "dimx = " << c << "\n";

        xmin = a;  xmax = b;  dimx = (size_t)c;
        if ( a > b )
        {
            out = printf("Error: inconsistent vector order: xmin = %lf, xmax = %lf\a\a\a\n", a, b);
            fclose(fid); exit(0);
        }

        tmp = fgets(tline, 256, fid);                                       //std::cout << tline;
        out = fscanf(fid, "%lf\n", &a);   tmp = fgets(tline, 256, fid);     //std::cout << "zmin = " << a << "\n";
        out = fscanf(fid, "%lf\n", &b);   tmp = fgets(tline, 256, fid);     //std::cout << "zmax = " << b << "\n";
        out = fscanf(fid, "%d \n", &c);   tmp = fgets(tline, 256, fid);     //std::cout << "dimz = " << c << "\n";

        zmin = a;  zmax = b;  dimz = (size_t)c;
        if ( a > b )
        {
            out = printf("Error: inconsistent vector order: zmin = %lf, zmax = %lf\a\a\a\n", a, b);
            fclose(fid); exit(0);
        }
        fclose(fid);

        if ( dimx > 1 )  { dx = (xmax-xmin)/(dimx-1); }
        else { dx = 0.0; }

        if ( dimz > 1 )  { dz = (zmax-zmin)/(dimz-1); }
        else { dz = 0.0; }

        (pGrid->dimx) = dimx;
        (pGrid->xmin) = xmin * (pGrid->L);
        (pGrid->xmax) = xmax * (pGrid->L);
        (pGrid->dx)   = dx   * (pGrid->L);

        (pGrid->dimz) = dimz;
        (pGrid->zmin) = zmin * (pGrid->L);
        (pGrid->zmax) = zmax * (pGrid->L);
        (pGrid->dz)   = dz   * (pGrid->L);

        /// self-control output
        std::cout << "grid settings file:\t"  << infile  << "\n";
        std::cout << "\nGrid Parameters:\n\n";
        std::cout << "spacedim\t"  << (pGrid->spacedim)  << "\n";
        std::cout << "scale\t\t"   << (pGrid->L)         << "\n";

        std::cout << "xmin\t\t"    << (pGrid->xmin)      << "\n";
        std::cout << "xmax\t\t"    << (pGrid->xmax)      << "\n";
        std::cout << "dx\t\t"      << (pGrid->dx)        << "\n";
        std::cout << "dimx\t\t"    << (pGrid->dimx)      << "\n";

        std::cout << "zmin\t\t"    << (pGrid->zmin)      << "\n";
        std::cout << "zmax\t\t"    << (pGrid->zmax)      << "\n";
        std::cout << "dz\t\t"      << (pGrid->dz)        << "\n";
        std::cout << "dimz\t\t"    << (pGrid->dimz)      << "\n\n";
    } // end of the function <read_grid_param>
};
#endif // CLASS_GRID_2D_H
