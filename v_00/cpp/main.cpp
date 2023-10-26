#include "class_bckg_yoonsem.h"
#include "class_bckg_kasko.h"

int main(int argc, char** argv)
{
    /*************************************************************************************************************/
    /******************************   THIS  BLOCK  IS  TO  STAY  IN  MAIN()   ************************************/
    /*************************************************************************************************************/

    /// -----------------------   CREATE   MODEL-SPECIFIED   BACKGROUND   OBJECT   -------------------------------

    BACKGROUND * pBCKG = NULL;         /// declare pointer to parent background object, defined in "class_bckg.h"
    BACKGROUND::INI_BCKG_SETT * pIBS = new BACKGROUND::INI_BCKG_SETT();  /// get initial settings

    switch (pIBS->MODEL)    /// cast parent background object to a model-specified child, get background parameters
    {
        case YOONSEM: { pBCKG = new BCKG_YOONSEM( pIBS ); }  break;      /// defined in "class_bckg_yoonsem.h"
        case KASKO:   { pBCKG = new BCKG_KASKO  ( pIBS ); }  break;      /// defined in "class_bckg_kasko.h"
        default:      { std::cout << "\n\nError: unknown background model: " << pIBS->ModelName <<"\n\n"; exit(0);}
    }   if ( pIBS )   { delete pIBS;    pIBS = NULL; }    // don't need it any more
    /*************************************************************************************************************/


    /*************************************************************************************************************/
    /*******************************   PLACE  THIS  BLOCK  WHEREVER  YOU  WANT   *********************************/
    /*************************************      include "class_bckg.h"      **************************************/
    /*******************************      pass argument <BACKGROUND* pBCKG>      *********************************/
    /*************************************************************************************************************/


    /// -----------------   SPECIFY   PARAMETERS   OF   THE   UNIFORM   TWO-DIMENSIONAL   GRID   -----------------

    size_t dimx;     double* x;       /// pointer to array of x-coordinates <double x[dimx]>.
    size_t dimz;     double* z;       /// pointer to array of z-coordinates <double z[dimz]>.
    GRID_2D * pGRID = NULL;           /// defined in "class_grid_2d.h". do not remove! required for text output

    if (1)  /// Option 1: generate 2d grid using setup file
    {
        pGRID = new GRID_2D( pBCKG->InOutFnames[ENUM_MNM_COUNT], (pBCKG->pNC)->L );

        dimx = (pGRID->pGrid)->dimx;     x = (pGRID->vX).data();    /// in user-defined units
        dimz = (pGRID->pGrid)->dimz;     z = (pGRID->vZ).data();    /// in user-defined units
    }

    else    /// Option 2: get quantities { dimx, dimz, x, z }  from somewhere else.  For example:
    {
        const size_t sx = 31;        double vx[sx];        dimx = (size_t)sx;
        for ( size_t j=0; j<dimx; ++j )     vx[j] = ( +4.0 + 0.2*j ) * ((pBCKG->pNC)->L);   x = &vx[0];

        const size_t sz = 41;        double vz[sz];        dimz = (size_t)sz;
        for ( size_t j=0; j<dimz; ++j )     vz[j] = ( -2.0 + 0.1*j ) * ((pBCKG->pNC)->L);   z = &vz[0];
    }
    /*************************************************************************************************************/


    /// -----------------   CALCULATE   BACKGROUND   MAGNETOPLASMA   QUANTITIES   QNT[q](x,z)   ------------------
    ///
    /// QNT[0](x,z) = Psi,      magnetic potential
    /// QNT[1](x,z) = Rho,      mass density = n*(Mp+Me)
    /// QNT[2](x,z) = Pg ,      gas pressure
    /// QNT[3](x,z) = Bx ,      magnetic field
    /// QNT[4](x,z) = By ,      magnetic field
    /// QNT[5](x,z) = Bz ,      magnetic field
    /// QNT[6](x,z) = Vx ,      plasma bulk velocity = (Vi*Mp + Ve*Me)/(Mp+Me)
    /// QNT[7](x,z) = Vy ,      plasma bulk velocity
    /// QNT[8](x,z) = Vz ,      plasma bulk velocity
    ///
    /// Reference system: GSE, rotated for 180 degrees around z axis


    size_t dimq = pBCKG->dimq;      /// <dimq> is debugging parameter. Normally, <dimq> is to be set to 9.

    double*** QNT = com::array3_generator<double>( (size_t) dimq, (size_t) dimx, (size_t) dimz ); /// create array

    pBCKG->calc_bckg_qnt( (size_t) dimq, (size_t) dimx, (size_t) dimz, (double*) x, (double*) z, (double***) QNT );
    /*************************************************************************************************************/


    /// -----------------------------------------------   OUTPUT   -----------------------------------------------

    /// save background in grid-format
    com::array3_grid_file_write<double>( (size_t) dimq, (size_t) dimx, (size_t) dimz,
        (double*) x, (double*) z, (double***) QNT, pBCKG->InOutFnames[ENUM_MNM_COUNT+1] );


    /// write run settings and grid+background parameters in text-file
    pBCKG->textfile_output( pGRID );
    /*************************************************************************************************************/


    /// ------------------------------------------   MEMORY   CLEANING   -----------------------------------------

    if ( pBCKG )  {  delete pBCKG;      pBCKG = NULL;  }
    if ( pGRID )  {  delete pGRID;      pGRID = NULL;  }
    if (  QNT  )  {  com::array3_destroyer( (double***) QNT, (size_t) dimq, (size_t) dimx );  }
    /*************************************************************************************************************/

    return 0;
}/****************************************************************************************************************/
/***************************************   END  OF  THE  FUNCTION  MAIN   ****************************************/
/*****************************************************************************************************************/
