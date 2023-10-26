#ifndef CLASS_BCKG_H
#define CLASS_BCKG_H

#include "class_grid_2d.h"

/// list of available models:

enum  MODEL_NAME  { EMPTY = 0 , YOONSEM , KASKO , ENUM_MNM_COUNT };


class BACKGROUND
{

public:

    /// structure, storing initial settings. <class BACKGROUND> has no member instances.
    struct INI_BCKG_SETT
    {
        char ModelName[256];                            /// background model name
        char units[16];                                 /// units of measurement
        int  MODEL;                                     /// background model number
        size_t dimq;                                    /// the number of magnetoplasma quantities
        std::string InOutFnames[ENUM_MNM_COUNT+4];      /// input/output filenames

        /// destructor
        ~INI_BCKG_SETT(void) {}

        /// constructor
        INI_BCKG_SETT(void)
        {
            assign_file_names();
            get_initial_settings();
        }


        /// assign input/output file names
        void assign_file_names(void)
        {
            char rundirsetup[256] = "setup_bckgr_rundir.dat";      /// file with run directory name
            int  sz = 0;

            sz++;   InOutFnames[EMPTY] = "setup_bckg.dat" ;                    /// initial settings

            sz++;   InOutFnames[YOONSEM] = "setup_bckg_yoonsem.dat" ;          /// yoonsem model settings
            sz++;   InOutFnames[KASKO]   = "setup_bckg_kasko.dat"   ;          /// kasko model settings

            sz++;   InOutFnames[ENUM_MNM_COUNT] = "setup_bckgr_grid.dat" ;     /// grid setup filename

            sz++;   InOutFnames[ENUM_MNM_COUNT+1] = "background.grid" ;        /// .grid-save-file-name
            sz++;   InOutFnames[ENUM_MNM_COUNT+2] = "background.info" ;        /// text output filename

            sz++;   InOutFnames[ENUM_MNM_COUNT+3] = "" ;           /// run directory name is so far unknown

            /// read run directory name
            std::string rundir = get_rundir( (const char*) rundirsetup );      /// get run directory name from input file

            if ( rundir.empty() )  { rundir = com::get_current_dir();   }      /// if empty, get current directory name

            char slash = com::slash();                                         /// get slash symbol

            if ( *rundir.rbegin() != slash ) { rundir.push_back(slash); }      /// force "rundir" to terminate with slash

            for (int j=0; j<sz; ++j)
            {  InOutFnames[j].insert( 0, rundir.data() );               }      /// insert "rundir" before filenames

            std::cout << "\n\nrun directory:\t" << InOutFnames[sz-1] << "\n";
        } // end of the function <assign_file_names>


        /// read run directory name
        std::string get_rundir(const char* infile)
        {
            std::string rundir;
            char tline[FILENAME_UMAX];           //std::cout<<"FILENAME_UMAX = "<<FILENAME_UMAX<<"\n\n";

            FILE *fid;
            if( (fid = fopen(infile, "r")) == NULL )
            {
                std::cout << "\n\nError: file <" << infile << "> was not opened. Use current directory.\n";
                return rundir;
            }

            fseek (fid, 0, SEEK_END);
            int fsize = (int) ftell(fid);

            if ( fsize == 0 )
            {
                std::cout << "\n\nWarning: file <" << infile << "> is empty. Use current directory.\n";
                fclose(fid);
                return rundir;
            }

            rewind(fid);
            char* tmp = fgets(tline, FILENAME_UMAX, fid);
            fclose(fid);

            tline[strcspn(tline, "\r\n")] = '\0';       /// remove new line character;
            rundir.assign((const char*)tline);
            return rundir;
        } // end of the function <get_rundir>


        /// get initial settings from the text-file
        void get_initial_settings(void)
        {
            char infile[FILENAME_UMAX];
            char* tmp = strcpy( infile, InOutFnames[EMPTY].c_str() );
            std::cout  << "\n\ninitial settings file\t\t:\t" << infile <<"\n";

            char tline[256];
            int out = 0;    int val = 0;
            size_t len = 0;

            FILE *fid;
            if( (fid = fopen(infile, "r")) == NULL )
            {
                out = printf("\nError: file %s was not opened\a\a\a\n", infile);
                exit(0);
            }

            for (int j=0; j<2; ++j) { tmp = fgets(tline, 256, fid); }     /// skip comments
            tmp = fgets(tline, 256, fid);                              /// background model name
            tline[strcspn(tline, "\r\n")] = '\0';                   /// remove new line character;

            /// cast model name to upper case:
            try {  len = com::strtoupper_autovec( (char*)ModelName, sizeof(ModelName), (const char*)tline );  }
            catch (const std::exception& e)
            {
                std::cout << "\na standard exception was caught, with message '"<< e.what() << "'\n";
                tmp = strcpy((char*)ModelName,(const char*)tline);
            }

            for (int j=0; j<3; ++j) { tmp = fgets(tline, 256, fid); }     /// skip comments
            out = fscanf(fid, "%d\n", &val);

            if ( val < 0 )
            {
                out = printf("\nError! Illegal value '%d' of the parameter 'dimq'!\a\a\a\n\n",val);
                fclose(fid);
                exit(0);
            }
            else  {  dimq = (size_t)val;    }                          /// number of background quantities

            for (int j=0; j<3; ++j) { tmp = fgets(tline, 256, fid); }     /// skip comments
            tmp = fgets(tline, 256, fid);                              /// units
            tline[strcspn(tline, "\r\n")] = '\0';                   /// remove new line character;
            fclose(fid);

            /// cast units to upper case:
            try {  len = com::strtoupper_autovec( (char*)units, sizeof(units), (const char*)tline );  }
            catch (const std::exception& e)
            {
                std::cout << "\na standard exception was caught, with message '"<< e.what() << "'\n";
                tmp = strcpy((char*)units,(const char*)tline);
            }

            MODEL = get_model_number((const char*)ModelName);       /// get model number from enum

            std::cout << "Background model name\t\t:\t"         << ModelName << "\n"  ;
            std::cout << "Background model number\t\t:\t"       << MODEL     << "\n"  ;
            std::cout << "Number of background quantities\t:\t" << dimq      << "\n\n";
            std::cout << "Units\t:\t"  << units << "\n\n";
        } // end of the function <get_initial_settings>


        /// get model number with using enum list
        int get_model_number(const char* tline)
        {
            int model = EMPTY;
            char buf1[256],  buf2[256];
            char* tmp = NULL;

              tmp = strcpy(buf1, "yoonsem");  tmp = strcpy(buf2, "YOONSEM");
            if (  ( strcmp(buf1, tline) == 0 ) || ( strcmp(buf2, tline) == 0 )  )
            {
                model = YOONSEM;
                return model;
            }

              tmp = strcpy(buf1, "kasko");    tmp = strcpy(buf2, "KASKO");
            if (  ( strcmp(buf1, tline) == 0 ) || ( strcmp(buf2, tline) == 0 )  )
            {
                model = KASKO;
                return model;
            }

            return model;
        } // end of the function <get_model_number>
    }; // end of the structure <INI_BCKG_SETT>


    /// structure, storing normalization coefficients. <class BACKGROUND> has a member instance.
    struct NORM_COEFF
    {
        char units[16];         /// units = { norm, si, cgs }
        double L;               /// scaling factor for: length
        double B0;              /// scaling factor for: magnetic field
        double N0;              /// scaling factor for: number density
        double P0;              /// scaling factor for: pressure
        double T0;              /// scaling factor for: temperature
        double Va;              /// scaling factor for: velocity
        double J0;              /// scaling factor for: current
        double Rho0;            /// scaling factor for: mass density
        double Psi0;            /// scaling factor for: magnetic potential
        double di;              /// ion inertial length = c/omega_p

        /// destructor
        ~NORM_COEFF(void) {}

        /// constructor
        NORM_COEFF(void)
        {
            set_default_norm();
        }

        /// set default values of normalization coefficients
        void set_default_norm(void)
        {
            char* tmp = strcpy(units,"norm");   //std::cout << "units\t:\t" << units << "\n";
            L  = 1.0;                           //std::cout << "L\t:\t"     << L     << "\n";
            B0 = 1.0;                           //std::cout << "B0\t:\t"    << B0    << "\n";
            N0 = 1.0;                           //std::cout << "N0\t:\t"    << N0    << "\n";
            T0 = 1.0;                           //std::cout << "T0\t:\t"    << T0    << "\n";
            Va = 1.0;                           //std::cout << "Va\t:\t"    << Va    << "\n";
            P0 = 1.0;                           //std::cout << "P0\t:\t"    << P0    << "\n";
            J0 = 1.0;                           //std::cout << "J0\t:\t"    << J0    << "\n";
            Rho0 = 1.0;                         //std::cout << "Rho0\t:\t"  << Rho0  << "\n";
            Psi0 = 1.0;                         //std::cout << "Psi0\t:\t"  << Psi0  << "\n";
            di = com::NaN<double>();            //std::cout << "di\t:\t"    << di    << "\n";
        } // end of the function <set_default_norm>

    }   NormCoeff;

    NORM_COEFF * pNC;           /// pointer to normalization coefficients structure

    char ModelName[256];        /// background model name

    int MODEL;                  /// background model number

    size_t dimq;                /// the number of magnetoplasma quantities to be calculated

    std::vector<std::string> InOutFnames;   /// vector of strings for input/output filenames

    std::vector<std::string> TextOut;       /// vector of strings for text-file output


    /// destructor
    virtual ~BACKGROUND(void){}


    /// constructor
    BACKGROUND(void)
    {
        INI_BCKG_SETT  IniSett = INI_BCKG_SETT();   /// create structure with initial settings
        INI_BCKG_SETT* pIniSett = &IniSett;         /// create pointer to this structure
        assign_bckg_members(pIniSett);              /// assign values to all class members
    }


    /// constructor
    BACKGROUND(INI_BCKG_SETT* pIniSett)
    {
        assign_bckg_members(pIniSett);              /// assign values to all class members
    }

    /// assign background object members
    void assign_bckg_members(INI_BCKG_SETT* pIniSett)
    {
        char* tmp = strcpy(ModelName,(const char*) pIniSett->ModelName);    /// copy model name from initial settings
        MODEL = pIniSett->MODEL;                                            /// copy model number from initial settings
        dimq  = pIniSett->dimq;                                             /// copy state vector rank from initial settings

        /// copy input/output filenames from initial settings
        size_t sz = sizeof(pIniSett->InOutFnames)/sizeof(pIniSett->InOutFnames[0]); //std::cout<<"sz = "<<sz<<"\n\n";
        for (size_t j=0; j<sz; ++j)
        {   InOutFnames.insert( InOutFnames.end(), pIniSett->InOutFnames[j] );  }

        pNC = &NormCoeff;                           /// pointer to the structure with default normalization coefficients
        tmp = strcpy(pNC->units, (const char*) pIniSett->units);   /// copy units from initial settings
        get_norm_coeff();                           /// read normalization coefficients from the text-file
    }


    /// get normalization coefficients from the text-file
    void get_norm_coeff(void)
    {
        if (  ( strcmp(pNC->units, "norm") != 0 ) && ( strcmp(pNC->units, "NORM") != 0 )  )
        {
            char infile[FILENAME_UMAX];
            char* tmp = strcpy( infile, InOutFnames[EMPTY].c_str() );
            //std::cout  << "\n\ninitial settings file:\t\t\t" << infile <<"\n";

            char tline[256];
            double val = 0.0;
            int out = 0;

            FILE *fid;
            if( (fid = fopen(infile, "r")) == NULL )
            {
                out = printf("\nError: file %s was not opened\a\a\a\n", infile);
                exit(0);
            }

            for (int j=0; j<14; ++j) { tmp = fgets(tline, 256, fid); }     /// skip header

            out = fscanf(fid, "%lf\n", &val);     tmp = fgets(tline, 256, fid);
            pNC->L = val;

            out = fscanf(fid, "%lf\n", &val);     tmp = fgets(tline, 256, fid);
            pNC->B0 = val;

            out = fscanf(fid, "%lf\n", &val);     tmp = fgets(tline, 256, fid);
            pNC->N0 = val;

            fclose(fid);
            cast_to_units();     /// cast coefficients to specified units of measurement
        }

        std::cout << "Normalization Factors:\n\n";
        std::cout << "L\t\t"      << (pNC->L  )   << "\n";
        std::cout << "di\t\t"     << (pNC->di )   << "\t\tion inertial length\n";
        std::cout << "B0\t\t"     << (pNC->B0 )   << "\n";
        std::cout << "N0\t\t"     << (pNC->N0 )   << "\n";
        std::cout << "T0\t\t"     << (pNC->T0 )   << "\n";
        std::cout << "Va\t\t"     << (pNC->Va )   << "\n";
        std::cout << "P0\t\t"     << (pNC->P0 )   << "\n";
        std::cout << "J0\t\t"     << (pNC->J0 )   << "\n";
        std::cout << "Rho0\t\t"   << (pNC->Rho0)  << "\n";
        std::cout << "Psi0\t\t"   << (pNC->Psi0 ) << "\n\n";
    } // end of the function <get_norm_coeff>


    /// cast scaling coefficient to the specified metric units
    void cast_to_units(void)
    {
        char units[16], buf1[16], buf2[16];
        char* tmp = strcpy(units,(const char*) pNC->units);          //std::cout << "units = " << units << "\n\n";

          tmp = strcpy(buf1, "norm");     tmp = strcpy(buf2, "NORM");
        if (  ( strcmp(buf1, units) == 0 ) || ( strcmp(buf2, units) == 0 )  )
        {   return;   }

        /// cast used-defined constants to SI, evaluate other constants

        double Pi  = com::Pi <double>();        // pi constant
        double ElC = com::ElC<double>();        // SI, elementary charge
        double Kb  = com::Kb <double>();        // SI, Boltzmann constant
        double Mp  = com::Mp <double>();        // SI, proton mass
        double Me  = com::Me <double>();        // SI, electron mass
        double MUo = com::MUo<double>();        // SI, vacuum permeability
        double c   = com::LS <double>();        // SI, speed of light

        pNC->L     = (pNC->L )*1.0e+6;                           /// m
        pNC->B0    = (pNC->B0)*1.0e-9;                           /// T
        pNC->N0    = (pNC->N0)*1.0e+6;                           /// m^(-3)

        pNC->P0    = (pNC->B0)*(pNC->B0)/MUo;                    /// J*m^(-3)
        pNC->T0    = (pNC->P0)/(pNC->N0)/Kb;                     /// K

        pNC->Rho0  = (pNC->N0)*(Mp+Me);                          /// kg*m^(-3)
        pNC->Va    = (pNC->B0)/sqrt( MUo*(pNC->Rho0) );          /// m*s^(-1)
        pNC->J0    = (pNC->B0)/(pNC->L)/MUo;                     /// A*m^(-2)
        pNC->Psi0  = -(pNC->B0)*(pNC->L);                        /// T*m

        pNC->di    = 1.0/ElC*sqrt(Mp/MUo/(pNC->N0));             /// m

          tmp = strcpy(buf1, "si");       tmp = strcpy(buf2, "SI");
        if (  ( strcmp(buf1, units) == 0 ) || ( strcmp(buf2, units) == 0 )  )
        {   return;   }

        /// cast constants from SI to CGS

          tmp = strcpy(buf1, "cgs");      tmp = strcpy(buf2, "CGS");
        if (  ( strcmp(buf1, units) == 0 ) || ( strcmp(buf2, units) == 0 )  )
        {
            pNC->L    = (pNC->L )*1.0e+2;                        /// cm
            pNC->B0   = (pNC->B0)*1.0e+4;                        /// G
            pNC->N0   = (pNC->N0)*1.0e-6;                        /// cm^(-3)
            /// Temperature units in CGS are not defined. Use Kelvin
            pNC->P0   = (pNC->P0)*1.0e+1;                        /// Erg*cm^(-3)
            pNC->Rho0 = (pNC->Rho0)*1.0e-3;                      /// g*cm^(-3)
            pNC->Va   = (pNC->Va)*1.0e+2;                        /// cm*s^(-1)
            //pNC->J0   = c*1.0e+2/4.0/Pi*(pNC->B0)/(pNC->L);      /// Fr*s^(-1)*cm^(-2)
            pNC->J0   = (pNC->J0)*c*1.0e-3;                      /// Fr*s^(-1)*cm^(-2)
            pNC->Psi0 = -(pNC->B0)*(pNC->L);                     /// G*cm
            pNC->di   = (pNC->di)*1.0e+2;                        /// cm
            return;
        }

        else
        {
            int out = printf("\nError: unknown system of units: %s\a\a\a\n",units);
            exit(0);
        }
    } // end of the function <cast_to_units>


    /// read background parameters
    virtual void read_bckg_config( const char* ){}


    /// evaluate temperature ratio and velocities in Harris-like models
    double eval_velocity_harris_family(double Ti, double Te, double Vi, double* Vel)
    {
        char units[16], buf1[16], buf2[16];
        char* tmp = strcpy(units,(const char*) pNC->units);          //std::cout << "units = " << units << "\n\n";
        int out = 0;

        double Eps = com::Eps<double>();
        double Inf = com::Inf<double>();
        double Mp  = com::Mp <double>();
        double Me  = com::Me <double>();
        double Ve, Vc, Vb, Ti2Te;

        if ( Te > 10.0*Eps )    Ti2Te = Ti/Te;
        else                    Ti2Te = Inf;

          tmp = strcpy(buf1, "norm");     tmp = strcpy(buf2, "NORM");
        if (  ( strcmp(buf1, units) != 0 ) && ( strcmp(buf2, units) != 0 )  )
        {
            double dinorm = (pNC->di)/(pNC->L);
            Ve = +2.0*dinorm*Te;                /// Te is normalized
            Vi = -2.0*dinorm*Ti;                /// Ti is normalized
            Vc = -2.0*dinorm*(Ti+Te);
            Vb = -2.0*dinorm*(Ti*Mp - Te*Me)/(Mp+Me);
        }

        else  /// in normalized units Vi is arbitrary, Ve = Ve(Te,Ti,Vi)
        {
            if ( Ti > 10.0*Eps )
            {
                Vi = -com::abs<double>(Vi);
                Ve = -Te/Ti*Vi;
            }

            else
            {
                out = printf("\nWarning: Ion temperature is zero, hence Vi = 0. Electron velocity is undefined.\n");
                Vi = 0.0;

                out = printf("\nEnter electron velocity:\t");
                try  {  std::cin >> Ve;  }
                catch (const std::exception& e)
                {
                    std::cout << "\na standard exception was caught, with message '"<< e.what() << "'\n";
                    out = printf("\nSet electron velocity to 1.\n");
                    Ve = 1.0;
                }
                Ve = +com::abs<double>(Ve);
            }

            Vc = Vi - Ve;
            Vb = ( Mp*Vi + Me*Ve ) / (Mp+Me);
        }

        Vb = com::round_decdig<double>(Vb);

        Vel[0] = Ve;  Vel[1] = Vi;  Vel[2] = Vc;  Vel[3] = Vb;

        return Ti2Te;
    } // end of the function <eval_velocity_harris_family>


    /// calculate background plasma quantities
    virtual void calc_bckg_qnt(size_t, size_t, size_t, double*, double*, double***){}


    /// control output in text file
    void textfile_output( GRID_2D* pGRID = NULL )
    {
        char buffer [256];  int counter;
        TextOut.insert( TextOut.end(), "PARAMETERS OF THE BACKGROUND CONFIGURATION:\n\n" );

        counter = sprintf (buffer, "Background model name:\t\t\t%s\n", ModelName);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "Background model number:\t\t%d\n", MODEL);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "Number of background quantities:\t%d\n\n", (int)dimq);
        TextOut.insert( TextOut.end(), buffer );


        counter = sprintf (buffer, "Units\t:\t%s\n\n" ,pNC->units);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "Normalization Factors:\n\n");
        TextOut.insert( TextOut.end(), buffer );

        if (std::isnan(pNC->di))
        {
            counter = sprintf (buffer, "L\t:\t%g\n" ,pNC->L);
            TextOut.insert( TextOut.end(), buffer );

            counter = sprintf (buffer, "di\t:\t%s\t\t\tion inertial length\n","NaN"  );
            TextOut.insert( TextOut.end(), buffer );
        }

        else
        {
            counter = sprintf (buffer, "L\t:\t%e\n" ,pNC->L);
            TextOut.insert( TextOut.end(), buffer );

            counter = sprintf (buffer, "di\t:\t%e\t\tion inertial length\n",pNC->di);
            TextOut.insert( TextOut.end(), buffer );
        }

        counter = sprintf (buffer, "B0\t:\t%g\n"    ,pNC->B0);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "N0\t:\t%g\n"    ,pNC->N0);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "T0\t:\t%g\n"    ,pNC->T0);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "Va\t:\t%g\n"    ,pNC->Va);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "P0\t:\t%g\n"    ,pNC->P0);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "J0\t:\t%g\n"    ,pNC->J0);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "Rho0\t:\t%g\n"  ,pNC->Rho0);
        TextOut.insert( TextOut.end(), buffer );

        counter = sprintf (buffer, "Psi0\t:\t%g\n\n",pNC->Psi0);
        TextOut.insert( TextOut.end(), buffer );


        if ( pGRID )
        {
            counter = sprintf (buffer, "\nGrid Parameters:\n\n");
            TextOut.insert( TextOut.end(), buffer );

            counter = sprintf (buffer, "spacedim:\t%d\n" ,(int)((pGRID->pGrid)->spacedim));
            TextOut.insert( TextOut.end(), buffer );

            counter = sprintf (buffer, "scale\t:\t%g\n"  ,(pGRID->pGrid)->L);
            TextOut.insert( TextOut.end(), buffer );

            counter = sprintf (buffer, "xmin\t:\t%g\n"   ,(pGRID->pGrid)->xmin);
            TextOut.insert( TextOut.end(), buffer );

            counter = sprintf (buffer, "xmax\t:\t%g\n"   ,(pGRID->pGrid)->xmax);
            TextOut.insert( TextOut.end(), buffer );

            counter = sprintf (buffer, "dx\t:\t%g\n"     ,(pGRID->pGrid)->dx);
            TextOut.insert( TextOut.end(), buffer );

            counter = sprintf (buffer, "dimx\t:\t%d\n"   ,(int)((pGRID->pGrid)->dimx));
            TextOut.insert( TextOut.end(), buffer );

            counter = sprintf (buffer, "zmin\t:\t%+g\n"  ,(pGRID->pGrid)->zmin);
            TextOut.insert( TextOut.end(), buffer );

            counter = sprintf (buffer, "zmax\t:\t%+g\n"  ,(pGRID->pGrid)->zmax);
            TextOut.insert( TextOut.end(), buffer );

            counter = sprintf (buffer, "dz\t:\t%g\n"     ,(pGRID->pGrid)->dz);
            TextOut.insert( TextOut.end(), buffer );

            counter = sprintf (buffer, "dimz\t:\t%d\n\n" ,(int)((pGRID->pGrid)->dimz));
            TextOut.insert( TextOut.end(), buffer );
        }

        std::string file_name = InOutFnames[ENUM_MNM_COUNT+2];              /// target filename
        textfile_write( (const std::string&)file_name );                    /// call file-writing function
        std::cout << "parameters are written to:\t" << file_name << "\n";   /// screen output
    } // end of the function <textfile_output>


    virtual void textfile_write( const std::string& ){}
};
#endif /// CLASS_BCKG_H
