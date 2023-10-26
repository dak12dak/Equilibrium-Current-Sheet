%%% create input file for Kasko model configuration
%clear all

infname = 'setup_bckg_kasko.dat';

a1 = 0;         a2  = 0;
b0 = 10;        phi = 0;       % phi in degrees!
f  = 0;
k  = 1;
n  = 1;
p_b = 0;
gdfd = 0;
rho_b = 0.1;
Ti = 0.25;
Vi = 0.1;
%%%------------------------------------------------------------------------

fid = fopen(infname,'w');
if (fid == -1) 
    disp([mfilename,' ERROR: Cannot open file <',infname,'>']);
    return,
end
ii = 0;

ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% setup background configuration of Kan-Semenov-Korovinskiy type [Ann. Geophys.,36, 641-653, 2018]');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% -----------------------  parameter sets  ---------------------------------------%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%                                        f  = 0,  n = 1,  k  =     0 -- Harris      %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','% a1 = 0, a2  = 0,  b0  = 0,  phi  = 0,  f  = 0,  n = 1,  k  = 0 | 1 -- Harris      %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','% a1 = 0, a2  = 0,  b0  = 0,  phi  = 0,  f <> 0,  n = 1,  k  = 0 | 1 -- Fadeev      %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','% a1 = 0, a2  = 0,  b0 <> 0,  phi  = 0,  f  = 0,  n = 1,  k  =     1 -- Kan         %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%         a2  = 0,  b0 <> 0,  phi  = 0,  f <> 0,  n = 1,  k  =     1 -- Manankova   %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%         a2  = 0,  b0 <> 0,  phi  = 0,           n = 1,  k ~= 0 | 1 -- Yoon        %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%         a2 ~= 0,  b0 <> 0,  phi ~= 0,           n = 1,      	     -- Semenov     %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%                                                 n ~=1,             -- Korovinskiy %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% --------------------------------------------------------------------------------%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% Parameters <a1> and <a2>     # shift in x and z directions');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%f\t%f\r\n',a1,a2);
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% Parameters <b0> and <phi>    # ( phi in DEGREES! )');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%f\t%f\r\n',b0,phi);
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% Parameters <f>, <k>, and <n>');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%f\t%f\t%f\r\n',f,k,n);
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% Parameters <rho_b> and <p_b> #extra background mass density and pressure');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%f\t%f\r\n',rho_b,p_b);
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% Parameter <gdfd>   #  guide field (By = constant)');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%f\r\n',gdfd);
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% Parameter <Ti>   #  0 <= Ti <= 1/2. Units: T0 = 2(Ti+Te).');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%f\t\t%s\r\n',Ti,'# For zero Vb set Ti = Me/2/(Me+Mp) = 2.721602872444568e-04');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% Parameter <Vi>   #  in Va units ');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%f\t\t\r\n',Vi);
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% NB: <a1> and <a2> have opposite signs as compared to Eq. 16 of Korovinskiy et al.(2018) Ann. Geophys.,36, 641-653');

for j = 1:1:length(str)
    fprintf(fid,'%s',str{j});
end
fclose(fid);
%%%------------------------------------------------------------------------
delete *.asv
%%%------------------------------------------------------------------------