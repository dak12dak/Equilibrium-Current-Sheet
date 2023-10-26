%%% create input file for Yoon-Sem model configuration
%clear all

infname = 'setup_bckg_yoonsem.dat';
    
a = 0;         alpha = 0;       % alpha in degrees!
b = 10;        beta  = 0;       % beta in degrees!
f  = 0;
k  = 1;
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
str{ii} = sprintf('%s\r\n','%%% setup background configuration of Yoon-Semenov type');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% -----------------------  parameter sets  --------------------------------');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','% a = 0, alpha  = 0,  b  = 0,  beta  = 0,  f  = 0,  k  = 0 | 1 -- Harris    %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','% a = 0, alpha  = 0,  b  = 0,  beta  = 0,  f <> 0,  k  = 0 | 1 -- Fadeev    %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','% a = 0, alpha  = 0,  b <> 0,  beta  = 0,  f  = 0,  k  =     1 -- Kan       %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%        alpha  = 0,  b <> 0,  beta  = 0,  f <> 0,  k  =     1 -- Manankova %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%        alpha  = 0,           beta  = 0,           k ~= 0 | 1 -- Yoon      %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%        alpha ~= 0,           beta ~= 0,                      -- Semenov   %');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% -------------------------------------------------------------------------');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% Parameters <a> and <alpha> (IN DEGREES)');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%f\t%f\r\n',a,alpha);
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% Parameters <b> and <beta> (IN DEGREES)');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%f\t%f\r\n',b,beta);
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%%% Parameters <f> and <k>');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','%');
ii = ii + 1;
str{ii} = sprintf('%f\t%f\r\n',f,k);
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
str{ii} = sprintf('%s\r\n','%%% Parameter <gdfd> #guide field (By = constant)');
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
str{ii} = sprintf('%s\r\n','%%% NB: <a> has opposite sign as compared to Eq. 18-19 of Yoon & Lui (2005) JGR 110, A01202');

for j = 1:1:length(str)
    fprintf(fid,'%s',str{j});
end
fclose(fid);
%%%------------------------------------------------------------------------
delete *.asv
%%%------------------------------------------------------------------------