%%% create input file containing the model type
clear all

writefname = 'setup_bckg.dat';

model = 'yoonsem';
model = 'kasko';

dimq = 9;       % 1 <= dimq <= 9

%%% scaling factors
L  = 3.1855e-00;         % length, 10^3 km 
B0 = 6.951295;           % magnetic field, nT
N0 = 0.01;              % number density, cm^(-3)
%%%------------------------------------------------------------------------

fid = fopen(writefname,'w');
if (fid == -1) 
    disp([mfilename,' ERROR: Cannot open file <',writefname,'>']);
    return,
end
ii = 0;

ii = ii + 1;
str{ii} = sprintf('%s\r\n','### model name: {yoonsem, kasko}');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','#');
ii = ii + 1;
str{ii} = sprintf('%s\r\n',model);
ii = ii + 1;
str{ii} = sprintf('%s\r\n','#');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','### number of output quantities QNT = {"Psi","Rho","Pg", "Bx","By","Bz", "Vx","Vy","Vz"}');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','#');
ii = ii + 1;
str{ii} = sprintf('%d\r\n',dimq);
ii = ii + 1;
str{ii} = sprintf('%s\r\n','#');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','### units: { norm, si, cgs }');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','#');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','norm');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','#');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','### scaling factors');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','#');
ii = ii + 1;
str{ii} = sprintf('%.6e\t\t%s\r\n',L,'# L    , length, 10^3 km');
ii = ii + 1;
str{ii} = sprintf('%.6e\t\t%s\r\n',B0,'# B0   , magnetic field, nT');
ii = ii + 1;
str{ii} = sprintf('%.6e\t\t%s\r\n',N0,'# N0   , number density, cm^(-3)');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','########################################################################################');

for j = 1:1:length(str)
    fprintf(fid,'%s',str{j});
end
fclose(fid);
%%%------------------------------------------------------------------------
delete *.asv
%%%------------------------------------------------------------------------