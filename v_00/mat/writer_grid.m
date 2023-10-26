%%% create input file with grid settings
%clear all

outfname = 'setup_bckgr_grid.dat';

spacedim = 2;

xmin = 04;		
xmax = 24;
dimx = 101;		

zmin = -3;		
zmax = +3;		
dimz = 121;		
%%%------------------------------------------------------------------------

fid = fopen(outfname,'w');
if (fid == -1) 
    disp([mfilename,' ERROR: Cannot open file <',outfname,'>']);
    return,
end
ii = 0;

ii = ii + 1;
str{ii} = sprintf('%d\t\t%s\r\n',spacedim,'# space dimension');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','#################################');
ii = ii + 1;
str{ii} = sprintf('%f\t\t%s\r\n',xmin,'# xmin');
ii = ii + 1;
str{ii} = sprintf('%f\t\t%s\r\n',xmax,'# xmax');
ii = ii + 1;
str{ii} = sprintf('%d\t\t%s\r\n',dimx,'# number of nodes');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','#################################');
ii = ii + 1;
str{ii} = sprintf('%f\t\t%s\r\n',zmin,'# zmin');
ii = ii + 1;
str{ii} = sprintf('%f\t\t%s\r\n',zmax,'# zmax');
ii = ii + 1;
str{ii} = sprintf('%d\t\t%s\r\n',dimz,'# number of nodes');
ii = ii + 1;
str{ii} = sprintf('%s\r\n','#################################');


for j = 1:1:length(str)
    fprintf(fid,'%s',str{j});
end
fclose(fid);
%%%------------------------------------------------------------------------
delete *.asv
%%%------------------------------------------------------------------------