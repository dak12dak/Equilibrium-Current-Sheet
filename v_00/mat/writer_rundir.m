%clear all
outfname = 'setup_bckgr_rundir.dat';

path1 = '';
path0 = pwd();  cd ..,
path1 = pwd();  cd(path0),   path1 = fullfile(path1,'run');

rundir = fullfile(path1,'');


fid = fopen(outfname,'w');
if (fid == -1) 
    disp([mfilename,' ERROR: Cannot open file <',outfname,'>']);
    return,
end

fprintf(fid,'%s',rundir);

fclose(fid);