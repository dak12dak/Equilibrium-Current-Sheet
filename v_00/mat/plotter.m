function varargout = plotter(varargin)
%%% load file with background configuration stored in grid-format
%%% output (optional): 
% array of background quantities val[dimq][dimx][dimz]; x; z.
    path1 = '';
    path0 = pwd();  cd ..,
    path1 = pwd();  cd(path0),   path1 = fullfile(path1,'run');

    if ( nargin > 0 ),  RunDir = varargin{1};
    else                RunDir = path1;
    end
     
    if ( nargin > 1 ),  BcgrFileName = varargin{2};   
    else                BcgrFileName = 'background.grid';   
    end
    
    BcgrFileName = fullfile(RunDir,BcgrFileName);

    [dim, grd, val] = f_grid_file_read(BcgrFileName);      %#ok<ASGLU>
    sz = size(val);
    if ( numel(sz) == 3 )
        x = grd(1:sz(2),1);
        z = grd(1:sz(3),2);
    else
        x = grd(1:sz(2),1);
        z = grd(1:sz(2),2);
    end
    [X,Z] = meshgrid(x,z);  %assignin('base','QNT',val);
    
    if (nargout>0),     varargout{1} = val;     end,
    if (nargout>1),     varargout{2} = x;       end,
    if (nargout>2),     varargout{3} = z;       end,
    
    % specify the number of the contour curves
    un0 = get(0,'Units');           set(0,'Units','centimeters');     
    scrsz = get(0,'ScreenSize');    set(0,'Units',un0);
    ContNum = round(scrsz(4));

    for j = 1:1:sz(1)
        
        No = j;      titstr = ['Quantity No ',num2str(No)];
        figure('name',['Self-Control: Q ',num2str(No)],...
            'units','normalized','position',[0.1, 0.1, 0.8, 0.8]); hold on;
        
        if ( dim(1) > 1 ) && ( dim(2) > 1 )
            
            qnt = squeeze(val(j,:,:));
             
            if (  abs( max(max(qnt)) - min(min(qnt)) ) < 3*eps  )
                mesh(X,Z,qnt.');              colorbar;  view(3);   grid on;
            else
                contourf(X,Z,qnt.',ContNum);  colorbar;  view(2);   grid off;
            end
            set(gca,'FontSize',16);     axis tight;          
            xlabel('x');  ylabel('z','rotation',0);  title(titstr);
            
        else
            
            qnt = val(j,:);
            
            if (  ( dim(1) == 1 ) && ( dim(2) == 1 )  )
                
                plot3(x,z,qnt,'Marker','*','MarkerSize',6); 
                set(gca,'FontSize',16);   axis tight;   grid on;  view(3);     
                xlabel('x'); ylabel('z','rotation',0); title(titstr);
                
            elseif (  ( dim(1) > 1 ) && ( dim(2) == 1)  )
                
                plot(x,qnt,'Linewidth',2); 
                set(gca,'FontSize',16);     axis tight;     grid on;     
                xlabel('x'); title(titstr);
                
            else
                
                plot(z,qnt,'Linewidth',2); 
                set(gca,'FontSize',16);     axis tight;     grid on;     
                xlabel('z'); title(titstr);
                
            end
        end
    end
    delete *.asv
end
%%%                 end of the function <plotter>
%%%------------------------------------------------------------------------

function [dim, grd, val] = f_grid_file_read(file,varargin)
%%% read text or binary file in grid format

    fid = fopen(file, 'r');

    sysdim = fscanf(fid, '%d', 1);      
        
    if ( ~isempty(sysdim) )             % text-file

        spacedim = fscanf(fid, '%d', 1);

        dim = fscanf(fid, '%d', spacedim);

        grd = zeros(max(dim), spacedim);

        for k=1:spacedim

            grd(1:dim(k), k) = fscanf(fid, '%f', dim(k));

        end

        val_raw = fscanf(fid, '%f'); 

    else                                % binary file
 
        specif = 'double';              % default
        if (~isempty(varargin)),    specif = varargin{1};       end,
        
        fseek(fid, 0, -1);          

        sysdim = fread(fid, 1, 'int');

        spacedim = fread(fid, 1, 'int');

        dim = fread(fid, spacedim, 'int');

        grd = zeros(max(dim), spacedim);

        for k=1:spacedim

            grd(1:dim(k), k) = fread(fid, dim(k), specif);

        end
        
        val_raw = fread(fid, inf, specif);

    end
    fclose(fid);  
    
    val = reshape(val_raw, [sysdim dim']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%  T H E   E N D   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%