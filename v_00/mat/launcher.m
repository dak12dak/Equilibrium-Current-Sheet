clear all,      %#ok<*NASGU> 

%MODEL = 'Erkaev';     % Harris sheet + linear Bz(x) - not implemented
MODEL = 'YoonSem';    % Kan-Fadeev-Manankova-Yoon-Semenov
MODEL = 'Kasko';      % Kan-Semenov-Korovinskiy

%%% setup uniform grid
xmin = +04;     xmax = +24;

zmin = -03;     zmax = +03;

grid_stp_pow = 3;     % grid step = (1/2)^grid_stp_pow

FigPlotMode = 1;      % plot background configuration
SelfControl = 1;      % plot numerical residual forces and saved quantities
SaveBckFile = 0;      % save background in grid-format file

%OutDataFormat = 'ascii';    
OutDataFormat = 'binary';

%%% PLASMA QUANTITIES TO BE CALCULATED SYMBOLICALLY:
%%%
QNT{1,1} = {'Psi','Rho','Pg','Bx','By','Bz','Vx','Vy','Vz','Jy','Fx','Fz'};
%%%
%%% REMOVE ANY NUMBER OF ITEMS FROM THE END, BUT DO NOT CHANGE THE ORDER!
%%% example: 
QNT{1,1} = {'Psi','Rho','Pg','Bx','By','Bz','Vx','Vy','Vz'};
%%%------------------------------------------------------------------------

switch upper(MODEL)
            
    case 'ERKAEV'
        fprintf('\nMESSAGE: model <ERKAEV> is not coded yet.\n');
        generator = [];
        return,
        
    case 'YOONSEM'
        
        generator = @backgr_generator_yoonsem; 
        
    case 'KASKO'
        
        generator = @backgr_generator_kasko;
        
    otherwise
        fprintf('\nERROR: Unknown model <%s>\n', MODEL);
        generator = [];
        return,
end
%%%------------------------------------------------------------------------

grid_stp = 1/2^grid_stp_pow;
x = (xmin:grid_stp:xmax).';
z = (zmin:grid_stp:zmax).';
[X,Z] = meshgrid(x,z);              % size [dimz, dimx]
%%%------------------------------------------------------------------------

[OutFile, BckParam, QNT{1,2}] = ...
    generator( length(QNT{1,1}), X, Z, SaveBckFile, OutDataFormat );
%%%------------------------------------------------------------------------

% specify the number of the contour curves
un0 = get(0,'Units');           set(0,'Units','centimeters');     
scrsz = get(0,'ScreenSize');    set(0,'Units',un0);
ContNum = round(scrsz(4));

%%% visualization
if      (  ( size(X,1) > 1) && ( size(X,2) > 1)  )
    DIM = 2;
    plotf = @(x,z,val) contourf(x,z,val,ContNum);XLabel = 'x';YLabel = 'z';
    
elseif  (  ( size(X,1) > 1) && ( size(X,2) == 1)  )
    DIM = 1;
    plotf = @(x,z,val) plot(z,val);         XLabel = 'z';   YLabel = ''; 
    
elseif  (  ( size(X,2) > 1) && ( size(X,1) == 1)  )
    DIM = 1;
    plotf = @(x,z,val) plot(x,val);         XLabel = 'x';   YLabel = '';
    
else
    DIM = 0;
    plotf = @(x,z,val) plot3(x,z,val,'Marker','*','MarkerSize',6);
                                            XLabel = 'x';   YLabel = 'z';
end

if FigPlotMode
    for j = 1:1:size(QNT{1,2},1)    
        figure('name',QNT{1,1}{j},'units','normalized','position',[0.1, 0.1, 0.8, 0.8]);
        qnt = squeeze(QNT{1,2}(j,:,:));
        if (  ( max(max(abs(qnt))) - min(min(abs(qnt))) < 3*eps ) && ...
              ( ~isempty(regexp(func2str(plotf),'contour','once')) )  )
            mesh(X,Z,qnt);      grid on;
        else
            plotf(X,Z,qnt);     grid off;
        end
        set(gca,'FontSize',16,'FontWeight','Normal');
        xlabel(XLabel); ylabel(YLabel,'Rotation',0); title(QNT{1,1}{j});
        axis tight; hold on; 
        if ( DIM == 2 ),    colorbar('FontSize',16);     end,
    end
end

%%% calculate numerically and plot residual forces
if  SelfControl && ( DIM > 1 ) && ( size(QNT{1,2},1) > 5 )
    
    dx = x(2)-x(1);     dz = z(2)-z(1);
    
    fcefind = @(substr) find(~cellfun(@isempty,regexp(QNT{1,1},substr)));
    
    ipg = fcefind('Pg');
    [DxPg,DzPg] = gradient(squeeze(QNT{1,2}(ipg,:,:)),dx,dz);
    
    ibx = fcefind('Bx');
    [DxBx,DzBx] = gradient(squeeze(QNT{1,2}(ibx,:,:)),dx,dz);
    
    ibz = fcefind('Bz');
    [DxBz,DzBz] = gradient(squeeze(QNT{1,2}(ibz,:,:)),dx,dz);
    
    numJy = +DzBx - DxBz;   
    numFx = -DxPg + numJy.*squeeze(QNT{1,2}(ibz,:,:));
    numFz = -DzPg - numJy.*squeeze(QNT{1,2}(ibx,:,:));
    numDivB = DxBx + DzBz;
    
    % cut edges
    ix = 1:length(x);       iz = 1:length(z);       titmess = '';
    ix = 2:length(x)-1;     iz = 2:length(z)-1;     titmess = '. Edges are cut.';

    figure('name','num Fx','units','normalized','position',[0.1, 0.1, 0.8, 0.8]);
    mesh(X(iz,ix),Z(iz,ix),numFx(iz,ix)); colorbar('FontSize',16); view(3);
    set(gca,'FontSize',16,'FontWeight','Normal');
    xlabel('x'); ylabel('z','Rotation',0); 
    title(['Numerical F_x. dx = ',num2str(dx),', dz = ',num2str(dz),titmess]);
    axis tight; hold on; grid on;

    figure('name','num Fz','units','normalized','position',[0.1, 0.1, 0.8, 0.8]);
    mesh(X(iz,ix),Z(iz,ix),numFz(iz,ix)); colorbar('FontSize',16); view(3);
    set(gca,'FontSize',16,'FontWeight','Normal');
    xlabel('x'); ylabel('z','Rotation',0); 
    title(['Numerical F_z. dx = ',num2str(dx),', dz = ',num2str(dz),titmess]);
    axis tight; hold on; grid on;
    
    figure('name','num div(B)','units','normalized','position',[0.1, 0.1, 0.8, 0.8]);
    mesh(X(iz,ix),Z(iz,ix),numDivB(iz,ix)); colorbar('FontSize',16); view(3);
    set(gca,'FontSize',16,'FontWeight','Normal');
    xlabel('x'); ylabel('z','Rotation',0); 
    title(['Numerical div(B). dx = ',num2str(dx),', dz = ',num2str(dz),titmess]);
    axis tight; hold on; grid on;
end

%%% control of the saved background
if (  SelfControl && ( ~isempty(OutFile) )  )
    plotter('',OutFile);
end
%%%-------------------------------------------------------------------------
delete *.asv
%%%-------------------------------------------------------------------------