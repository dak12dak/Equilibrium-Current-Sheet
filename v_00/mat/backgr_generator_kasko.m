function [OutFile, BckParam, VAL] = ...
    backgr_generator_kasko(dimq,X,Z,SaveBckFile,OutDataFormat)
%%% generator of background configuration to be translated in C++
%%% INPUT:
%%% 1. dimq - number of output quantities
%%% 2,3. (X,Z)- coordinates, numeric, 2D arrays (1D, 0D are allowed also). 
%%%    Reference system is GSM, rotated for 180 degrees around Z-axis
%%%    Requirements:
%%%    a). size(X) = size(Z) = [dimz, dimx]
%%%    b). all points must belong to 1 or/and 4 quadrants (except origin).
%%% 4. SaveBckFile - {1,0} = {save, dont save} background configuration 
%%% 5. OutDataFormat - {'ascii','binary'} = save in {text, binary} file
%%% OUTPUT:
%%% 1. OutFile - the name of file, where configuration is saved.
%%%    If input parameter <SaveBckFile> = 0, OutFile = ''.
%%% 2. BckParam - structure, containing background parameters. 
%%%    Phase angles <alpha, beta> are converted to radians.
%%% 3. VAL - background configuration, numeric, size(VAL)=[dimq,dimz,dimx].

    InFile  = 'setup_bckg_kasko.dat';  % file with configuration parameters
    OutFile = 'background.grid';       % file to save configuration

    BckParam = f_read_config_Kasko(InFile);   % read input file

    x = X(1,:);     dimx = length(x);       % NB: size(X) = [dimz,dimx]
    z = Z(:,1);     dimz = length(z);       % NB: size(Z) = [dimz,dimx]

    VAL = NaN*zeros(dimq,dimz,dimx);        % Allocate output array

%     % calculate plasma quantities point-by-point
%     for ind1 = 1:1:dimz
%         for ind2 = 1:1:dimx
%             x0 = X(ind1,ind2);
%             z0 = Z(ind1,ind2);
%             VAL(:,ind1,ind2) = f_bckgr_funcs_Kasko(dimq,x0,z0,BckParam);    
%         end
%     end

    % calculate plasma quantities at 2D grid
    VAL = f_bckgr_funcs_Kasko(dimq,X,Z,BckParam); 

    % save background in grid-file if grid is 2-dimensional
    if SaveBckFile %&& ( dimx > 1 ) && ( dimz > 1 )
        tVal = zeros(size(VAL,1),size(VAL,3),size(VAL,2));  % tVal = VAL.' 
        sz = size(VAL);
        if ( numel(sz) == 3)
            for j = 1:1:size(VAL,1)
                tVal(j,:,:) = (squeeze(VAL(j,:,:))).';
            end
        else 
            for j = 1:1:size(VAL,1)
                tVal(j,:) = (VAL(j,:)).';
            end            
        end
        f_grid_file_write(x, z, [], tVal, OutFile, OutDataFormat);
        fprintf('\nBackground configuration is saved in <%s>\n',OutFile);
    else
        OutFile = '';
    end       
end
%%%------------------------------------------------------------------------
%%%   E N D   O F   T H E   F U N C T I O N   <backgr_kasko_generator>
%%%------------------------------------------------------------------------

%%% read background parameters (Kan-Semenov-Korovinskiy) from text-file
function BckParam = f_read_config_Kasko(infile)
%%% INPUT: file-name. OUTPUT: structure, containing background parameters.
    %#ok<*NASGU>
    
    % open file with configuration parameters
    [fid, status] = fopen(infile,'r');
    if ( fid < 0 )
        fprintf('\nERROR: Unable open file <%s>. Status: %s',infile,status);
        return,
    end
    
    % skip header
    for j = 1:1:14
        tline = fgetl(fid);       %fprintf('\n%s',tline);
    end
    
    % read <a1>, <a2>
    tval = fscanf(fid,'%f',2);
    a1 = tval(1); a2 = tval(2); %fprintf('\na1 = %g\ta2 = %g',a1,a2);
    tline = fgets(fid);
    
    % skip comments
    for j = 1:1:3
        tline = fgetl(fid);       %fprintf('\n%s',tline);
    end
    
    % read <b0>, <phi>
    tval = fscanf(fid,'%f',2);
    b0 = tval(1); phi = tval(2);  %fprintf('\nb0 = %g\tphi = %g',b0,phi);
    tline = fgets(fid);
    
    % skip comments
    for j = 1:1:3
        tline = fgetl(fid);       %fprintf('\n%s',tline);
    end
    
    % read <f>, <k>, <n>
    tval = fscanf(fid,'%f',3);
    f = tval(1); k = tval(2); n = tval(3);   %fprintf('\nf = %g\tk = %g\tn = %g',f,k,n);
    tline = fgets(fid);
    
    % skip comments
    for j = 1:1:3
        tline = fgetl(fid);       %fprintf('\n%s',tline);
    end
    
    % read <rho_b>, <p_b>
    tval = fscanf(fid,'%f',2);
    rho_b = tval(1); p_b = tval(2);  %fprintf('\nrho_bf = %g\tp_b = %g',rho_b,p_b);
    tline = fgets(fid);
    
    % skip comments
    for j = 1:1:3
        tline = fgetl(fid);       %fprintf('\n%s',tline);
    end
    
    % read <gdfd>
    gdfd = fscanf(fid,'%f',1);    %fprintf('\ngdfd = %g',gdfd);
    tline = fgets(fid); 
    
    % skip comments
    for j = 1:1:3
        tline = fgetl(fid);       %fprintf('\n%s',tline);
    end
    
    % read <Ti>
    Ti = fscanf(fid,'%f',1);    %fprintf('\nTi = %g',Ti);
    tline = fgets(fid);
    
    % skip comments
    for j = 1:1:3
        tline = fgetl(fid);       %fprintf('\n%s',tline);
    end
    
    % read <Vi>
    Vi = fscanf(fid,'%f',1);    %fprintf('\nVi = %g',Vi);
    tline = fgets(fid);
    
    fclose(fid);
    
    if (  ( Ti < 0.0 ) || ( Ti > 0.5 )  )
        disp(['Error! Illegal value of Ti = ',num2str(Ti),'. Legal range is [0, 0.5].']);
        return
    end

    % self-control output
    fprintf('\nBackground Parameters:\n');
    fprintf('a1\t\t%g\n',a1);         fprintf('a2\t\t%g\n' , a2);
    fprintf('b0\t\t%g\n',b0);         fprintf('phi\t\t%g\t\t%g\n', phi, phi/180*pi);
    fprintf('f\t\t%g\n',f);           fprintf('k\t\t%g\n',k);      fprintf('n\t\t%g\n',n);
    fprintf('rho_b\t%g\n',rho_b);
    fprintf(  'p_b\t\t%g\n',  p_b);
    fprintf( 'gdfd\t%g\n', gdfd);
    fprintf( 'Ti\t\t%g\n', Ti);
    fprintf( 'Vi\t\t%g\n', Vi);
        
    phi = phi/180*pi;       % cast phi  to radians

    b1 = b0*cos(phi);       b2 = b0*sin(phi);
    
    BckParam = struct('a1',a1,'a2',a2,'b0',b0,'phi',phi, ...
                      'f',f,'k',k,'n',n,...
                      'rho_b',rho_b,'p_b',p_b,'gdfd',gdfd,'Ti',Ti,'Vi',Vi);
end
%%%-----------  end of the function <f_read_config_Kasko>   -------------

%%% calculate background plasma qunatities (Kasko) in 1 point
function VAL = f_bckgr_funcs_Kasko(dimq,x,z,BckParam)
%%% INPUT: 
%%% 1. dimq - length of vector of background plasma quantities.
%%% 2,3. (x,z) - coordinates of 1 point, numeric.
%%% 4. BckParam - structure, containing background parameters.
%%% OUTPUT: vector of background plasma quantities.

    a1 = BckParam.a1;         a2 = BckParam.a2;
    b0 = BckParam.b0;         phi = BckParam.phi;
    f  = BckParam.f;          k  = BckParam.k;      n  = BckParam.n;
       
    p_b = BckParam.p_b;   
    gdfd = BckParam.gdfd;
    rho_b = BckParam.rho_b;
    Ti = BckParam.Ti;
    Vi = BckParam.Vi;
    
    Mp = 1.6726219e-27;
    Me = 9.10938356e-31;
    Te = 0.5 - Ti;
    if ( Ti > 10.0*eps )       
        Vi = -abs(Vi);
        Ve = -Te/Ti*Vi;
    else
        fprintf('\nWarning: Ion temperature is zero, hence Vi = 0. Electron velocity is undefined.\n');
        Vi = 0.0;
        Ve = input('\nEnter electron velocity:   ');    Ve = +abs(Ve);
    end
    Vc = Vi - Ve;
    Vb = ( Mp*Vi + Me*Ve ) / (Mp+Me);

    sx = size(x);   % sz = size(z);
    if (  ( numel(x) == 1 ) && ( numel(z) == 1 )  ) 
        VAL = NaN*zeros(dimq,1,1);
        iz = 1; ix = 1;
    else
        VAL = NaN*zeros(dimq,sx(1),sx(2));
        iz = 1:sx(1);   ix = 1:sx(2);
    end
    
        
    if ( dimq > 00)        % VAL(01) = Psi
        VAL(01,iz,ix) = log((cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2).^(k./2)).*(f.^2 + 1).^(1./2) + f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2).^(k./2)))./(n.^2.*(x.^2 + z.^2).^(n - 1) + (b0.^2.*k.^2)./((a1 + x).^2 + (a2 + z).^2).^(k + 1) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)).^(1./2));
    end
    
    if ( dimq > 01)        % VAL(02) = Rho
        VAL(02,iz,ix) = rho_b + exp(-2.0.*VAL(01,iz,ix));
        %VAL(02,iz,ix) = rho_b + (n.^2.*((x.^2 + z.^2).^(1./2)).^(2.*n - 2) + (b0.^2.*k.^2)./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(2.*k + 2) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*((x.^2 + z.^2).^(1./2)).^(n - 1))./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(k + 1))./(f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)) + cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2)).^2;
    end
    
    if ( dimq > 02)        % VAL(03) = Pg     
        VAL(03,iz,ix) = p_b + 0.5.*exp(-2.0.*VAL(01,iz,ix));
        %VAL(03,iz,ix) = p_b + (n.^2.*((x.^2 + z.^2).^(1./2)).^(2.*n - 2) + (b0.^2.*k.^2)./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(2.*k + 2) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*((x.^2 + z.^2).^(1./2)).^(n - 1))./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(k + 1))./(2.*(f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)) + cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2)).^2);
    end    
    
    if ( dimq > 03)        % VAL(04) = Bx     
        VAL(04,iz,ix) = (((f.*sin(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(n.*z.*cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) - n.*x.*sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) + (b0.*k.*cos(phi - k.*atan((a2 + z)./(a1 + x))).*(2.*a2 + 2.*z))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1)) - (b0.*k.*sin(phi - k.*atan((a2 + z)./(a1 + x))).*(a1 + x))./(((a1 + x).^2 + (a2 + z).^2).^(k./2).*(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2))) - sinh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2).*(n.*x.*cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) + n.*z.*sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) + (b0.*k.*sin(phi - k.*atan((a2 + z)./(a1 + x))).*(2.*a2 + 2.*z))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1)) + (b0.*k.*cos(phi - k.*atan((a2 + z)./(a1 + x))).*(a1 + x))./(((a1 + x).^2 + (a2 + z).^2).^(k./2).*(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2))))./(n.^2.*(x.^2 + z.^2).^(n - 1) + (b0.^2.*k.^2)./((a1 + x).^2 + (a2 + z).^2).^(k + 1) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)).^(1./2) + ((cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2)).*(f.^2 + 1).^(1./2) + f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2))).*(n.^2.*z.*(x.^2 + z.^2).^(n - 2).*(2.*n - 2) - (b0.^2.*k.^2.*(2.*k + 2).*(2.*a2 + 2.*z))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k + 2)) + (2.*b0.*k.*n.*sin(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2).*((x.*(n - 1))./(x.^2 + z.^2) + ((a1 + x).*(k + 1))./(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2)))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2) - (b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2).*(k + 1).*(2.*a2 + 2.*z))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 3./2) + (2.*b0.*k.*n.*z.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 3./2).*(n - 1))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)))./(2.*(n.^2.*(x.^2 + z.^2).^(n - 1) + (b0.^2.*k.^2)./((a1 + x).^2 + (a2 + z).^2).^(k + 1) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)).^(3./2))).*(n.^2.*(x.^2 + z.^2).^(n - 1) + (b0.^2.*k.^2)./((a1 + x).^2 + (a2 + z).^2).^(k + 1) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)).^(1./2))./(cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2)).*(f.^2 + 1).^(1./2) + f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2)));
    end
    
    if ( dimq > 04)        % VAL(05) = By     
        VAL(05,iz,ix) = sym(gdfd);        
    end
    
    if ( dimq > 05)        % VAL(06) = Bz     
        VAL(06,iz,ix) = -(((f.*sin(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(n.*x.*cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) + n.*z.*sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) + (b0.*k.*cos(phi - k.*atan((a2 + z)./(a1 + x))).*(2.*a1 + 2.*x))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1)) + (b0.*k.*sin(phi - k.*atan((a2 + z)./(a1 + x))).*(a2 + z))./(((a1 + x).^2 + (a2 + z).^2).^(k./2).*(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2))) + sinh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2).*(n.*z.*cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) - n.*x.*sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) - (b0.*k.*sin(phi - k.*atan((a2 + z)./(a1 + x))).*(2.*a1 + 2.*x))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1)) + (b0.*k.*cos(phi - k.*atan((a2 + z)./(a1 + x))).*(a2 + z))./(((a1 + x).^2 + (a2 + z).^2).^(k./2).*(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2))))./(n.^2.*(x.^2 + z.^2).^(n - 1) + (b0.^2.*k.^2)./((a1 + x).^2 + (a2 + z).^2).^(k + 1) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)).^(1./2) - ((cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2)).*(f.^2 + 1).^(1./2) + f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2))).*((b0.^2.*k.^2.*(2.*k + 2).*(2.*a1 + 2.*x))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k + 2)) - n.^2.*x.*(x.^2 + z.^2).^(n - 2).*(2.*n - 2) + (2.*b0.*k.*n.*sin(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2).*((z.*(n - 1))./(x.^2 + z.^2) + ((a2 + z).*(k + 1))./(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2)))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2) + (b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2).*(k + 1).*(2.*a1 + 2.*x))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 3./2) - (2.*b0.*k.*n.*x.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 3./2).*(n - 1))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)))./(2.*(n.^2.*(x.^2 + z.^2).^(n - 1) + (b0.^2.*k.^2)./((a1 + x).^2 + (a2 + z).^2).^(k + 1) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)).^(3./2))).*(n.^2.*(x.^2 + z.^2).^(n - 1) + (b0.^2.*k.^2)./((a1 + x).^2 + (a2 + z).^2).^(k + 1) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)).^(1./2))./(cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2)).*(f.^2 + 1).^(1./2) + f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2)));
    end
    
    if ( dimq > 06)        % VAL(07) = Vx     
        VAL(07,iz,ix) = sym(0);        
    end
    
    if ( dimq > 07)        % VAL(08) = Vy
        VAL(08,iz,ix) = sym(Vb);
    end
    
    if ( dimq > 08)        % VAL(09) = Vz
        VAL(09,iz,ix) = sym(0);
    end
    
    if ( dimq > 09)        % VAL(10) = Jy
        VAL(10,iz,ix) = -exp(-2.0.*VAL(01,iz,ix));
        %VAL(10,iz,ix) = -(n.^2.*((x.^2 + z.^2).^(1./2)).^(2.*n - 2) + (b0.^2.*k.^2)./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(2.*k + 2) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*((x.^2 + z.^2).^(1./2)).^(n - 1))./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(k + 1))./(f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)) + cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2)).^2;
    end
    
    if ( dimq > 10)        % VAL(11) = Fx
        VAL(11,iz,ix) = ((b0.^2.*k.^2.*(2.*k + 2).*(2.*a1 + 2.*x))./(2.*((a1 + x).^2 + (a2 + z).^2).^(1./2).*(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(2.*k + 3)) - (n.^2.*x.*(2.*n - 2).*((x.^2 + z.^2).^(1./2)).^(2.*n - 3))./(x.^2 + z.^2).^(1./2) + (2.*b0.*k.*n.*sin(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*((x.^2 + z.^2).^(1./2)).^(n - 1).*(((a2 + z).*(k + 1))./((a1 + x).^2.*((a2 + z).^2./(a1 + x).^2 + 1)) + (z.*(n - 1))./(x.^2.*(z.^2./x.^2 + 1))))./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(k + 1) + (b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*((x.^2 + z.^2).^(1./2)).^(n - 1).*(k + 1).*(2.*a1 + 2.*x))./(((a1 + x).^2 + (a2 + z).^2).^(1./2).*(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(k + 2)) - (2.*b0.*k.*n.*x.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*((x.^2 + z.^2).^(1./2)).^(n - 2).*(n - 1))./((x.^2 + z.^2).^(1./2).*(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(k + 1)))./(2.*(f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)) + cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2)).^2) - ((f.*sin(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(n.*x.*cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) + (b0.*k.*cos(phi - k.*atan((a2 + z)./(a1 + x))).*(2.*a1 + 2.*x))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1)) + (n.*z.*sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2))./(x.^2.*(z.^2./x.^2 + 1)) + (b0.*k.*sin(phi - k.*atan((a2 + z)./(a1 + x))).*(a2 + z))./((a1 + x).^2.*((a1 + x).^2 + (a2 + z).^2).^(k./2).*((a2 + z).^2./(a1 + x).^2 + 1))) - sinh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2).*(n.*x.*sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) + (b0.*k.*sin(phi - k.*atan((a2 + z)./(a1 + x))).*(2.*a1 + 2.*x))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1)) - (n.*z.*cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2))./(x.^2.*(z.^2./x.^2 + 1)) - (b0.*k.*cos(phi - k.*atan((a2 + z)./(a1 + x))).*(a2 + z))./((a1 + x).^2.*((a1 + x).^2 + (a2 + z).^2).^(k./2).*((a2 + z).^2./(a1 + x).^2 + 1)))).*(n.^2.*((x.^2 + z.^2).^(1./2)).^(2.*n - 2) + (b0.^2.*k.^2)./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(2.*k + 2) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*((x.^2 + z.^2).^(1./2)).^(n - 1))./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(k + 1)))./(f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)) + cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2)).^3 + (((f.*sin(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(n.*x.*cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) + n.*z.*sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) + (b0.*k.*cos(phi - k.*atan((a2 + z)./(a1 + x))).*(2.*a1 + 2.*x))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1)) + (b0.*k.*sin(phi - k.*atan((a2 + z)./(a1 + x))).*(a2 + z))./(((a1 + x).^2 + (a2 + z).^2).^(k./2).*(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2))) + sinh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2).*(n.*z.*cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) - n.*x.*sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) - (b0.*k.*sin(phi - k.*atan((a2 + z)./(a1 + x))).*(2.*a1 + 2.*x))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1)) + (b0.*k.*cos(phi - k.*atan((a2 + z)./(a1 + x))).*(a2 + z))./(((a1 + x).^2 + (a2 + z).^2).^(k./2).*(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2))))./(n.^2.*(x.^2 + z.^2).^(n - 1) + (b0.^2.*k.^2)./((a1 + x).^2 + (a2 + z).^2).^(k + 1) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)).^(1./2) - ((cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2)).*(f.^2 + 1).^(1./2) + f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2))).*((b0.^2.*k.^2.*(2.*k + 2).*(2.*a1 + 2.*x))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k + 2)) - n.^2.*x.*(x.^2 + z.^2).^(n - 2).*(2.*n - 2) + (2.*b0.*k.*n.*sin(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2).*((z.*(n - 1))./(x.^2 + z.^2) + ((a2 + z).*(k + 1))./(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2)))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2) + (b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2).*(k + 1).*(2.*a1 + 2.*x))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 3./2) - (2.*b0.*k.*n.*x.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 3./2).*(n - 1))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)))./(2.*(n.^2.*(x.^2 + z.^2).^(n - 1) + (b0.^2.*k.^2)./((a1 + x).^2 + (a2 + z).^2).^(k + 1) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)).^(3./2))).*(n.^2.*(x.^2 + z.^2).^(n - 1) + (b0.^2.*k.^2)./((a1 + x).^2 + (a2 + z).^2).^(k + 1) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)).^(1./2).*(n.^2.*((x.^2 + z.^2).^(1./2)).^(2.*n - 2) + (b0.^2.*k.^2)./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(2.*k + 2) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*((x.^2 + z.^2).^(1./2)).^(n - 1))./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(k + 1)))./((cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2)).*(f.^2 + 1).^(1./2) + f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2))).*(f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)) + cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2)).^2);
    end
    
    if ( dimq > 11)        % VAL(12) = Fz
        VAL(12,iz,ix) = (((f.*sin(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(n.*z.*cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) - n.*x.*sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) + (b0.*k.*cos(phi - k.*atan((a2 + z)./(a1 + x))).*(2.*a2 + 2.*z))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1)) - (b0.*k.*sin(phi - k.*atan((a2 + z)./(a1 + x))).*(a1 + x))./(((a1 + x).^2 + (a2 + z).^2).^(k./2).*(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2))) - sinh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2).*(n.*x.*cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) + n.*z.*sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) + (b0.*k.*sin(phi - k.*atan((a2 + z)./(a1 + x))).*(2.*a2 + 2.*z))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1)) + (b0.*k.*cos(phi - k.*atan((a2 + z)./(a1 + x))).*(a1 + x))./(((a1 + x).^2 + (a2 + z).^2).^(k./2).*(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2))))./(n.^2.*(x.^2 + z.^2).^(n - 1) + (b0.^2.*k.^2)./((a1 + x).^2 + (a2 + z).^2).^(k + 1) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)).^(1./2) + ((cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2)).*(f.^2 + 1).^(1./2) + f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2))).*(n.^2.*z.*(x.^2 + z.^2).^(n - 2).*(2.*n - 2) - (b0.^2.*k.^2.*(2.*k + 2).*(2.*a2 + 2.*z))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k + 2)) + (2.*b0.*k.*n.*sin(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2).*((x.*(n - 1))./(x.^2 + z.^2) + ((a1 + x).*(k + 1))./(a1.^2 + 2.*a1.*x + a2.^2 + 2.*a2.*z + x.^2 + z.^2)))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2) - (b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2).*(k + 1).*(2.*a2 + 2.*z))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 3./2) + (2.*b0.*k.*n.*z.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 3./2).*(n - 1))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)))./(2.*(n.^2.*(x.^2 + z.^2).^(n - 1) + (b0.^2.*k.^2)./((a1 + x).^2 + (a2 + z).^2).^(k + 1) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)).^(3./2))).*(n.^2.*(x.^2 + z.^2).^(n - 1) + (b0.^2.*k.^2)./((a1 + x).^2 + (a2 + z).^2).^(k + 1) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*(x.^2 + z.^2).^(n./2 - 1./2))./((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1./2)).^(1./2).*(n.^2.*((x.^2 + z.^2).^(1./2)).^(2.*n - 2) + (b0.^2.*k.^2)./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(2.*k + 2) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*((x.^2 + z.^2).^(1./2)).^(n - 1))./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(k + 1)))./((cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2)).*(f.^2 + 1).^(1./2) + f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./(2.*a1.*x + 2.*a2.*z + a1.^2 + a2.^2 + x.^2 + z.^2).^(k./2))).*(f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)) + cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2)).^2) - ((f.*sin(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(n.*z.*cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) - (n.*sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2))./(x.*(z.^2./x.^2 + 1)) + (b0.*k.*cos(phi - k.*atan((a2 + z)./(a1 + x))).*(2.*a2 + 2.*z))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1)) - (b0.*k.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).*((a1 + x).^2 + (a2 + z).^2).^(k./2).*((a2 + z).^2./(a1 + x).^2 + 1))) - sinh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2).*(n.*z.*sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2 - 1) + (b0.*k.*sin(phi - k.*atan((a2 + z)./(a1 + x))).*(2.*a2 + 2.*z))./(2.*((a1 + x).^2 + (a2 + z).^2).^(k./2 + 1)) + (n.*cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2))./(x.*(z.^2./x.^2 + 1)) + (b0.*k.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).*((a1 + x).^2 + (a2 + z).^2).^(k./2).*((a2 + z).^2./(a1 + x).^2 + 1)))).*(n.^2.*((x.^2 + z.^2).^(1./2)).^(2.*n - 2) + (b0.^2.*k.^2)./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(2.*k + 2) + (2.*b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*((x.^2 + z.^2).^(1./2)).^(n - 1))./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(k + 1)))./(f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)) + cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2)).^3 - ((n.^2.*z.*(2.*n - 2).*((x.^2 + z.^2).^(1./2)).^(2.*n - 3))./(x.^2 + z.^2).^(1./2) - (b0.^2.*k.^2.*(2.*k + 2).*(2.*a2 + 2.*z))./(2.*((a1 + x).^2 + (a2 + z).^2).^(1./2).*(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(2.*k + 3)) + (2.*b0.*k.*n.*sin(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*((k + 1)./((a1 + x).*((a2 + z).^2./(a1 + x).^2 + 1)) + (n - 1)./(x.*(z.^2./x.^2 + 1))).*((x.^2 + z.^2).^(1./2)).^(n - 1))./(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(k + 1) - (b0.*k.*n.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*((x.^2 + z.^2).^(1./2)).^(n - 1).*(k + 1).*(2.*a2 + 2.*z))./(((a1 + x).^2 + (a2 + z).^2).^(1./2).*(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(k + 2)) + (2.*b0.*k.*n.*z.*cos(phi - atan(z./x).*(n - 1) - atan((a2 + z)./(a1 + x)).*(k + 1)).*((x.^2 + z.^2).^(1./2)).^(n - 2).*(n - 1))./((x.^2 + z.^2).^(1./2).*(((a1 + x).^2 + (a2 + z).^2).^(1./2)).^(k + 1)))./(2.*(f.*cos(cos(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*cos(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)) + cosh(sin(n.*atan(z./x)).*(x.^2 + z.^2).^(n./2) - (b0.*sin(phi - k.*atan((a2 + z)./(a1 + x))))./((a1 + x).^2 + (a2 + z).^2).^(k./2)).*(f.^2 + 1).^(1./2)).^2);
    end
    
    if (dimq > 12)
        fprintf('WARNING! Number of calculated quantities is less then requested!\n');
    end
end
%%%-----------  end of the function <f_bckgr_funcs_Kasko>   -------------

%%% save data in single grid-format file 
function f_grid_file_write(x, y, z, val, file, format)
%%% INPUT:
%%% 1,2,3. (x,y,z) - coordinates of the grid points.
%%% 4. val - array of quatities. size = [dimq,dimx,dimy,dimz].
%%% 5. file - file name, string.
%%% 6. format - {'ascii','binary'}, file type.

    sysdim = size(val, 1);
    spacedim = ~isempty(x) + ~isempty(y) + ~isempty(z);

    fid = fopen(file, 'w');

    if strcmp(format, 'ascii')

        fprintf(fid, '%d\n', sysdim);

        fprintf(fid, '%d\n', spacedim);

        if ~isempty(x), fprintf(fid, '%d\n', length(x)); end
        if ~isempty(y), fprintf(fid, '%d\n', length(y)); end
        if ~isempty(z), fprintf(fid, '%d\n', length(z)); end

        if ~isempty(x), fprintf(fid, '%.15e\n', x); end
        if ~isempty(y), fprintf(fid, '%.15e\n', y); end
        if ~isempty(z), fprintf(fid, '%.15e\n', z); end

        frm = ''; for k=1:sysdim, frm = strcat(frm,'%.15e\t'); end; frm = strcat(frm,'\n');

        fprintf(fid, frm, val);

    else

        fwrite(fid, sysdim, 'int');
        fwrite(fid, spacedim, 'int');

        if ~isempty(x),     fwrite(fid, length(x), 'int');  end,
        if ~isempty(y),     fwrite(fid, length(y), 'int');  end,
        if ~isempty(z),     fwrite(fid, length(z), 'int');  end,

        if ~isempty(x),     fwrite(fid, x, 'double');       end,
        if ~isempty(y),     fwrite(fid, y, 'double');       end,
        if ~isempty(z),     fwrite(fid, z, 'double');       end,

        fwrite(fid, val, 'double');
        
    end
    fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%  T H E   E N D   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%