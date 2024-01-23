% maketissue_18apr17.m
% maketissue.m
%   Creates a cube of optical property pointers,T(y,x,z), saved in
%       myname_T.bin = a tissue structure file
%   which specifies a complex tissue for use by mcxyz.c.
%
%   Also prepares a listing of the optical properties at chosen wavelength
%   for use by mcxyz.c, [mua, mus, g], for each tissue type specified
%   in myname_T.bin. This listing is saved in
%       myname_H.mci = the input file for use by mcxyz.c.
%
%   Will generate a figure illustrating the tissue with its various
%   tissue types and the beam being launched.
%
%   Uses
%       makeTissueList.m
%
%   To use, 
%       1. Prepare makeTissueList.m so that it contains the tissue
%   types desired.
%       2. Specify the USER CHOICES.
%       2. Run this program, maketissue.m.
%
%   Note: mcxyz.c can use optical properties in cm^-1 or mm^-1 or m^-1,
%       if the bin size (binsize) is specified in cm or mm or m,
%       respectively.
%
%  Steven L. Jacques. updated Aug 21, 2014.
%       

clear
format compact
clc
home

%%% USER CHOICES %%%%%%%% <-------- You must set these parameters ------
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci 
                        % 0 = don't save. Just check the program.

myname      = 'test';% name for files: myname_T.bin, myname_H.mci  
time_min    = 10;      	% time duration of the simulation [min] <----- run time -----
nm          = 450;   	% desired wavelength of simulation
Nbins       = 400;    	% # of bins in each dimension of cube 
binsize     = 0.00264; 	% size of each bin, eg. [cm] or [mm]

% Set Monte Carlo launch flags
mcflag      = 0;     	% launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt. 
                        % 3 = rectangular beam (use xfocus,yfocus for x,y halfwidths)
launchflag  = 0;        % 0 = let mcxyz.c calculate launch trajectory
                        % 1 = manually set launch vector.
boundaryflag = 1;       % 0 = no boundaries, 1 = escape at boundaries
                        % 2 = escape at surface only. No x, y, bottom z
                        % boundaries

% Sets position of source
xs          = 0;      	% x of source
ys          = 0;        % y of source
zs          = 0.12276;  	% z of source 

% Set position of focus, so mcxyz can calculate launch trajectory
xfocus      = 0;        % set x,position of focus
yfocus      = 0;        % set y,position of focus
zfocus      = inf;    	% set z,position of focus (=inf for collimated beam)

% only used if mcflag == 0 or 1 or 3 (not 2=isotropic pt.)
radius      = 0.0500;   % 1/e radius of beam at tissue surface
waist       = 0.0500;  	% 1/e radius of beam at focus

% only used if launchflag == 1 (manually set launch trajectory):
ux0         = 0;      % trajectory projected onto x axis
uy0         = 0;      % trajectory projected onto y axis
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% 
% Prepare Monte Carlo 
%%%%%%%%%%

% Create tissue properties
tissue = makeTissueList(nm); % also --> global tissue(1:Nt).s

Nt = length(tissue);
for i=1:Nt
    muav(i)  = tissue(i).mua;
    musv(i)  = tissue(i).mus;
    gv(i)    = tissue(i).g;
end

% Specify Monte Carlo parameters    
Nx = Nbins;
Ny = Nbins;
Nz = Nbins;
dx = binsize;
dy = binsize;
dz = binsize;
x  = ([1:Nx]'-Nx/2)*dx;
y  = ([1:Ny]'-Ny/2)*dy;
z  = [1:Nz]'*dz;
zmin = min(z);
zmax = max(z);
xmin = min(x);
xmax = max(x);

if isinf(zfocus), zfocus = 1e12; end

%%%%%%
% CREATE TISSUE STRUCTURE T(y,x,z)
%   Create T(y,x,z) by specifying a tissue type (an integer)
%   for each voxel in T.
%
%   Note: one need not use every tissue type in the tissue list.
%   The tissue list is a library of possible tissue types.

T = double(zeros(Ny,Nx,Nz)); 

T = T + 2; %start by filling entire model with CSF


%Pig T10 Dimensions
TS_DV = 0.26; %thecal sac dorsal ventral radius (cm) %0.378
TS_Width = 0.38; %thecal sac transverse radius (cm) %0.369
dura_thickness = 0.027; %thickness in cm


SC_DV = 0.280; %spinal cord dorsal ventral radius(cm)
SC_Width = 0.342; % spinal cord transverse radius (cm)
AP_ecc = 0.104; %eccentrecity of spinal cord within the canal. ranges from -1 to 1. 0 is prefectly centered

%{ 
Human T10 Dimensions
% Declare parameters for spinal cord geometry (all units in cm). 
TS_DV = 0.630;%;  %thecal sac dorsal-ventral (or anteroposterior radius);
TS_APTR = 0.796; % anteroposterior to transverse width ratio for thecal sac.
TS_Width = TS_DV/TS_APTR; %thecal sac transverse radius
dura_thickness = 0.031;

SC_DV = 0.33;  %spinal cord dorsal-ventral radius
SC_Width = 0.42; %spinal cord transverse radius
AP_ecc = 0.385; %positive means cord is close to dorsal side. %should be 0.385
%} 

Dura_DV = TS_DV + dura_thickness; %dura DV radius
Dura_Width = TS_Width + dura_thickness; %Dura TV radius

%alright, first lets start by loading an image that dictates the placementcof gray matter
gray_map = imread('Lession.png'); 

%now, lets figure out where the gray_map starts in relation to the tissue model
xc=0; zc= 0.528; %specified center position of dural sac
zcord = zc + (TS_DV-SC_DV)*AP_ecc; %center of spinal cord
xc_idx = Nx/2; %assumes xc is 0.

zc_idx = round(zcord/binsize);
half_w = round(size(gray_map,2)/2);
half_h = round(size(gray_map,1)/2);
z_start_idx = round((Nz - size(gray_map,2)/2 + (zcord-SC_DV)/binsize)); 

for jz=1:size(gray_map,1)
    for jx=1:size(gray_map,2)
        if gray_map(jz,jx,1) > 100 %if red channel high: white
            T(:,xc_idx - half_w + jx, zc_idx - half_h +jz) = 4; %4 dictates white matter
        end
        if gray_map(jz,jx,2) > 100 %if green channel high: gray matter
            T(:,xc_idx - half_w + jx, zc_idx - half_h +jz) = 5; %5 dictates gray matter
        end
        if gray_map(jz,jx,3) > 100 %if blue channel high: hydrogel
            T(:,xc_idx - half_w + jx, zc_idx - half_h +jz) = 6; %6 dictates Hydrogel
        end
    end
end


%now we have the gray matter and white matter set. Everything outside cord
%diameter, lets set as CSF, then Dura, then air.

for iz=1:Nz % for every depth z(iz)

    for ix=1:Nx
        
        %calculate if outside spinal cord - then fill with csf
        xd = x(ix) - xc; zd = z(iz) - zcord;
        r = sqrt(xd^2/SC_Width^2 + zd^2/SC_DV^2); %equation for ellipse
        if (r > 1)   %fill everything outside cord as CSF
                T(:,ix,iz) = 2; %2 dictates CSF
        end

        %calculate inside boundary of dura mater - fill everything outside with dura
        zd = z(iz) - zc;
        r = sqrt((xd)^2/TS_Width^2 + (zd)^2/TS_DV^2); %equation for ellipse
        if (r > 1)   %fill with Dura
                T(:,ix,iz) = 3;
        end

        %calculate if outside dura - then fill with air
        r = sqrt(xd^2/Dura_Width^2 + zd^2/Dura_DV^2); %equation for ellipse
        if (r > 1)  
                T(:,ix,iz) = 1;
        end
        
        
    end

    
end % iz


%%
if SAVEON
    tic
    % convert T to linear array of integer values, v(i)i = 0;
    v = uint8(reshape(T,Ny*Nx*Nz,1));

    %% WRITE FILES
    % Write myname_H.mci file
    %   which contains the Monte Carlo simulation parameters
    %   and specifies the tissue optical properties for each tissue type.
    commandwindow
    disp(sprintf('--------create %s --------',myname))
    filename = sprintf('%s_H.mci',myname);
    fid = fopen(filename,'w');
        % run parameters
        fprintf(fid,'%0.2f\n',time_min);
        fprintf(fid,'%d\n'   ,Nx);
        fprintf(fid,'%d\n'   ,Ny);
        fprintf(fid,'%d\n'   ,Nz);
        fprintf(fid,'%0.4f\n',dx);
        fprintf(fid,'%0.4f\n',dy);
        fprintf(fid,'%0.4f\n',dz);
        % launch parameters
        fprintf(fid,'%d\n'   ,mcflag);
        fprintf(fid,'%d\n'   ,launchflag);
        fprintf(fid,'%d\n'   ,boundaryflag);
        fprintf(fid,'%0.4f\n',xs);
        fprintf(fid,'%0.4f\n',ys);
        fprintf(fid,'%0.4f\n',zs);
        fprintf(fid,'%0.4f\n',xfocus);
        fprintf(fid,'%0.4f\n',yfocus);
        fprintf(fid,'%0.4f\n',zfocus);
        fprintf(fid,'%0.4f\n',ux0); % if manually setting ux,uy,uz
        fprintf(fid,'%0.4f\n',uy0);
        fprintf(fid,'%0.4f\n',uz0);
        fprintf(fid,'%0.4f\n',radius);
        fprintf(fid,'%0.4f\n',waist);
        % tissue optical properties
        fprintf(fid,'%d\n',Nt);
        for i=1:Nt
            fprintf(fid,'%0.4f\n',muav(i));
            fprintf(fid,'%0.4f\n',musv(i));
            fprintf(fid,'%0.4f\n',gv(i));
        end
    fclose(fid);

    %% write myname_T.bin file
    filename = sprintf('%s_T.bin',myname);
    disp(['create ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,v,'uint8');
    fclose(fid);

    toc
end % SAVEON


%% Look at structure of Tzx at iy=Ny/2
Txzy = shiftdim(T,1);   % Tyxz --> Txzy
Tzx  = Txzy(:,:,Ny/2)'; % Tzx


%%
figure(1); clf
sz = 12;  fz = 10; 
imagesc(x,z,Tzx,[1 Nt])
hold on
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
colorbar
cmap = makecmap(Nt);
colormap(cmap)
set(colorbar,'fontsize',1)
% label colorbar
zdiff = zmax-zmin;
%%%

for i=1:Nt
    yy = (Nt-i)/(Nt-1)*Nz*dz;
    text(max(x)*1.2,yy, tissue(i).name,'fontsize',fz)
end

text(xmax,zmin - zdiff*0.06, 'Tissue types','fontsize',fz)
axis equal image
axis([xmin xmax zmin zmax])
% 
%%% draw launch
N = 7; % # of beam rays drawn
switch mcflag
    case 0 % uniform
        for i=0:N
            plot((-radius + 2*radius*i/N)*[1 1],[zs max(z)],'r-')
        end

    case 1 % Gaussian
        for i=0:N
            plot([(-radius + 2*radius*i/N) xfocus],[zs zfocus],'r-')
        end

    case 2 % iso-point
        for i=1:N
            th = (i-1)/19*2*pi;
            xx = Nx/2*cos(th) + xs;
            zz = Nx/2*sin(th) + zs;
            plot([xs xx],[zs zz],'r-')
        end
        
    case 3 % rectangle
        zz = max(z);
        for i=1:N
            xx = -radius + 2*radius*i/20;
            plot([xx xx],[zs zz],'r-')
        end
end

disp('done')

