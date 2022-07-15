% FDTD 3D Cavity with CPML (Measuring S-Parameters)
% Author:shayan dodge
% Email address:dodgeshayan@gmail.com
%% initialize the matlab workspace
close all;
clear all;
clc;
complex_number=j;
%% frequency domain parameters
frequency_domain.start=2.85e9;
frequency_domain.end=3.3e9;
frequency_domain.gam=6400;
frequency_domain.step=(frequency_domain.end-frequency_domain.start)/frequency_domain.gam;
FER=[frequency_domain.start:frequency_domain.step:frequency_domain.end];
%% some constants
mu_0 = 1.2566370614359173e-06;
eps_0= 8.8541878176203892e-12;
c=299792458.0;% speed of light
%% wave definition
waveforms.cosine_modulated_gaussian(1).bandwith=4e9;%0.55
waveforms.cosine_modulated_gaussian(1).modulation_frequency=frequency_domain.end;%+(frequency_domain.end-frequency_domain.start)/2;
T=1/waveforms.cosine_modulated_gaussian(1).modulation_frequency;
amptidute=100;
lambda=(c*T);
omega=2*pi*waveforms.cosine_modulated_gaussian(1).modulation_frequency;
%% FDTD variables
number_of_cells_per_wavelength=20;%.*(frequency_domain.start./index);
dx=0.9*lambda/number_of_cells_per_wavelength;
dy=0.9*lambda/number_of_cells_per_wavelength;
dz=0.9.*lambda/number_of_cells_per_wavelength;
totalTime=500*T;
courant_factor=1;
dt=1/(c*sqrt((1/dx^2)+(1/dy^2)+(1/dz^2)));
dt=courant_factor*dt;
totalTimeStep=floor(totalTime/dt);
number_of_time_steps=500000;
%% Waveguide parameters
a=100.*10^-3;
b=100.*10^-3;
L=200.*10^-3;
t=(0.1*L)+1*(dx+dy);
%% boundary conditions
boundary.type_xn='cpml';
boundary.air_buffer_number_of_cells_xn=0; 
boundary.cpml_number_of_cells_xn=-ceil(0.1*L/dx);

boundary.type_xp = 'cpml';
boundary.air_buffer_number_of_cells_xp=0;
boundary.cpml_number_of_cells_xp=-ceil(0.1*L/dx);

boundary.type_yn = 'cpml';
boundary.air_buffer_number_of_cells_yn=0;
boundary.cpml_number_of_cells_yn=-ceil(0.1*L/dy);

boundary.type_yp = 'cpml';
boundary.air_buffer_number_of_cells_yp=0;
boundary.cpml_number_of_cells_yp=-ceil(0.1*L/dy);

boundary.type_zn = 'cpml';
boundary.air_buffer_number_of_cells_zn=ceil(0.5*L/dz);
boundary.cpml_number_of_cells_zn=ceil(0.25*L/dz);

boundary.type_zp='cpml';
boundary.air_buffer_number_of_cells_zp=0;
boundary.cpml_number_of_cells_zp=-ceil(0.25*L/dz);

boundary.cpml_order = 4;
boundary.cpml_sigma_max = 1;
boundary.cpml_kappa_max = 15;
boundary.cpml_alpha_order = 1; 
boundary.cpml_alpha_max = 0.24;
boundary.cpml_eps_R= 1;
%% materialtype
%here define and initialize the arrays of material types
%air
material_type(1).eps_r=1;
material_type(1).mu_r=1;
material_type(1).sigma_e=0;
material_type(1).sigma_m=0;
material_type(1).color=[1 1 1];
%a dielectric_1
material_type(2).eps_r=1;
material_type(2).mu_r=1;
material_type(2).sigma_e=omega.*eps_0*0.1;
material_type(2).sigma_m=omega.*mu_0*0.1;
material_type(2).color=[1 1 1];
%a pec
material_type(3).eps_r=1;
material_type(3).mu_r=1;
material_type(3).sigma_e=10^20;
material_type(3).sigma_m=0;
material_type(3).color=[1 1 1];
%indices of material types defining air, pec, and pmc
material_type_index_air=1;
material_type_index_dielectric_1=2;
material_type_index_pec=3;
%% define_geometry
%define a brick
% %y direction
brick(1).min_x=-((a./2)+(t));
brick(1).min_y= b./2;
brick(1).min_z= -0.5*L;
brick(1).max_x= ((a./2)+(t));
brick(1).max_y= ((b./2)+(t));
brick(1).max_z= ((1.75.*L));
brick(1).material_type=3;

brick(2).min_x=-((a./2)+(t));
brick(2).min_y=-((b./2)+(t));
brick(2).min_z= -0.5*L;
brick(2).max_x= ((a./2)+(t));
brick(2).max_y=-b./2;
brick(2).max_z=((1.75.*L));
brick(2).material_type=3;
%x direction
brick(3).min_x=-((a./2)+(t));
brick(3).min_y=-((b./2)+(t));
brick(3).min_z=-0.5*L;
brick(3).max_x=-(a./2);
brick(3).max_y= ((b./2)+(t));
brick(3).max_z=((1.75.*L));
brick(3).material_type=3;

brick(4).min_x= a./2;
brick(4).min_y=-((b./2)+(t));
brick(4).min_z=-0.5*L;
brick(4).max_x= ((a./2)+(t));
brick(4).max_y= ((b./2)+(t));
brick(4).max_z=((1.75.*L));
brick(4).material_type=3;

%% calculating the domain size
number_of_brick=size(brick,2);
%find the minimum and maximum coordinates of a box encapsulating the object
number_of_objects=1;

for i=1:number_of_brick
    min_x(number_of_objects)=brick(i).min_x;
    min_y(number_of_objects)=brick(i).min_y;
    min_z(number_of_objects)=brick(i).min_z;
    max_x(number_of_objects)=brick(i).max_x;
    max_y(number_of_objects)=brick(i).max_y; 
    max_z(number_of_objects)=brick(i).max_z;
    number_of_objects=number_of_objects+1;
end
fdtd_domain.min_x=min(min_x);
fdtd_domain.min_y=min(min_y);    
fdtd_domain.min_z=min(min_z);    
fdtd_domain.max_x=max(max_x);    
fdtd_domain.max_y=max(max_y);    
fdtd_domain.max_z=max(max_z);     
% Determine the problem space boundaries including air buffers
fdtd_domain.min_x = fdtd_domain.min_x-dx *...
    boundary.air_buffer_number_of_cells_xn;
fdtd_domain.min_y = fdtd_domain.min_y-dy *...
    boundary.air_buffer_number_of_cells_yn ;
fdtd_domain.min_z = fdtd_domain.min_z-dz*...
    boundary.air_buffer_number_of_cells_zn ; 
fdtd_domain.max_x = fdtd_domain.max_x+dx *...
    boundary.air_buffer_number_of_cells_xp; 
fdtd_domain.max_y = fdtd_domain. max_y+dy *...
    boundary.air_buffer_number_of_cells_yp ;
fdtd_domain. max_z= fdtd_domain. max_z+dz *...
    boundary.air_buffer_number_of_cells_zp;
% Determine the problem space boundaries including cpml layers
if strcmp (boundary.type_xn,'cpml') &&(boundary.cpml_number_of_cells_xn>0)
    fdtd_domain.min_x = fdtd_domain.min_x -dx *...
        boundary.cpml_number_of_cells_xn;
end
if strcmp (boundary.type_xp, 'cpml') &&(boundary.cpml_number_of_cells_xp >0)
    fdtd_domain.max_x = fdtd_domain. max_x + dx *...
        boundary.cpml_number_of_cells_xp;
end
if strcmp( boundary.type_yn, 'cpml') &&(boundary.cpml_number_of_cells_yn >0)
    fdtd_domain .min_y = fdtd_domain. min_y - dy *...
        boundary.cpml_number_of_cells_yn;
end
if strcmp (boundary.type_yp, 'cpml') &&(boundary.cpml_number_of_cells_yp >0)
   fdtd_domain. max_y = fdtd_domain. max_y+ dy *...
       boundary.cpml_number_of_cells_yp ;
end
if strcmp( boundary.type_zn, 'cpml') &&(boundary.cpml_number_of_cells_zn >0)
    fdtd_domain. min_z = fdtd_domain.min_z - dz *...
        boundary.cpml_number_of_cells_zn ;
end
if strcmp (boundary.type_zp, 'cpml') && (boundary.cpml_number_of_cells_zp>0) 
    fdtd_domain. max_z = fdtd_domain. max_z + dz*...
        boundary.cpml_number_of_cells_zp;
end
%detemining the problem space size
fdtd_domain.size_x=fdtd_domain.max_x-fdtd_domain.min_x;
fdtd_domain.size_y=fdtd_domain.max_y-fdtd_domain.min_y;
fdtd_domain.size_z=fdtd_domain.max_z-fdtd_domain.min_z;
%number of cells in x, y, and z directions
nx=round(fdtd_domain.size_x/dx);
ny=round(fdtd_domain.size_y/dy);
nz=round(fdtd_domain.size_z/dz);
%adjust domain size by snapping to cells
fdtd_domain.size_x=nx*dx;
fdtd_domain.size_y=ny*dy;
fdtd_domain.size_z=nz*dz;
fdtd_domain.max_x=fdtd_domain.min_x+fdtd_domain.size_x;
fdtd_domain.max_y=fdtd_domain.min_y+fdtd_domain.size_y;
fdtd_domain.max_z=fdtd_domain.min_z+fdtd_domain.size_z;
%some frequently used auxiliary parametrs
nxp1=nx+1;  nyp1=ny+1;  nzp1=nz+1;
nxm1=nx-1;  nym1=ny-1;  nzm1=nz-1;
nxm2=nx-2;  nym2=ny-2;  nzm2=nz-2;
%create arrays storing the center coordinates of the cells in
fdtd_domain.cell_center_coordinates_x=zeros(nx,ny,nz);
fdtd_domain.cell_center_coordinates_y=zeros(nx,ny,nz);
fdtd_domain.cell_center_coordinates_z=zeros(nx,ny,nz);
for ind =1:nx
    fdtd_domain.cell_center_coordinates_x(ind,:,:)=...
        (ind-0.5)*dx+fdtd_domain.min_x;
end
for ind =1:ny
    fdtd_domain.cell_center_coordinates_y(:,ind,:)=...
        (ind-0.5)*dy+fdtd_domain.min_y;
end
for ind =1:nz
    fdtd_domain.cell_center_coordinates_z(:,:,ind)=...
        (ind-0.5)*dz+fdtd_domain.min_z;
end
xcoor=linspace (fdtd_domain. min_x, fdtd_domain.max_x, nxp1);
ycoor=linspace (fdtd_domain. min_y, fdtd_domain.max_y, nyp1);
zcoor=linspace (fdtd_domain. min_z, fdtd_domain.max_z, nzp1);
%% material_3d_space
material_3d_space=ones(nx,ny,nz);
%% creating_brick
%creat the 3d object in problem space by
for ind=1:number_of_brick
%%convert brick end coordinates to node indices
blx=round((brick(ind).min_x-fdtd_domain.min_x)/dx)+1;
bly=round((brick(ind).min_y-fdtd_domain.min_y)/dy)+1;
blz=round((brick(ind).min_z-fdtd_domain.min_z)/dz)+1;
bux=round((brick(ind).max_x-fdtd_domain.min_x)/dx)+1;
buy=round((brick(ind).max_y-fdtd_domain.min_y)/dy)+1;
buz=round((brick(ind).max_z-fdtd_domain.min_z)/dz)+1;
%%assign material type of brick to the cells
material_3d_space(blx:bux-1,bly:buy-1,blz:buz-1)=brick(ind).material_type;
end
cx=fdtd_domain.cell_center_coordinates_x;
cy=fdtd_domain.cell_center_coordinates_y;
cz=fdtd_domain.cell_center_coordinates_z;
distance=zeros(nx,ny,nz);


for ind=number_of_brick
%%convert brick end coordinates to node indices
blx=round((brick(ind).min_x-fdtd_domain.min_x)/dx)+1;
bly=round((brick(ind).min_y-fdtd_domain.min_y)/dy)+1;
blz=round((brick(ind).min_z-fdtd_domain.min_z)/dz)+1;
bux=round((brick(ind).max_x-fdtd_domain.min_x)/dx)+1;
buy=round((brick(ind).max_y-fdtd_domain.min_y)/dy)+1;
buz=round((brick(ind).max_z-fdtd_domain.min_z)/dz)+1;
%%assign material type of brick to the cells
material_3d_space(blx:bux-1,bly:buy-1,blz:buz-1)=brick(ind).material_type;
end

eps_r_x=ones(nx,nyp1,nzp1);
eps_r_y=ones(nxp1,ny,nzp1);
eps_r_z=ones(nxp1,nyp1,nz);
mu_r_x=ones(nxp1,ny,nz);
mu_r_y=ones(nx,nyp1,nz);
mu_r_z=ones(nx,ny,nzp1);
sigma_e_x=zeros(nx,nyp1,nzp1);
sigma_e_y=zeros(nxp1,ny,nzp1);
sigma_e_z=zeros(nxp1,nyp1,nz);
sigma_m_x=zeros(nxp1,ny,nz);
sigma_m_y=zeros(nx,nyp1,nz);
sigma_m_z=zeros(nx,ny,nzp1);
%% calculate_material_component_values
%calculate material component values by averaging
for ind=1:size(material_type,2)
    t_eps_r(ind)=material_type(ind).eps_r;
    t_mu_r(ind)=material_type(ind).mu_r;
    t_sigma_e(ind)=material_type(ind).sigma_e;
    t_sigma_m(ind)=material_type(ind).sigma_m;
end
%assign negligibly small values to t_mu_r and t_sigma_m where they are zero
%in order to prevent division by zerro error
t_mu_r(find(t_mu_r==0))=1e-20;
t_sigma_m(find(t_sigma_m==0.0000))=1e-20;
disp('calculating_eps_r_x');
%eps_r_x(i,j,k)is average of four cells(i,j,k)(i,j-1,k)(i,j,k-1)(i,j-1,k-1)
eps_r_x(1:nx,2:ny,2:nz)=0.25*(t_eps_r(material_3d_space(1:nx,2:ny,2:nz))+...
    t_eps_r(material_3d_space(1:nx,1:ny-1,2:nz))+...
    t_eps_r(material_3d_space(1:nx,2:ny,1:nz-1))+...
    t_eps_r(material_3d_space(1:nx,1:ny-1,1:nz-1)));
disp('calculating_eps_r_y');
%%eps_r_y(i,j,k)is average of four cells(i,j,k)(i-1,j,k)(i,j,k-1)(i-1,j,k-1)
eps_r_y(2:nx,1:ny,2:nz)=0.25*(t_eps_r(material_3d_space(2:nx,1:ny,2:nz))+...
    t_eps_r(material_3d_space(1:nx-1,1:ny,2:nz))+...
    t_eps_r(material_3d_space(2:nx,1:ny,1:nz-1))+...
    t_eps_r(material_3d_space(1:nx-1,1:ny,1:nz-1)));
disp('calculating_eps_r_z');
%eps_r_z(i,j,k)is average of four cells(i,j,k)(i-1,j,k)(i,j-1,k)(i-1,j-1,k)
eps_r_z(2:nx,2:ny,1:nz)=0.25*(t_eps_r(material_3d_space(2:nx,1:ny-1,1:nz))+...
    t_eps_r(material_3d_space(1:nx-1,2:ny,1:nz))+...
    t_eps_r(material_3d_space(2:nx,2:ny,1:nz))+...
    t_eps_r(material_3d_space(1:nx-1,1:ny-1,1:nz)));
disp('calculating_sigma_e_x');
%%sigma_e_x(i,j,k)is average of four cells(i,j,k)(i,j-1,k)(i,j,k-1)(i,j-1,k-1)
sigma_e_x(1:nx,2:ny,2:nz)=0.25*(t_sigma_e(material_3d_space(1:nx,2:ny,2:nz))+...
    t_sigma_e(material_3d_space(1:nx,1:ny-1,2:nz))+...
    t_sigma_e(material_3d_space(1:nx,2:ny,1:nz-1))+...
    t_sigma_e(material_3d_space(1:nx,1:ny-1,1:nz-1)));
disp('calculating_sigma_e_y');
%%sigma_e_y(i,j,k)is average of four cells(i,j,k)(i-1,j,k)(i,j,k-1)(i-1,j,k-1)
sigma_e_y(2:nx,1:ny,2:nz)=0.25*(t_sigma_e(material_3d_space(2:nx,1:ny,2:nz))+...
    t_sigma_e(material_3d_space(1:nx-1,1:ny,2:nz))+...
    t_sigma_e(material_3d_space(2:nx,1:ny,1:nz-1))+...
    t_sigma_e(material_3d_space(1:nx-1,1:ny,1:nz-1)));
disp('calculating_sigma_e_z');
%sigma_e_z(i,j,k)is average of four cells(i,j,k)(i-1,j,k)(i,j-1,k)(i-1,j-1,k)
sigma_e_z(2:nx,2:ny,1:nz)=0.25*(t_sigma_e(material_3d_space(2:nx,1:ny-1,1:nz))+...
    t_sigma_e(material_3d_space(1:nx-1,2:ny,1:nz))+...
    t_sigma_e(material_3d_space(2:nx,2:ny,1:nz))+...
    t_sigma_e(material_3d_space(1:nx-1,1:ny-1,1:nz)));
disp('calculating_sigma_m_x');
%%sigma_m_x(i,j,k)is average of two cells(i,j,k)(i-1,j,k)
sigma_m_x(2:nx,1:ny,1:nz)=0.5*(t_sigma_m(material_3d_space(2:nx,1:ny,1:nz))+...
    t_sigma_m(material_3d_space(1:nx-1,1:ny,1:nz)));
disp('calculating_sigma_m_y');
%sigma_e_y(i,j,k)is average of two cells(i,j,k)(i,j-1,k)
sigma_m_y(1:nx,2:ny,1:nz)=0.5*(t_sigma_m(material_3d_space(1:nx,2:ny,1:nz))+...
    t_sigma_m(material_3d_space(1:nx,1:ny-1,1:nz)));
disp('calculating_sigma_m_z');
%sigma_e_z(i,j,k)is average of two cells(i,j,k)(i,j,k-1)
sigma_m_z(1:nx,1:ny,2:nz)=0.5*(t_sigma_m(material_3d_space(1:nx,1:ny,2:nz))+...
    t_sigma_m(material_3d_space(1:nx,1:ny,1:nz-1)));
disp('calculating_mu_r_x');
%mu_r_x(i,j,k)is average of two cells(i,j,k)(i-1,j,k)
mu_r_x(2:nx,1:ny,1:nz)=0.5*(t_mu_r(material_3d_space(2:nx,1:ny,1:nz))+...
    t_mu_r(material_3d_space(1:nx-1,1:ny,1:nz)));
disp('calculating_mu_r_y');
%mu_r_y(i,j,k)is average of two cells(i,j,k)(i,j-1,k)
mu_r_y(1:nx,2:ny,1:nz)=0.5*(t_mu_r(material_3d_space(1:nx,2:ny,1:nz))+...
    t_mu_r(material_3d_space(1:nx,1:ny-1,1:nz)));
disp('calculating_mu_r_z');
%mu_r_z(i,j,k)is average of two cells(i,j,k)(i,j,k-1)
mu_r_z(1:nx,1:ny,2:nz)=0.5*(t_mu_r(material_3d_space(1:nx,1:ny,2:nz))+...
    t_mu_r(material_3d_space(1:nx,1:ny,1:nz-1)));

%% creating_field_array
%create and initialize field and current arrays
Hx=zeros(nxp1,ny,nz);
Hy=zeros(nx,nyp1,nz);
Hz=zeros(nx,ny,nzp1);
Ex=zeros(nx,nyp1,nzp1);
Ey=zeros(nxp1,ny,nzp1);
Ez=zeros(nxp1,nyp1,nz);
%% define_output_parameters
% define sampled voltage
sampled_voltages(1).min_x=-(a/2);
sampled_voltages(1).min_y=-(b/2);
sampled_voltages(1).min_z=-L/4;
sampled_voltages(1).max_x=(a/2);
sampled_voltages(1).max_y=(b/2);
sampled_voltages(1).max_z=-L/4;
sampled_voltages(1).direction='yp';

sampled_voltages(2).min_x=-(a/2);
sampled_voltages(2).min_y=-(b/2);
sampled_voltages(2).min_z=L+(2e-3);
sampled_voltages(2).max_x=(a/2);
sampled_voltages(2).max_y=(b/2);
sampled_voltages(2).max_z=L+(2e-3);
sampled_voltages(2).direction='yp';

sampled_currents(1).min_x=-(a/2);
sampled_currents(1).min_y=-(b/2);
sampled_currents(1).min_z=-L/4;
sampled_currents(1).max_x=(a/2);
sampled_currents(1).max_y=(b/2);
sampled_currents(1).max_z=-L/4;
sampled_currents(1).direction='yp';

sampled_currents(2).min_x=-(a/2);
sampled_currents(2).min_y=-(b/2);
sampled_currents(2).min_z=L+(2e-3);
sampled_currents(2).max_x=(a/2);
sampled_currents(2).max_y=(b/2);
sampled_currents(2).max_z=L+(2e-3);
sampled_currents(2).direction='yp';

%define ports
ports(1).sampled_voltage_index=1;
ports(1).sampled_current_index=1;
ports(1).impedance=50;
ports(1).is_source_port=true;

ports(2).sampled_voltage_index=2;
ports(2).sampled_current_index=2;
ports(2).impedance=50;
ports(2).is_source_port=false;
%% initializing_sources_and _lumped_elements
%initialize sinusoidal waveform
time=dt*[0:number_of_time_steps-1].';
for ind=1:size(waveforms.cosine_modulated_gaussian,2)
frequency=...
waveforms.cosine_modulated_gaussian(ind).modulation_frequency;
tau=0.966/waveforms.cosine_modulated_gaussian(ind).bandwith;
waveforms.cosine_modulated_gaussian(ind).tau=tau;
t_0=4.5*waveforms.cosine_modulated_gaussian(ind).tau;
waveforms.cosine_modulated_gaussian(ind).t_0=t_0;
waveforms.cosine_modulated_gaussian(ind).waveform=...
cos(2*pi*frequency*(time-t_0)).*exp(-((time-t_0)/tau).^2);
end
%% initialize_updating_coefficients
% Coeffiecients updating Ex 
Cexe=(2*eps_r_x*eps_0-dt*sigma_e_x)./(2*eps_r_x*eps_0+dt*sigma_e_x);
Cexhz=(2*dt/dy)./(2*eps_r_x*eps_0+dt*sigma_e_x);
Cexhy=-(2*dt/dz)./(2*eps_r_x*eps_0+dt*sigma_e_x);
% Coeffiecients updating Ey 
Ceye=(2*eps_r_y*eps_0-dt*sigma_e_y)./(2*eps_r_y*eps_0+dt*sigma_e_y);
Ceyhx=(2*dt/dz)./(2*eps_r_y*eps_0+dt*sigma_e_y);
Ceyhz=-(2*dt/dx)./(2*eps_r_y*eps_0+dt*sigma_e_y);
% Coeffiecients updating Ez
Ceze =(2*eps_r_z*eps_0-dt*sigma_e_z)./(2*eps_r_z.*eps_0+dt*sigma_e_z);
Cezhy=(2*dt/dx)./(2*eps_r_z*eps_0+dt*sigma_e_z);
Cezhx=-(2*dt/dy)./(2*eps_r_z.*eps_0+dt*sigma_e_z);
%general magnetic field updating coefficients
% Coeffiecients updating Hx
Chxh =(2*mu_r_x*mu_0-dt*sigma_m_x)./(2*mu_r_x*mu_0+dt*sigma_m_x);
Chxez=-(2*dt/dy)./(2*mu_r_x*mu_0+dt*sigma_m_x);
Chxey=(2*dt/dz)./(2*mu_r_x*mu_0+dt*sigma_m_x);
% Coeffiecients updating Hy
Chyh=(2*mu_r_y*mu_0-dt*sigma_m_y)./(2*mu_r_y.*mu_0+dt*sigma_m_y);
Chyex=(-2*dt/dz)./(2*mu_r_y*mu_0+dt*sigma_m_y);
Chyez=(2*dt/dx)./(2*mu_r_y.*mu_0+dt*sigma_m_y);
% Coeffiecients updating Hz
Chzh =(2*mu_r_z*mu_0-dt*sigma_m_z)./(2*mu_r_z*mu_0+dt*sigma_m_z);
Chzey=-(2*dt/dx)./(2*mu_r_z*mu_0+dt*sigma_m_z);
Chzex=(2*dt/dy)./(2*mu_r_z*mu_0+dt*sigma_m_z);

%% initialize_boundary_conditions_3d
% define logical parameters for the conditions that will be used often
n_cpml_xn = abs( boundary.cpml_number_of_cells_xn);
n_cpml_xp = abs( boundary.cpml_number_of_cells_xp);
n_cpml_yn = abs( boundary.cpml_number_of_cells_yn);
n_cpml_yp = abs( boundary.cpml_number_of_cells_yp);
n_cpml_zn= abs ( boundary.cpml_number_of_cells_zn); 
n_cpml_zp = abs( boundary.cpml_number_of_cells_zp);
% Call CPML initialization routine if any side is CPML
% Initialize CPML boundary condition
pml_order = boundary.cpml_order;% order of the polynomial distribution 
sigma_max = boundary.cpml_sigma_max;
kappa_max = boundary.cpml_kappa_max;
alpha_order = boundary.cpml_alpha_order;
alpha_max = boundary.cpml_alpha_max;
eps_R = boundary.cpml_eps_R;
% Initialize cpml for xn region
sigma_opt=sigma_max*(n_cpml_xn+1)/(sqrt(eps_R)*150*pi*dx);
rho_e=((n_cpml_xn:-1:1)-0.75)/n_cpml_xn;
rho_m=((n_cpml_xn:-1:1)-0.25)/n_cpml_xn;
sigma_ex_xn=sigma_opt*abs(rho_e).^pml_order;
sigma_mx_xn=sigma_opt*abs(rho_m).^pml_order;
kappa_ex_xn=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_mx_xn=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ex_xn=alpha_max*abs(rho_e).^pml_order;
alpha_mx_xn=alpha_max*abs(rho_m).^pml_order;

cpml_b_ex_xn=exp((-dt/eps_0)...
    *((sigma_ex_xn./kappa_ex_xn)+alpha_ex_xn));
cpml_a_ex_xn=(1/dx)*(cpml_b_ex_xn-1).*sigma_ex_xn ...
    ./(kappa_ex_xn.*(sigma_ex_xn+kappa_ex_xn.*alpha_ex_xn));
cpml_b_mx_xn=exp((-dt/eps_0)...
    *((sigma_mx_xn./kappa_mx_xn)+alpha_mx_xn));
cpml_a_mx_xn=(1/dx)*(cpml_b_mx_xn-1).*sigma_mx_xn ...
    ./(kappa_mx_xn.*(sigma_mx_xn+kappa_mx_xn.*alpha_mx_xn));

Psi_eyx_xn = zeros (n_cpml_xn , ny, nzp1);
Psi_ezx_xn = zeros (n_cpml_xn , nyp1, nz);
Psi_hyx_xn = zeros (n_cpml_xn , nyp1, nz);
Psi_hzx_xn = zeros (n_cpml_xn , ny, nzp1 );

CPsi_eyx_xn = Ceyhz (2:n_cpml_xn +1,:,:)* dx;
CPsi_ezx_xn = Cezhy (2:n_cpml_xn +1,:,:)* dx;
CPsi_hyx_xn = Chyez (1:n_cpml_xn , : , :)*dx;
CPsi_hzx_xn = Chzey (1:n_cpml_xn,:,:) * dx ;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for i = 1: n_cpml_xn
Ceyhz (i+1,:,:) = Ceyhz (i +1,:,:)/ kappa_ex_xn (i); 
Cezhy (i+1,:,:) = Cezhy(i+1,:,:)/kappa_ex_xn (i);
Chyez(i,:,:) = Chyez(i,:,:)/ kappa_mx_xn(i);
Chzey(i,:,:) = Chzey(i, :,:)/ kappa_mx_xn(i);
end
% Initialize cpml for xp region
sigma_opt=sigma_max*(n_cpml_xp+1)/(sqrt(eps_R)*150*pi*dx);
rho_e=((1:1:n_cpml_xp)-0.75)/n_cpml_xp;
rho_m=((1:1:n_cpml_xp)-0.25)/n_cpml_xp;
sigma_ex_xp=sigma_opt*abs(rho_e).^pml_order;
sigma_mx_xp=sigma_opt*abs(rho_m).^pml_order;
kappa_ex_xp=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_mx_xp=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ex_xp=alpha_max*abs(rho_e).^pml_order;
alpha_mx_xp=alpha_max*abs(rho_m).^pml_order;

cpml_b_ex_xp=exp((-dt/eps_0)...
    *((sigma_ex_xp./kappa_ex_xp)+alpha_ex_xp));
cpml_a_ex_xp=(1/dx)*(cpml_b_ex_xp-1).*sigma_ex_xp ...
    ./(kappa_ex_xp.*(sigma_ex_xp+kappa_ex_xp.*alpha_ex_xp));
cpml_b_mx_xp=exp((-dt/eps_0)...
    *((sigma_mx_xp./kappa_mx_xp)+alpha_mx_xp));
cpml_a_mx_xp=1/dx*(cpml_b_mx_xp-1).*sigma_mx_xp ...
    ./(kappa_mx_xp.*(sigma_mx_xp+kappa_mx_xp.*alpha_mx_xp));

Psi_eyx_xp = zeros (n_cpml_xp , ny, nzp1);
Psi_ezx_xp = zeros (n_cpml_xp , nyp1, nz);
Psi_hyx_xp = zeros (n_cpml_xp , nyp1, nz);
Psi_hzx_xp = zeros (n_cpml_xp , ny, nzp1 );

CPsi_eyx_xp = Ceyhz (nxp1-n_cpml_xp:nx,:,:)* dx;
CPsi_ezx_xp = Cezhy (nxp1-n_cpml_xp:nx,:,:)* dx;
CPsi_hyx_xp = Chyez (nxp1-n_cpml_xp:nx,:,:)* dx;
CPsi_hzx_xp = Chzey (nxp1-n_cpml_xp:nx,:,:)* dx ;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for i = 1:n_cpml_xp
Ceyhz (nx-n_cpml_xp+i,:,:) = Ceyhz(nx-n_cpml_xp+i,:,:)/ kappa_ex_xp (i);
Cezhy (nx-n_cpml_xp+i,:,:) = Cezhy(nx-n_cpml_xp+i,:,:)/ kappa_ex_xp (i);
Chyez (nx-n_cpml_xp+i,:,:) = Chyez(nx-n_cpml_xp+i,:,:)/ kappa_mx_xp (i);
Chzey (nx-n_cpml_xp+i,:,:) = Chzey(nx-n_cpml_xp+i,:,:)/ kappa_mx_xp (i);
end

% Initialize cpml for yn region
sigma_opt=sigma_max*(n_cpml_yn+1)/(sqrt(eps_R)*150*pi*dy);
rho_e=((n_cpml_yn:-1:1)-0.75)/n_cpml_yn;
rho_m=((n_cpml_yn:-1:1)-0.25)/n_cpml_yn;
sigma_ey_yn=sigma_opt*abs(rho_e).^pml_order;
sigma_my_yn=sigma_opt*abs(rho_m).^pml_order;
kappa_ey_yn=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_my_yn=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ey_yn=alpha_max*abs(rho_e).^pml_order;
alpha_my_yn=alpha_max*abs(rho_m).^pml_order;

cpml_b_ey_yn=exp((-dt/eps_0)...
    *((sigma_ey_yn./kappa_ey_yn)+alpha_ey_yn));
cpml_a_ey_yn=1/dy*(cpml_b_ey_yn-1).*sigma_ey_yn ...
    ./(kappa_ey_yn.*(sigma_ey_yn+kappa_ey_yn.*alpha_ey_yn));
cpml_b_my_yn=exp((-dt/eps_0)...
    *((sigma_my_yn./kappa_my_yn)+alpha_my_yn));
cpml_a_my_yn=1/dy*(cpml_b_my_yn-1).*sigma_my_yn ...
    ./(kappa_my_yn.*(sigma_my_yn+kappa_my_yn.*alpha_my_yn));

Psi_exy_yn = zeros (nx,n_cpml_yn,nzp1); 
Psi_ezy_yn = zeros (nxp1,n_cpml_yn,nz);
Psi_hxy_yn = zeros (nxp1,n_cpml_yn,nz);
Psi_hzy_yn = zeros(nx,n_cpml_yn, nzp1 );
% Create and initialize 2D cpml convolution coefficients 
% Notice that Ey (1,:,:) and Ez (1,:,:) are not updated by cmpl 
CPsi_exy_yn = Cexhz (:,2:n_cpml_yn+1,:)*dy;
CPsi_ezy_yn = Cezhx (:,2:n_cpml_yn+1,:)*dy;
CPsi_hxy_yn = Chxez (:,1:n_cpml_yn  ,:)*dy;
CPsi_hzy_yn = Chzex (:,1:n_cpml_yn  ,:)*dy;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for j = 1: n_cpml_yn
Cexhz (:,j+1,:) = Cexhz (:,j+1,:)/ kappa_ey_yn (j);
Cezhx (:,j+1,:) = Cezhx (:,j+1,:)/ kappa_ey_yn (j);
Chxez (:,j  ,:) = Chxez (:,j  ,:)/ kappa_my_yn (j);
Chzex (:,j  ,:) = Chzex (:,j  ,:)/ kappa_my_yn (j);
end
% Initialize cpml for yp region
sigma_opt=sigma_max*(n_cpml_yp+1)/(sqrt(eps_R)*150*pi*dy);
rho_e=((1:1:n_cpml_yp)-0.75)/n_cpml_yp;
rho_m=((1:1:n_cpml_yp)-0.25)/n_cpml_yp;
sigma_ey_yp=sigma_opt*abs(rho_e).^pml_order;
sigma_my_yp=sigma_opt*abs(rho_m).^pml_order;
kappa_ey_yp=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_my_yp=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ey_yp=alpha_max*abs(rho_e).^pml_order;
alpha_my_yp=alpha_max*abs(rho_m).^pml_order;

cpml_b_ey_yp=exp((-dt/eps_0)...
    *((sigma_ey_yp./kappa_ey_yp)+alpha_ey_yp));
cpml_a_ey_yp=1/dx*(cpml_b_ey_yp-1).*sigma_ey_yp ...
    ./(kappa_ey_yp.*(sigma_ey_yp+kappa_ey_yp.*alpha_ey_yp));
cpml_b_my_yp=exp((-dt/eps_0)...
    *((sigma_my_yp./kappa_my_yp)+alpha_my_yp));
cpml_a_my_yp=1/dy*(cpml_b_my_yp-1).*sigma_my_yp ...
    ./(kappa_my_yp.*(sigma_my_yp+kappa_my_yp.*alpha_my_yp));

Psi_exy_yp = zeros (nx,n_cpml_yp,nzp1); 
Psi_ezy_yp = zeros (nxp1,n_cpml_yp,nz);
Psi_hxy_yp = zeros (nxp1,n_cpml_yp,nz);
Psi_hzy_yp = zeros (nx,n_cpml_yp,nzp1);
% Create and initialize 2D cpml convolution coefficients
% Notice that Ey (nxp1,:,:) and Ez (nxp1 ,:,:) are not updated by cmpl 
CPsi_exy_yp = Cexhz (:,nyp1-n_cpml_yp:ny,:)* dy;
CPsi_ezy_yp = Cezhx (:,nyp1-n_cpml_yp:ny,:)* dy;
CPsi_hxy_yp = Chxez (:,nyp1-n_cpml_yp:ny,:)* dy;
CPsi_hzy_yp = Chzex (:,nyp1-n_cpml_yp:ny,:)* dy;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (nxp1,:,:) and Ez (nxp1 ,:,:) are not updated by cmpl 
for j = 1:n_cpml_yp
Cexhz (:,ny-n_cpml_yp+j,:) = Cexhz (:,ny-n_cpml_yp+j,:)/ kappa_ey_yp (j);
Cezhx (:,ny-n_cpml_yp+j,:) = Cezhx (:,ny-n_cpml_yp+j,:)/ kappa_ey_yp (j);
Chxez (:,ny-n_cpml_yp+j,:) = Chxez (:,ny-n_cpml_yp+j,:)/ kappa_my_yp (j);
Chzex (:,ny-n_cpml_yp+j,:) = Chzex (:,ny-n_cpml_yp+j,:)/ kappa_my_yp (j);
end

% Initialize cpml for zn region
sigma_opt=sigma_max*(n_cpml_zn+1)/(sqrt(eps_R)*150*pi*dz);
rho_e=((n_cpml_zn:-1:1)-0.75)/n_cpml_zn;
rho_m=((n_cpml_zn:-1:1)-0.25)/n_cpml_zn;
sigma_ez_zn=sigma_opt*abs(rho_e).^pml_order;
sigma_mz_zn=sigma_opt*abs(rho_m).^pml_order;
kappa_ez_zn=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_mz_zn=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ez_zn=alpha_max*abs(rho_e).^pml_order;
alpha_mz_zn=alpha_max*abs(rho_m).^pml_order;

cpml_b_ez_zn=exp((-dt/eps_0)...
    *((sigma_ez_zn./kappa_ez_zn)+alpha_ez_zn));
cpml_a_ez_zn=1/dz*(cpml_b_ez_zn-1).*sigma_ez_zn ...
    ./(kappa_ez_zn.*(sigma_ez_zn+kappa_ez_zn.*alpha_ez_zn));
cpml_b_mz_zn=exp((-dt/eps_0)...
    *((sigma_mz_zn./kappa_mz_zn)+alpha_mz_zn));
cpml_a_mz_zn=1/dz*(cpml_b_mz_zn-1).*sigma_mz_zn ...
    ./(kappa_mz_zn.*(sigma_mz_zn+kappa_mz_zn.*alpha_mz_zn));

Psi_eyz_zn = zeros (nxp1,ny,n_cpml_zn);
Psi_exz_zn = zeros (nx,nyp1,n_cpml_zn);
Psi_hyz_zn = zeros (nx,nyp1,n_cpml_zn);
Psi_hxz_zn = zeros (nxp1,ny,n_cpml_zn );
% Create and initialize 2D cpml convolution coefficients 
% Notice that Ey (1,:,:) and Ez (1,:,:) are not updated by cmpl 
CPsi_eyz_zn = Ceyhx (:,:,2: n_cpml_zn+1)* dz;
CPsi_exz_zn = Cexhy (:,:,2: n_cpml_zn+1)* dz;
CPsi_hyz_zn = Chyex (:,:,1: n_cpml_zn)* dz;
CPsi_hxz_zn = Chxey (:,:,1: n_cpml_zn)* dz ;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for j = 1: n_cpml_zn
Cexhy (:,:,j+1) = Cexhy (:,:,j+1)/ kappa_ez_zn (j); 
Ceyhx (:,:,j+1) = Ceyhx(:,:,j+1)/kappa_ez_zn (j);
Chxey (:,:,j) = Chxey(:,:,j)/ kappa_mz_zn(j);
Chyex (:,:,j) = Chyex(:,:,j)/ kappa_mz_zn(j);
end
% Initialize cpml for zp region
sigma_opt=sigma_max*(n_cpml_zp+1)/(sqrt(eps_R)*150*pi*dz);
rho_e=((1:1:n_cpml_zp)-0.75)/n_cpml_zp;
rho_m=((1:1:n_cpml_zp)-0.25)/n_cpml_zp;
sigma_ez_zp=sigma_opt*abs(rho_e).^pml_order;
sigma_mz_zp=sigma_opt*abs(rho_m).^pml_order;
kappa_ez_zp=1+(kappa_max-1)*abs(rho_e).^pml_order;
kappa_mz_zp=1+(kappa_max-1)*abs(rho_m).^pml_order;
alpha_ez_zp=alpha_max*abs(rho_e).^pml_order;
alpha_mz_zp=alpha_max*abs(rho_m).^pml_order;

cpml_b_ez_zp=exp((-dt/eps_0)...
    *((sigma_ez_zp./kappa_ez_zp)+alpha_ez_zp));
cpml_a_ez_zp=1/dz*(cpml_b_ez_zp-1).*sigma_ez_zp ...
    ./(kappa_ez_zp.*(sigma_ez_zp+kappa_ez_zp.*alpha_ez_zp));
cpml_b_mz_zp=exp((-dt/eps_0)...
    *((sigma_mz_zp./kappa_mz_zp)+alpha_mz_zp));
cpml_a_mz_zp=1/dz*(cpml_b_mz_zp-1).*sigma_mz_zp ...
    ./(kappa_mz_zp.*(sigma_mz_zp+kappa_mz_zp.*alpha_mz_zp));

Psi_eyz_zp = zeros (nxp1,ny,n_cpml_zp);
Psi_exz_zp = zeros (nx,nyp1,n_cpml_zp);
Psi_hyz_zp = zeros (nx,nyp1,n_cpml_zp);
Psi_hxz_zp = zeros (nxp1,ny,n_cpml_zp);
% Create and initialize 2D cpml convolution coefficients 
% Notice that Ey (1,:,:) and Ez (1,:,:) are not updated by cmpl 
CPsi_eyz_zp = Ceyhx (:,:,nzp1-n_cpml_zp:nz)* dz;
CPsi_exz_zp = Cexhy (:,:,nzp1-n_cpml_zp:nz)* dz;
CPsi_hyz_zp = Chyex (:,:,nzp1-n_cpml_zp:nz)* dz;
CPsi_hxz_zp = Chxey (:,:,nzp1-n_cpml_zp:nz)* dz;
% Adjust FDTD coefficients in the CPML region 
% Notice that Ey (1 ,:,:) and Ez (1,:,:) are not updated by cmpl
for i = 1: n_cpml_zp
Cexhy(:,:,nz-n_cpml_zp+i) = Cexhy(:,:,nz-n_cpml_zp+i)/kappa_ez_zp(i); 
Ceyhx(:,:,nz-n_cpml_zp+i) = Ceyhx(:,:,nz-n_cpml_zp+i)/kappa_ez_zp(i);
Chxey(:,:,nz-n_cpml_zp+i) = Chxey(:,:,nz-n_cpml_zp+i)/kappa_mz_zp(i);
Chyex(:,:,nz-n_cpml_zp+i) = Chyex(:,:,nz-n_cpml_zp+i)/kappa_mz_zp(i);
end
%% initialize_output_parameters_3d
number_of_sampled_voltages=size(sampled_voltages,2);
number_of_sampled_currents=size(sampled_currents,2);
number_of_ports=size(ports,2);
%initial frequency domain parameters
frequency_domain.frequencies=[frequency_domain.start:...
    frequency_domain.step:frequency_domain.end];
frequency_domain.number_of_frequencies=...
    size(frequency_domain.frequencies,2);
%initialize sampled voltage terms
for ind=1:number_of_sampled_voltages
    is=round((sampled_voltages(ind).min_x-fdtd_domain.min_x)/dx);
    js=round((sampled_voltages(ind).min_y-fdtd_domain.min_y)/dy)+1;
    ks=round((sampled_voltages(ind).min_z-fdtd_domain.min_z)/dz)+1;
    ie=round((sampled_voltages(ind).max_x-fdtd_domain.min_x)/dx)+1;
    je=round((sampled_voltages(ind).max_y-fdtd_domain.min_y)/dy)+1;
    ke=round((sampled_voltages(ind).max_z-fdtd_domain.min_z)/dz)+1;
    sampled_voltages(ind).is=is;
    sampled_voltages(ind).js=js;
    sampled_voltages(ind).ks=ks;
    sampled_voltages(ind).ie=ie;
    sampled_voltages(ind).je=je;
    sampled_voltages(ind).ke=ke;
    sampled_voltages(ind).samples_value=zeros(1,number_of_time_steps);

    switch (sampled_voltages(ind).direction(1))
        case 'x'
            fi=create_linear_index_list(Ex,is:ie-1,js:je-1,ks:ke);
            sampled_voltages(ind).Csvf=-dx/((je-js+1)*(ke-ks+1));
        case 'y'
            fi=create_linear_index_list(Ey,is:ie,js:je-1,ks:ke);
            sampled_voltages(ind).Csvf=-dy/((ke-ks+1)*(ie-is+1));            
    end
    if strcmp(sampled_voltages(ind).direction(2),'n')
        sampled_voltages(ind).Csvf=-sampled_voltages(ind).Csvf;
    end
    sampled_voltages(ind).field_indices=fi;
    sampled_voltages(ind).time=([1:number_of_time_steps])*dt;
end

%initialize sampled current terms
for ind=1:number_of_sampled_currents
    is=round((sampled_currents(ind).min_x-fdtd_domain.min_x)/dx)+1;
    js=round((sampled_currents(ind).min_y-fdtd_domain.min_y)/dy)+1;
    ks=round((sampled_currents(ind).min_z-fdtd_domain.min_z)/dz)+1;
    ie=round((sampled_currents(ind).max_x-fdtd_domain.min_x)/dx)+1;
    je=round((sampled_currents(ind).max_y-fdtd_domain.min_y)/dy)+1;
    ke=round((sampled_currents(ind).max_z-fdtd_domain.min_z)/dz)+1;
    sampled_currents(ind).is=is;
    sampled_currents(ind).js=js;
    sampled_currents(ind).ks=ks;
    sampled_currents(ind).ie=ie;
    sampled_currents(ind).je=je;
    sampled_currents(ind).ke=ke;
    sampled_currents(ind).samples_value=zeros(1,number_of_time_steps);
    sampled_currents(ind).time=([1:number_of_time_steps]-0.5)*dt;
 end    

Hx_sample=zeros(nxp1,ny,nz);
Hy_sample=zeros(nx-21,nyp1-21,nz-81);
Hz_sample=zeros(nx-21,ny-21,nzp1-81);
Ex_sample=zeros(nx-21,nyp1-21,nzp1-81);
Ey_sample=zeros(nxp1,ny,nzp1);
Ez_sample=zeros(nxp1-21,nyp1-21,nz-81);

figure(1)
for ind=1:number_of_brick
hold on
patch([brick(ind).min_y brick(ind).min_y brick(ind).max_y ...
    brick(ind).max_y],[brick(ind).min_x brick(ind).max_x...
    brick(ind).max_x brick(ind).min_x],[brick(ind).min_z...
    brick(ind).min_z brick(ind).min_z brick(ind).min_z],...
    material_type(+1).color)
patch([brick(ind).min_y brick(ind).min_y brick(ind).max_y...
    brick(ind).max_y], [brick(ind).min_x brick(ind).max_x...
    brick(ind).max_x brick(ind).min_x], [brick(ind).max_z...
    brick(ind).max_z brick(ind).max_z brick(ind).max_z],...
    material_type(+1).color)
patch([brick(ind).min_y brick(ind).min_y brick(ind).min_y...
    brick(ind).min_y], [brick(ind).min_x brick(ind).max_x...
    brick(ind).max_x brick(ind).min_x], [brick(ind).min_z...
    brick(ind).min_z brick(ind).max_z brick(ind).max_z],...
    material_type(+1).color)
patch([brick(ind).max_y brick(ind).max_y brick(ind).max_y...
    brick(ind).max_y], [brick(ind).min_x brick(ind).max_x...
    brick(ind).max_x brick(ind).min_x], [brick(ind).min_z...
    brick(ind).min_z brick(ind).max_z brick(ind).max_z],...
    material_type(+1).color)
patch([brick(ind).min_y brick(ind).min_y brick(ind).max_y...
    brick(ind).max_y], [brick(ind).min_x brick(ind).min_x...
    brick(ind).min_x brick(ind).min_x], [brick(ind).min_z...
    brick(ind).max_z brick(ind).max_z brick(ind).min_z],...
    material_type(+1).color)
patch([brick(ind).min_y brick(ind).min_y brick(ind).max_y...
    brick(ind).max_y], [brick(ind).max_x brick(ind).max_x...
    brick(ind).max_x brick(ind).max_x], [brick(ind).min_z...
    brick(ind).max_z brick(ind).max_z brick(ind).min_z],...
    material_type(+1).color)
alpha(.009)
xlim([-0.09 0.09])
ylim([-0.09 0.09])
zlim([-0.2 0.6])
view(3)
end
drawnow

cpml_b_mx_xn=cpml_b_mx_xn';
cpml_a_mx_xn=cpml_a_mx_xn';
cpml_b_mx_xp=cpml_b_mx_xp';
cpml_a_mx_xp=cpml_a_mx_xp';
cpml_b_ex_xn=cpml_b_ex_xn';
cpml_a_ex_xn=cpml_a_ex_xn';
cpml_b_ex_xp=cpml_b_ex_xp';
cpml_a_ex_xp=cpml_a_ex_xp';

cpml_b_mz_zn=permute(cpml_b_mz_zn,[1 3 2]);
cpml_a_mz_zn=permute(cpml_a_mz_zn,[1 3 2]);
cpml_b_mz_zp=permute(cpml_b_mz_zp,[1 3 2]);
cpml_a_mz_zp=permute(cpml_a_mz_zp,[1 3 2]);

cpml_b_ez_zn=permute(cpml_b_ez_zn,[1 3 2]);
cpml_a_ez_zn=permute(cpml_a_ez_zn,[1 3 2]);
cpml_b_ez_zp=permute(cpml_b_ez_zp,[1 3 2]);
cpml_a_ez_zp=permute(cpml_a_ez_zp,[1 3 2]);

tic
%% run_fdtd_time_marching_loop
for time_step=1:number_of_time_steps
% update_magnetic_fields_CPML
Psi_hyx_xn(:,:,:) = cpml_b_mx_xn.*Psi_hyx_xn+cpml_a_mx_xn.*...
    ( Ez (2:n_cpml_xn+1,:,:) - Ez(1:n_cpml_xn,:,:));
Psi_hzx_xn(:,:,:) = cpml_b_mx_xn.* Psi_hzx_xn+cpml_a_mx_xn.*...
    ( Ey (2:n_cpml_xn+1,:,:) - Ey(1:n_cpml_xn ,:,:));
Psi_hyx_xp (:,:,:) = cpml_b_mx_xp.* Psi_hyx_xp+...
    cpml_a_mx_xp.*(Ez (nx-n_cpml_xp+2:nx+1,:,:)-Ez(nx-n_cpml_xp+1:nx,:,:));
Psi_hzx_xp (:,:,:) = cpml_b_mx_xp.* Psi_hzx_xp+...
    cpml_a_mx_xp.*(Ey (nx-n_cpml_xp+2:nx+1,:,:)-Ey(1+nx-n_cpml_xp:nx,:,:));

Psi_hxy_yn(:,:,:)=cpml_b_my_yn.*Psi_hxy_yn+cpml_a_my_yn.*...
    (Ez(:,2:n_cpml_yn+1,:) - Ez(:,1:n_cpml_yn,:));
Psi_hzy_yn(:,:,:)=cpml_b_my_yn.*Psi_hzy_yn+cpml_a_my_yn.*...
    (Ex(:,2:n_cpml_yn+1,:) - Ex(:,1:n_cpml_yn,:));
Psi_hxy_yp(:,:,:)=cpml_b_my_yp.*Psi_hxy_yp(:,:,:)+cpml_a_my_yp.*...
    (Ez(:,ny-n_cpml_yp+2:ny+1,:)-Ez(:,ny-n_cpml_yp+1:ny,:));
Psi_hzy_yp(:,:,:)=cpml_b_my_yp.*Psi_hzy_yp(:,:,:)+cpml_a_my_yp.*...
    (Ex(:,ny-n_cpml_yp+2:ny+1,:)-Ex(:,1+ny-n_cpml_yp:ny,:));

Psi_hxz_zn(:,:,:)=(cpml_b_mz_zn).*Psi_hxz_zn(:,:,:)+(cpml_a_mz_zn).*...
    (Ey(:,:,2:n_cpml_zn+1) - Ey(:,:,1:n_cpml_zn));
Psi_hyz_zn(:,:,:)=(cpml_b_mz_zn).*Psi_hyz_zn(:,:,:)+(cpml_a_mz_zn).*...
    (Ex(:,:,2:n_cpml_zn+1) - Ex(:,:,1:n_cpml_zn));
Psi_hyz_zp (:,:,:) = (cpml_b_mz_zp).*Psi_hyz_zp(:,:,:)+(cpml_a_mz_zp).*...
    (Ex(:,:,nz-n_cpml_zp+2:nz+1)-Ex(:,:,nz-n_cpml_zp+1:nz));
Psi_hxz_zp (:,:,:) = (cpml_b_mz_zp).*Psi_hxz_zp(:,:,:)+(cpml_a_mz_zp).*...
    (Ey(:,:,nz-n_cpml_zp+2:nz+1)-Ey(:,:,nz-n_cpml_zp+1:nz));

Hy(1:n_cpml_xn,:,:)=Hy(1:n_cpml_xn,:,:)+CPsi_hyx_xn(:,:,:).*Psi_hyx_xn(:,:,:);
Hz(1:n_cpml_xn,:,:)=Hz(1:n_cpml_xn,:,:)+CPsi_hzx_xn(:,:,:).*Psi_hzx_xn(:,:,:);

Hy(nx-n_cpml_xp+1:nx,:,:)=Hy(nx-n_cpml_xp+1:nx,:,:)+CPsi_hyx_xp(:,:,:).*Psi_hyx_xp(:,:,:);
Hz(nx-n_cpml_xp+1:nx,:,:)=Hz(nx-n_cpml_xp+1:nx,:,:)+CPsi_hzx_xp(:,:,:).*Psi_hzx_xp(:,:,:);

Hx(:,1:n_cpml_yn,:)=Hx(:,1:n_cpml_yn,:)+CPsi_hxy_yn(:,:,:).*Psi_hxy_yn(:,:,:);
Hz(:,1:n_cpml_yn,:)=Hz(:,1:n_cpml_yn,:)+CPsi_hzy_yn(:,:,:).*Psi_hzy_yn(:,:,:);

Hx(:,ny-n_cpml_yp+1:ny,:)=Hx(:,ny-n_cpml_yp+1:ny,:)+CPsi_hxy_yp(:,:,:).*Psi_hxy_yp(:,:,:);
Hz(:,ny-n_cpml_yp+1:ny,:)=Hz(:,ny-n_cpml_yp+1:ny,:)+CPsi_hzy_yp(:,:,:).*Psi_hzy_yp(:,:,:);

Hx(:,:,1:n_cpml_zn)=Hx(:,:,1:n_cpml_zn)+CPsi_hxz_zn(:,:,:).*Psi_hxz_zn(:,:,:);
Hy(:,:,1:n_cpml_zn)=Hy(:,:,1:n_cpml_zn)+CPsi_hyz_zn(:,:,:).*Psi_hyz_zn(:,:,:);

Hx(:,:,nz-n_cpml_zp+1:nz)=Hx(:,:,nz-n_cpml_zp+1:nz)+CPsi_hxz_zp(:,:,:).* Psi_hxz_zp(:,:,:);
Hy(:,:,nz-n_cpml_zp+1:nz)=Hy(:,:,nz-n_cpml_zp+1:nz)+CPsi_hyz_zp(:,:,:).* Psi_hyz_zp(:,:,:);
%% update_magnetic_fields
Hx = Chxh.* Hx+Chxey .*(Ey (1:nxp1,1:ny,2:nzp1)- Ey (1:nxp1,1:ny,1:nz))+...
    Chxez.*(Ez(1:nxp1,2:nyp1,1:nz)- Ez (1:nxp1,1:ny,1:nz));
Hy = Chyh.* Hy+Chyez .*(Ez (2:nxp1,1:nyp1,1:nz)- Ez (1:nx,1:nyp1,1:nz))+...
    Chyex.*(Ex(1:nx,1:nyp1,2:nzp1)- Ex (1:nx,1:nyp1,1:nz));
Hz = Chzh.* Hz+Chzex .*(Ex (1:nx,2:nyp1,1:nzp1)- Ex (1:nx,1:ny,1:nzp1))+...
    Chzey.*(Ey(2:nxp1,1:ny,1:nzp1)- Ey (1:nx,1:ny,1:nzp1));
%% update_electric_fields_for_CPML
Psi_eyx_xn(:,:,:)= cpml_b_ex_xn.* Psi_eyx_xn(:,:,:)+ cpml_a_ex_xn.*...
    (Hz(2:n_cpml_xn+1,:,:)-Hz(1:n_cpml_xn,:,:));
Psi_ezx_xn(:,:,:)= cpml_b_ex_xn.* Psi_ezx_xn(:,:,:)+ cpml_a_ex_xn.*...
    (Hy(2:n_cpml_xn+1,:,:)-Hy(1:n_cpml_xn,:,:));

Ey (2:n_cpml_xn+1,:,:) = Ey (2:n_cpml_xn+1,:,:)+CPsi_eyx_xn.*Psi_eyx_xn;
Ez (2:n_cpml_xn+1,:,:) = Ez (2:n_cpml_xn+1,:,:)+CPsi_ezx_xn.*Psi_ezx_xn; 

Psi_eyx_xp(:,:,:) = cpml_b_ex_xp.*Psi_eyx_xp(:,:,:)+cpml_a_ex_xp.*...
    (Hz(1+nx-n_cpml_xp:nx ,:,:) - Hz(nx-n_cpml_xp:nx-1,:,:));
Psi_ezx_xp(:,:,:) = cpml_b_ex_xp.*Psi_ezx_xp(:,:,:)+cpml_a_ex_xp.*...
    (Hy(1+nx-n_cpml_xp:nx,:,:) - Hy(nx-n_cpml_xp:nx-1,:,:));

Ey(nx-n_cpml_xp+1:nx,:,:)=Ey(nx-n_cpml_xp+1:nx,:,:)+CPsi_eyx_xp.*Psi_eyx_xp;
Ez(nx-n_cpml_xp+1:nx,:,:)=Ez(nx-n_cpml_xp+1:nx,:,:)+CPsi_ezx_xp.*Psi_ezx_xp;

Psi_exy_yn(:,:,:) = cpml_b_ey_yn.* Psi_exy_yn(:,:,:)+ cpml_a_ey_yn.*...
    (Hz(:,2:n_cpml_yn+1,:) - Hz(:,1:n_cpml_yn,:));
Psi_ezy_yn(:,:,:) = cpml_b_ey_yn.* Psi_ezy_yn(:,:,:)+ cpml_a_ey_yn.*...
    (Hx(:,2:n_cpml_yn+1,:) - Hx(:,1:n_cpml_yn,:));

Ex (:,2:n_cpml_yn+1,:)=Ex(:,2:n_cpml_yn +1,:)+ CPsi_exy_yn .* Psi_exy_yn;
Ez (:,2:n_cpml_yn+1,:)=Ez(:,2:n_cpml_yn +1,:)+ CPsi_ezy_yn .* Psi_ezy_yn; 

Psi_exy_yp(:,:,:) = cpml_b_ey_yp.*Psi_exy_yp(:,:,:)+cpml_a_ey_yp.*...
    (Hz(:,1+ny-n_cpml_yp:ny,:) - Hz(:,ny-n_cpml_yp:ny-1,:));
Psi_ezy_yp(:,:,:) = cpml_b_ey_yp.*Psi_ezy_yp(:,:,:)+cpml_a_ey_yp.*...
    (Hx(:,1+ny-n_cpml_yp:ny,:) - Hx(:,ny-n_cpml_yp:ny-1,:));

Ex(:,ny-n_cpml_yp+1:ny,:)=Ex(:,ny-n_cpml_yp+1:ny,:)+CPsi_exy_yp.*Psi_exy_yp;
Ez(:,ny-n_cpml_yp+1:ny,:)=Ez(:,ny-n_cpml_yp+1:ny,:)+CPsi_ezy_yp.*Psi_ezy_yp;

Psi_exz_zn(:,:,:) = cpml_b_ez_zn.* Psi_exz_zn(:,:,:)+ cpml_a_ez_zn.*...
    (Hy(:,:,2:n_cpml_zn+1) - Hy(:,:,1:n_cpml_zn));
Psi_eyz_zn(:,:,:) = cpml_b_ez_zn.* Psi_eyz_zn(:,:,:)+ cpml_a_ez_zn.*...
    (Hx(:,:,2:n_cpml_zn+1) - Hx(:,:,1:n_cpml_zn));

Ex (:,:,2:n_cpml_zn+1) = Ex(:,:,2:n_cpml_zn+1)+ CPsi_exz_zn .* Psi_exz_zn;
Ey (:,:,2:n_cpml_zn +1)= Ey(:,:,2:n_cpml_zn+1)+ CPsi_eyz_zn .* Psi_eyz_zn; 

Psi_exz_zp(:,:,:) = cpml_b_ez_zp .*Psi_exz_zp(:,:,:)+cpml_a_ez_zp.*...
    (Hy(:,:,1+nz-n_cpml_zp:nz) - Hy(:,:,nz-n_cpml_zp:nz-1));
Psi_eyz_zp(:,:,:) = cpml_b_ez_zp.*Psi_eyz_zp(:,:,:)+cpml_a_ez_zp.*...
    (Hx(:,:,1+nz-n_cpml_zp:nz) - Hx(:,:,nz-n_cpml_zp:nz-1));

Ex(:,:,nz-n_cpml_zp+1:nz)=Ex(:,:,nz-n_cpml_zp+1:nz)+CPsi_exz_zp.*Psi_exz_zp;
Ey(:,:,nz-n_cpml_zp+1:nz)=Ey(:,:,nz-n_cpml_zp+1:nz)+CPsi_eyz_zp.*Psi_eyz_zp;
%% update_electric_fields
Ex(1:nx,2:ny,2:nz)=Cexe(1:nx,2:ny,2:nz).*Ex(1:nx,2:ny,2:nz)+...
    Cexhz(1:nx,2:ny,2:nz).*(Hz(1:nx,2:ny,2:nz)-Hz(1:nx,1:ny-1,2:nz))+...
    Cexhy(1:nx,2:ny,2:nz).*(Hy(1:nx,2:ny,2:nz)-Hy(1:nx,2:ny,1:nz-1));
Ey(2:nx,1:ny,2:nz)=Ceye(2:nx,1:ny,2:nz).*Ey(2:nx,1:ny,2:nz)+...
    Ceyhx(2:nx,1:ny,2:nz).*(Hx(2:nx,1:ny,2:nz)-Hx(2:nx,1:ny,1:nz-1))+...
    Ceyhz(2:nx,1:ny,2:nz).*(Hz(2:nx,1:ny,2:nz)-Hz(1:nx-1,1:ny,2:nz));
Ez (2:nx,2:ny,1:nz)=Ceze(2:nx,2:ny,1:nz).*Ez(2:nx,2:ny,1:nz)+...
    Cezhy(2:nx,2:ny,1:nz).*(Hy(2:nx,2:ny,1:nz)-Hy(1:nx-1,2:ny,1:nz))+...
    Cezhx(2:nx,2:ny,1:nz).*(Hx(2:nx,2:ny,1:nz)-Hx(2:nx,1:ny-1,1:nz));
%% update_voltage_sources
Ey(:,:,n_cpml_zn+ceil(0.1*L/dz))=Ey(:,:,n_cpml_zn+ceil(0.1*L/dz))+...
                cos(2*pi*frequency*(time_step.*dt-t_0)).*exp(-((time_step.*dt-t_0)/tau).^2);%  sin(2*pi*waveforms.sinusoidal(1).frequency*time_step.*dt);%
%% capture_and_display_sampled_fields
if time_step>0
% capturing sampled currents
for ind=1:number_of_sampled_currents
    is=sampled_currents(ind).is;
    js=sampled_currents(ind).js;
    ks=sampled_currents(ind).ks;
    ie=sampled_currents(ind).ie;
    je=sampled_currents(ind).je;
    ke=sampled_currents(ind).ke;
%     switch (sampled_currents(ind).direction(1))
%         case 'x'
%             sampled_value=...
%                 +dy*sum(sum(sum(Hy(ie-1,js:je,ks-1))))...
%                 +dz*sum(sum(sum(Hz(ie-1,je,ks:ke))))...
%                 -dy*sum(sum(sum(Hy(ie-1,js:je,ke))))...
%                 -dz*sum(sum(sum(Hz(ie-1,js-1,ks:ke))));
%         case 'y'
            sampled_value=...
                +dz*sum(sum(sum(Hz(is-1,je-1,ks:ke))))...
                +dx*sum(sum(sum(Hx(is:ie,je-1,ke))))...
                -dz*sum(sum(sum(Hz(ie,je-1,ks:ke))))...
                -dx*sum(sum(sum(Hx(is:ie,je-1,ks-1))));            
%         case 'z'
%             sampled_value=...
%                 +dx*sum(sum(sum(Hx(ie:ie,js-1,ke-1))))...
%                 +dy*sum(sum(sum(Hy(ie,js:je,ke-1))))...
%                 -dx*sum(sum(sum(Hx(is:ie,je,ke-1))))...
%                 -dy*sum(sum(sum(Hy(is-1,js:je,ke-1))));   
%     end
%     if strcmp(sampled_currents(ind).direction(2),'n')
%         sampled_value=-sampled_value;
%     end
    sampled_currents(ind).sampled_value(time_step)=sampled_value;
end

%capturing sampled voltages
for ind=1:number_of_sampled_voltages
    fi=sampled_voltages(ind).field_indices;
    Csvf=sampled_voltages(ind).Csvf;
%     switch (sampled_voltages(ind).direction(1))
%         case 'x'
%             sampled_value=Csvf*sum(Ex(fi));
%         case 'y'
            sampled_value=Csvf*sum(Ey(fi));
%     end
    sampled_voltages(ind).sampled_value(time_step)=sampled_value;
end

if mod(time_step,ceil(50))==0||time_step==totalTimeStep%ceil(number_of_time_steps/2)
    
%calculate_frequency_domain_output
frequency_array=frequency_domain.frequencies;
%sampled currents in frequency domain
for ind=1:number_of_sampled_currents
    x=sampled_currents(ind).sampled_value;
    time_shift=-dt/2;
    [X]=time_to_frequency_domain(x,dt,frequency_array,time_shift,complex_number);
    sampled_currents(ind).frequency_domain_value=X;
    sampled_currents(ind).frequencies=frequency_array;
end
    
%sampled voltage in frequency domain
for ind=1:number_of_sampled_voltages
    x=sampled_voltages(ind).sampled_value;
    time_shift=0;
    [X]=time_to_frequency_domain(x,dt,frequency_array,time_shift,complex_number);
    sampled_voltages(ind).frequency_domain_value=X;
    sampled_voltages(ind).frequencies=frequency_array;
end    

%calculation of S_parameters
for ind=1:number_of_ports
    svi=ports(ind).sampled_voltage_index;
    sci=ports(ind).sampled_current_index;
    Z=ports(ind).impedance;
    V=sampled_voltages(svi).frequency_domain_value;
    I=sampled_currents(sci).frequency_domain_value;
    ports(ind).a=0.5*(V+Z.*I)./sqrt(real(Z));
    ports(ind).b=0.5*(V-conj(Z).*I)./sqrt(real(Z));
    ports(ind).frequencies=frequency_array;
end
for ind=1:number_of_ports
    if ports(ind).is_source_port==true
        for oind=1:number_of_ports
            ports(ind).S(oind).values=ports(oind).b./ports(ind).a;
        end
    end
end
%% figure for S_parameters
figure(2)
for ind=1:number_of_ports
    if ports(ind).is_source_port==true
        frequencies=ports(ind).frequencies*1e-9;
        for oind=number_of_ports
            S=ports(ind).S(oind).values;
            sdb=20*log10(abs(S));
            Sphase=angle(S)*180/pi;
         
            subplot(2,1,1);
            plot(frequencies, sdb, 'b-','linewidth',1.5);
            title(['S' num2str(oind) num2str(ind)],'fontsize',12);
            xlabel('frequency (GHz)','fontsize',12);
            ylabel('magnitude (dB)','fontsize',12);
%             legend('S 21')
            grid on
            subplot(2,1,2);
            plot(frequencies,(Sphase),'r-','linewidth',1.5);
            xlabel('frequency (GHz)','fontsize',12);
            ylabel('phase (degrees)','fontsize',12);
            grid on
            drawnow
        end
    end
end

Hx_sample=Hx(1:end,1:end,1:end);
Hy_sample=Hy(1:end,1:end,1:end);
Hz_sample=Hz(1:end,1:end,1:end);
Ex_sample=Ex(1:end,1:end,1:end);
Ey_sample=(Ey(1:end,1:end,1:end));
Ez_sample=Ez(1:end,1:end,1:end);

%% Ey plot
figure(1)
colormap(jet)
colorbar
[x,y,z]=meshgrid(ycoor(2:end),xcoor(1:end),zcoor(1:end));
h=slice(x,y,z,Ey_sample,[],0,[])
h_p=slice(x,y,z,Ey_sample,0,[],[])
% h_pp=slice(x,y,z,Ey_sample,[],[],0)
set(h, 'edgecolor', 'none')
set(h_p, 'edgecolor', 'none')
% set(h_pp, 'edgecolor', 'none')
caxis([-1.3 1.3])
title("Ey")
drawnow update

disp(time_step)
toc
end
end
end


