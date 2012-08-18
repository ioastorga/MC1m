%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For runing a Monte Carlo in C %%
%%     settings and storage      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
 
N = 700000;
%conc = input('ICG concentration from 6.5 uM to 1290 uM '); 
%ktimes = input('Number of runs per curve ');  
file_n = input('Enter file number ','s');
% ml_spheres = input('Volume of spheres per sample ');
% ml_sample = input('Sample volume');
% fv = ml_sphere*0.0952/ml_sample;   % given the density of polystirene 1.05g/ml

pol_dens = 1.05;
wa_dens = 1;
pc_solids = 0.1; % 10%
A = (wa_dens/pol_dens)*pc_solids;
sample_ml = 7.5;
ml_spheres = 0.3544;
fv = ml_spheres*A/sample_ml;

sd_d = 0.1;  % [um]
m_d = 1;       % [um]

sd_size = sd_d/2;  % [um]
m_size = m_d/2;       % [um]
max_size = (m_d+3*sd_d)/2;   % [um]
min_size = (m_d-3*sd_d)/2;  % [um]


% sd_size = 0.021/2;  % [um]
% max_size = 1.031/2;   % [um]
% min_size = 0.905/2;  % [um]
% m_size = 0.968/2;       % [um]

par_dir = date;
file_data = strcat('data',file_n,'_','int');
file_array = strcat('array',file_n,'_','int');
mkdir(par_dir,file_data)
mkdir(par_dir,file_array)

tic;

%Tissue radial geometry (for the array)
r_max = 8; % [mm] 
x2_max = 10;  % [mm] depth of 

%Fiber geometry
r3 = 0.125/2;                    % [mm] outer cladding
r2 = 0.105/2;                    % [mm] inner clading
r1 = 0.008/2;                    % [mm] core 
fiber_radii = [r1 r2 r3];
sensing_area = (pi*(r2^2-r1^2)); % Detection area pf the fiber

%Fiber numerical aperture
na_core = 0.12;
na_clad = 0.46;
na_fiber = [na_core,na_clad];

%Refractive index: Tissue
n_tissue = 1.4;          
%Refractive index: Fiber
n_core = 1.5;
n_clad = 1.46;
n_fiber = [n_core,n_clad];

%Extinction of the medium & operation bandwidth

min_wavelength = 790;    % [nm]  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_wavelength = 860;    % [nm]   BW reduced
num_wavelengths = 10;    % [nm]  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavelengths = linspace(min_wavelength,max_wavelength,num_wavelengths);
%   [mus_mie,g_mie] = microspheres_mie(fv, wavelengths);

%Tissue Image Array 
grid = 50;
grid0 = grid*r_max;
grid2 = grid*x2_max;
norm_matr_1 = linspace(0,r_max,grid0+1);
norm_matr_1=norm_matr_1(1:grid0); 
norm_matr = repmat(norm_matr_1,grid2,1);
width = 1/grid;
area = pi*width.*(2.*norm_matr+width);

%ICG concentrations
conc = 0;
        num_scatt=zeros(size(wavelengths))';
ka=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for w=wavelengths 
            ka=ka+1;
            kc = num2str(w);

            % Acceptance cone
            th_rand = rand(N,1);
            psi_rand = rand(N,1);
            th_max = asin(na_core/n_tissue);

            r_initial = sqrt((r1^2)*rand(N,1));
            psi = 2*pi*rand(N,1);
            th = acos(1-(1-cos(th_max))*rand(N,1));
            
            %Initial position
            x0 = r_initial.*cos(psi);
            x1 = r_initial.*sin(psi);
            x2 = zeros(N,1);

            %Initial direction
            u0 = sin(th).*cos(psi); 
            u1 = sin(th).*sin(psi);
            u2 = cos(th);

            %Optical paramenters of tissue
            %Optical paramenters of tissue
           [mu_a_H2O]=getMuaH2O(w);
           [Ua_icg]=getMuaICG_exp(conc,w);
           mua0 = mu_a_H2O + Ua_icg;   % [mm^-1]    
           
           nsteps = 7200;
           [mus_mie,g_mie, th, ph_ini] = microspheres_mie_ph_v2(fv, w,min_size,max_size,m_size,sd_size,nsteps);
           % **cumulative table** 
           % Toublanc D. Henyey-greenstein and mie phase functions 
           % in monte carlo radiative transfer computations. Appl Opt. 1996;35(18):3270-4.
           ph_tot = sum(ph_ini);
           ph0 = ph_ini/ph_tot;
           ph = cumsum(ph0);          
           
           
            %%%%%%%%% Monte Carlo en .c %%%%%%%%%   
            [tissue,fiber_tissue,r,depth,weight, path, z,num_scatt]=mc_1m_ph_sf(N,mua0,mus_mie,g_mie,x0,x1,x2,u0,u1,u2,...
                na_fiber,r_max,x2_max,n_tissue,n_fiber,fiber_radii,grid,th,ph); 
%            [mu_a_H2O]=getMuaH2O(w);
% %            if conc == 0 
% %                Ua_icg = 0;
% %            else
% %                [Ua_icg]=getMuaICG_conc2(conc,w);               
% %            end
%            [Ua_icg]=getMuaICG_water1(conc,w); 
%            mua_icg_w = mu_a_H2O + Ua_icg;   % [mm^-1]            
%            [mus_mie,g_mie] = microspheres_mie(fv, w,min_size,max_size,m_size,sd_size);  
%            
           mua(ka) =mua0;
           mus(ka) = mus_mie;
           g(ka) = g_mie;
% %           %%%%%%%%% Monte Carlo en .c %%%%%%%%%   
%             [tissue,fiber_tissue,r,depth,weight, path, z,num_scatt]=mc_1m_fc(N,mua_icg_w,mus_mie,g_mie,x0,x1,x2,u0,u1,u2,na_fiber,r_max,x2_max,n_tissue,n_fiber,fiber_radii,grid); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            inds0 = find(r);                    %Find detected photons        
            r_into_fiber = r(inds0);            %Find radial position of detected photons
            depth_into_fiber = depth(inds0);    %Find depth of detected photons
            weight_into_fiber = weight(inds0);  %Find weight of detected photons
            path_into_fiber = path(inds0);      %Find pathlength of detected photons
            z_into_fiber = z(inds0);            %Number of photons detected
            scatt_into_fiber = num_scatt(inds0);                    %k increments with wavelenght          
                        
            tissue_area = tissue./area;
            fiber_tissue_area = fiber_tissue./area;
            
            figure
            subplot(1,2,1)
            imagesc(0:r_max, 0:x2_max, log(tissue_area))
            set(gca,'DataAspectRatio',[1 1 1])
            colormap(hot)
            ylabel('Z [mm]')
            xlabel('r [mm]')
            wl = num2str(w);
            wav = strcat(wl,' nm');
            title(wav) 
            subplot(1,2,2)
            imagesc(0:r_max, 0:x2_max, log(fiber_tissue_area))
            set(gca,'DataAspectRatio',[1 1 1])
            colormap(hot)
            ylabel('Z [mm]')
            xlabel('r[mm]')
            wl = num2str(w);
            wav = strcat(wl,' nm');
            title(wav) 
            set(gcf, 'Visible', 'off') 

            file_tissue = strcat(par_dir,'/',file_array,'/tissue_',kc);    
            saveas(gcf, file_tissue, 'jpg')
            
            figure
            hist(scatt_into_fiber,1000);
            ylabel('N')
            xlabel('Number of scattering events')
            title('Scattering events') 
            set(gcf, 'Visible', 'off') 

            file_scatt = strcat(par_dir,'/',file_array,'/hist_',kc);    
            saveas(gcf, file_scatt, 'jpg')
    
            file_r_into_fiber = strcat(par_dir,'/',file_data,'/file_r_into_fiber_',kc);
            file_weight_into_fiber= strcat(par_dir,'/',file_data,'/weight_into_fiber',kc);
            file_scatt_into_fiber = strcat(par_dir,'/',file_data,'/file_scatt_into_fiber_',kc);

            save(file_r_into_fiber, 'r_into_fiber', '-ascii');
            save(file_weight_into_fiber,'weight_into_fiber','-ascii');
            save(file_scatt_into_fiber, 'scatt_into_fiber','-ascii');    


       end    % for wavelengths
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       file_wavelengths = strcat(par_dir,'/',file_data,'/file_wavelengths');
       save(file_wavelengths, 'wavelengths', '-ascii');
       file_mua_zero = strcat(par_dir,'/',file_data,'/file_mua');
       save(file_mua_zero, 'mua', '-ascii');
       file_mus_zero = strcat(par_dir,'/',file_data,'/file_mus');
       save(file_mus_zero, 'mus', '-ascii');
       file_g_zero = strcat(par_dir,'/',file_data,'/file_g');
       save(file_g_zero, 'g', '-ascii');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t=toc
