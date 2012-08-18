%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For runing a Monte Carlo in C %%
%%     settings and storage      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
 

N = input('Enter number of photons ');
%conc = input('ICG concentration from 6.5 uM to 1290 uM '); 
ktimes = input('Number of runs per curve ');  
file_n = input('Enter file number ','s');

par_dir = date;
file_data = strcat('data',file_n,'_','icg');
file_array = strcat('array',file_n,'_','icg');
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
fv=0.0045;
%fv=0.0085;

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
conc_icg = [0 10];
%conc_icg = [0 10 40 100 200];
%conc_icg = [0];
% c0 = 1200;
% c1 = 0.75*c0;
% c2 = 0.5*c0;
% c3 = 0.25*c0;
% conc_icg = [0 c3 c2 c1 c0];

% load mie_1u_780_880.mat
% scatt_w = mie_1u_780_880(:,1);
% scatt_mus = mie_1u_780_880(:,2);
% scatt_g = mie_1u_780_880(:,3);
% sd_size = 0.021/2;  % [um]
% max_size = 1.031/2;   % [um]
% min_size = 0.905/2;  % [um]
% m_size = 0.968/2;       % [um]
sd_d = 0.1;  % [um]
m_d = 1;       % [um]

sd_size = sd_d/2;  % [um]
m_size = m_d/2;       % [um]
max_size = (m_d+3*sd_d)/2;   % [um]
min_size = (m_d-3*sd_d)/2;  % [um]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for c = 1:1:length(conc_icg)
    conc=conc_icg(c);   
    kb = num2str(conc);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k_times = 1:ktimes
        k = 0;
        photons_area = zeros(length(wavelengths),1); 
        depth_area = zeros(length(wavelengths),1); 
        path_area = zeros(length(wavelengths),1); 
        ka = num2str(k_times);
        %ka = 3;
        out=zeros(size(wavelengths))';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for w=wavelengths 
            k=k+1;
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
           [mu_a_H2O]=getMuaH2O(w);
           [Ua_icg]=getMuaICG_exp(conc,w);
          %[Ua_icg]=getMuaICG_water1(conc,w);
           mua = mu_a_H2O + Ua_icg;   % [mm^-1]    
           
           nsteps = 7200;
           [mus_mie,g_mie, th, ph_ini] = microspheres_mie_ph_v2(fv, w,min_size,max_size,m_size,sd_size,nsteps);
           % **cumulative table** 
           % Toublanc D. Henyey-greenstein and mie phase functions 
           % in monte carlo radiative transfer computations. Appl Opt. 1996;35(18):3270-4.
           ph_tot = sum(ph_ini);
           ph0 = ph_ini/ph_tot;
           ph = cumsum(ph0);          
           
           
            %%%%%%%%% Monte Carlo en .c %%%%%%%%%   
            [tissue,fiber_tissue,r,depth,weight, path, z,num_scatt]=mc_1m_ph1(N,mua,mus_mie,g_mie,x0,x1,x2,u0,u1,u2,...
                na_fiber,r_max,x2_max,n_tissue,n_fiber,fiber_radii,grid,th,ph); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            inds0 = find(r);                     %Find detected photons        
            r_into_fiber = r(inds0);             %Find radial position of detected photons
            depth_into_fiber = depth(inds0);     %Find depth of detected photons
            weight_into_fiber = weight(inds0);   %Find weight of detected photons
            path_into_fiber = path(inds0);       %Find pathlength of detected photons
            z_into_fiber = z(inds0);             %Number of photons detected
            scatt_into_fiber = num_scatt(inds0); % k increments with wavelenght

            photons_r = 0;
            for k1=1:1:length(r_into_fiber) 
                photons_r = photons_r + weight_into_fiber(k1);           
            end

            photons_area(k) = photons_r/sensing_area;
            depth_r(k) = mean(depth_into_fiber);
            path_r(k) = mean(path_into_fiber);

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

            file_tissue = strcat(par_dir,'/',file_array,'/tissue_',kb,'_',ka,'_',kc);    
           % file_fiber_tissue = strcat(par_dir,'/array/fiber_tissue_',kb,'_',ka,'_',kc);    
           % save(file_tissue, 'tissue_area','-ascii');  
           % save(file_fiber_tissue, 'fiber_tissue_area','-ascii');  
            saveas(gcf, file_tissue, 'jpg')
            
            figure
            hist(scatt_into_fiber,100);
            ylabel('N')
            xlabel('Number of scattering events')
            title('Scattering events') 
            set(gcf, 'Visible', 'off') 

            file_scatt = strcat(par_dir,'/',file_array,'/hist_',kb,'_',ka,'_',kc);    
            saveas(gcf, file_scatt, 'jpg')

       end    % for wavelengths
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    file_reflectance = strcat(par_dir,'/',file_data,'/file_reflectance_',kb,'_',ka);
    file_depth= strcat(par_dir,'/',file_data,'/file_depth',kb,'_',ka);
    file_path = strcat(par_dir,'/',file_data,'/file_path_',kb,'_',ka);
    file_out = strcat(par_dir,'/',file_data,'/file_out_',kb,'_',ka);
    reflectance = photons_area/N;
    save(file_reflectance, 'reflectance', '-ascii');
    save(file_depth,'depth_r','-ascii');
    save(file_path, 'path_r','-ascii');    
    save(file_out, 'out','-ascii'); 

    end     % for k_times 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end    % for c:conc_icg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=toc
