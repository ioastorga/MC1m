
%%%%%mejor 3!

clc
clear
close all

min_wavelength = 790;    % [nm]  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_wavelength = 860;    % [nm]   BW reduced
% min_wavelength = 795;    % [nm]  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% max_wavelength = 865;    % [nm]   BW reduced
num_wavelengths = 10; 
wavelengths = linspace(min_wavelength,max_wavelength,num_wavelengths);
times = 3; 
directory = '15-Aug-2012/data3_icg';

tot_depth=[];
tot_depth0=[];
tot_path=[];
tot_path0=[];

j=['k','r','m','b','c','g','k-'];
%ICG concentrations
conc_icg = [0 40]; %data3_icg

r3 = 0.125/2;                    % [mm] outer cladding
r2 = 0.105/2;                    % [mm] inner clading
r1 = 0.008/2;                    % [mm] core 
fiber_radii = [r1 r2 r3];
sensing_area = (pi*(r2^2-r1^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for c = 1:1:length(conc_icg)
    conc=conc_icg(c);
    reflectance0 = zeros(num_wavelengths,times);
    kb = num2str(conc);      
        for k_times = 1:times
            ka = num2str(k_times);
          
            file_reflectance = strcat(directory,'/file_reflectance_',kb,'_',ka);
            file_depth= strcat(directory,'/file_depth',kb,'_',ka);
            file_path = strcat(directory,'/file_path_',kb,'_',ka);

            reflectance0(:,k_times) = load (file_reflectance);
            depth(:,k_times) = load (file_depth);
            path(:,k_times) = load (file_path);
        end                   
    mean_reflectance(:,c) = mean(reflectance0,2);
    e_ref (:,c)= std(reflectance0,0,2);
    mean_depth(:,c) = mean(depth,2);
    e_depth(:,c) = std(depth,0,2);
    mean_path (:,c)= mean(path,2);
    e_path(:,c) = std(path,0,2);        
end

    wav=wavelengths;
    ref_w = 800;      
    val = ref_w; % from wavelengths 
   
    [dum,bin]=histc(val,wavelengths); 
    index=bin+1; 
    if abs(val-wavelengths(bin))<abs(val-wavelengths(bin+1)) 
        fclosest=find(wavelengths(bin)) ;
        index=bin; 
    else 
        fclosest=find(wavelengths(index)) ;
    end 

   reflectance =  mean_reflectance/mean_reflectance(index,1); 
   % reflectance =  mean_reflectance; 
   for c = 1:1:length(conc_icg)
    figure(1)
    errorbar(wav,reflectance(:,c),e_ref(:,c),j(c));
    legend('off')
    ylabel('Reflectance')
    xlabel('Wavelength [nm]')
    axis([min_wavelength max_wavelength 0 1.1])
    grid on
    hold on

    figure(2)
    errorbar(wav,mean_depth(:,c),e_depth(:,c),j(c));
    ylabel('Mean depth[mm]')
    xlabel('Wavelength [nm]')
    %axis([790 870 1.68 1.86])
    grid on
    hold on

    figure(3)
    errorbar(wav,mean_path(:,c),e_path(:,c),j(c));
    ylabel('Mean Pathlength [mm]')
    xlabel('Wavelength [nm]')
    %axis([790 870 1.68 1.86])
    grid on
    hold on
   end
