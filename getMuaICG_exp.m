function [Ua_icg]=getMuaICG_exp(conc,wavelength)
%WATER
%790-870 nm
load mua.mat;
if conc == 0 
   Ua_icg = 0;
else
wavelengths = mua(:,1);
mu_a_icg_10 = interp1(wavelengths,mua(:,2),wavelength);
mu_a_icg_40 = interp1(wavelengths,mua(:,3),wavelength);
mu_a_icg_100 = interp1(wavelengths,mua(:,4),wavelength);
mu_a_icg_200 = interp1(wavelengths,mua(:,5),wavelength);

def_conc = [0 10 40 100 200];
mu_a_icg =[0 mu_a_icg_10 mu_a_icg_40 mu_a_icg_100 mu_a_icg_200];
Ua_icg = interp1(def_conc, mu_a_icg, conc);
end
% end