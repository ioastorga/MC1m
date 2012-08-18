function [Ua_icg]=getMuaICG_water1(conc,wavelength)
%WATER
%600 - 900 nm
load icg_extinction.mat;
if conc == 0 
   Ua_icg = 0;
else
wavelengths=icg_extinction(:,1);
mu_a_icg_6_5=0.1*log(10)*0.0000065*interp1(wavelengths,icg_extinction(:,2),wavelength);
mu_a_icg_65=0.1*log(10)*0.000065*interp1(wavelengths,icg_extinction(:,3),wavelength);
mu_a_icg_650=0.1*log(10)*0.00065*interp1(wavelengths,icg_extinction(:,4),wavelength);
mu_a_icg_1290=0.1*log(10)*0.001290*interp1(wavelengths,icg_extinction(:,5),wavelength);

def_conc = [0 6.5 65 650 1290];
mu_a_icg =[0 mu_a_icg_6_5 mu_a_icg_65 mu_a_icg_650 mu_a_icg_1290];
Ua_icg = interp1(def_conc, mu_a_icg, conc);
end
% end