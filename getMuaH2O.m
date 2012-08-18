
function [ mu_a_H2O]=getMuaH2O(wavelength)


load H2O_extinction.mat

mu_as=2.303*H2O_extinction(:,2);
wavelengths2=H2O_extinction(:,1);

mu_a_H2O=0.1*interp1(wavelengths2,mu_as(:,1),wavelength);

