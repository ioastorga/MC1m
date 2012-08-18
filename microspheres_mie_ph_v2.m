function [mu_s,g,th, ph_ini] = microspheres_mie_ph_v2(fv, wavelength, min_size,max_size,m_size,sd_size,nsteps)

n_med= 1.3316;
lambda = wavelength/1000;

p_size = linspace(min_size,max_size,100);
prob_fv1 = normpdf(p_size,m_size,sd_size);
ya=sum(prob_fv1);

% Polystirene
% http://refractiveindex.info/?group=PLASTICS&material=PS
% S. N. Kasarova et al. Analysis of the dispersion of optical plastic materials,
% Optical Materials 29, 1481-1490 (2007) doi:10.1016/j.optmat.2006.07.010
n_scat=sqrt(2.610025-6.143673*10^(-2)*lambda^2-1.312267*10^(-1)*lambda^(-2)+6.865432*10^(-2)*...
    lambda^(-4)-1.295968*10^(-2)*lambda^(-6)+9.055861*10^(-4)*lambda^(-8));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = n_scat/n_med;
ca =0;
gi = zeros(size(p_size));
mu_si = zeros(size(p_size));
phfi = zeros([nsteps,length(p_size)]);
for a = p_size
    ca = ca+1;
    prob_fv = normpdf(a,m_size,sd_size);
    rho = (prob_fv/ya*fv)/((4/3)*pi*a^3);  %volume density
    x = a*2*pi/(lambda/n_med);
    result= Mie(m, x);  
    gi(ca) = result(8);
    qsca = result(5);
    A = pi*a^2;
    sigma_s = qsca*A;
    musica = rho*sigma_s;
    mu_si(ca) = musica;
    th_phi = Mie_tetascanio(m, x,nsteps);
    ph_inter = th_phi(1:nsteps,2);
    phfi(:,ca) = ph_inter;

end
th = th_phi(1:nsteps,1);


%%Schmitt JM, Kumar G. Optical scattering properties of soft tissue: A discrete particle model. 
%%Appl Opt. 1998;37(13):2788-97.
num = sum(mu_si.*gi);
den = sum(mu_si);
g = num/den;
mu_s= den*1000;
ph_ini = sum(phfi,2)/den;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






