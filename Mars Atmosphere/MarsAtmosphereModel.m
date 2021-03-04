clear all; close all; clc;

% import Mars data
% make sure excel file is in same folder as this matlab file 
marsdata = xlsread('Mars_data_all.xlsx');
mars_altitude = marsdata(1:end,2);       % km
mars_densities = marsdata(1:end,4);    % kg/m^3 
mars_temps = marsdata(1:end,5); % Kelvin
mars_pressures = marsdata(1:end,6); %N/m^2
%MARSGRAM data. tpdmsy11.txt & tpdloy11.txt Ls=180, Lat =7.5

n = 215000;
altitudes = zeros(n,1);
temperatures = zeros(n,1);
densities = zeros(n,1);
for i=1:length(altitudes)
    altitudes(i) = (i-1)/1000; %alt in km
end
for i=1:length(altitudes)
    alt = altitudes(i);
   if alt>220
       alt_m=alt*1000;  % convert km-->m
       rho = exp(2.65472e-11*(alt_m^5) - 2.45558e-08*(alt_m^4) + 6.31410e-06*(alt_m^3) + 4.73359e-04*(alt_m^2) - 0.443712*alt_m + 23.79408);
   else 
        rho = interp1(mars_altitude,mars_densities,alt,'makima','extrap'); 
   end
   T_CO2 = interp1(mars_altitude,mars_temps,alt,'makima','extrap');
   P_inf = interp1(mars_altitude,mars_pressures,alt,'makima','extrap');
   densities(i)=rho;
   temperatures(i)=T_CO2;
end
figure;
yyaxis left
semilogy(altitudes,densities)
ylabel('Density (kg/m^3)')
xlabel('Altitude (km)')
hold on
yyaxis right
plot(altitudes,temperatures,'--')
ylabel('Temperature (K)')
legend('Density','Temperature','Location','east')