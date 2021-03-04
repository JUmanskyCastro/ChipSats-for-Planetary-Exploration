clear all; close all; clc;

% import Titan data
% make sure excel file is in same folder as this matlab file 
titandata = xlsread('yel_rec_excel.xlsx');
titan_altitude = titandata(1:end,1);       % km
titan_densities = titandata(1:end,2).*(10^3);    % g/cm^3 --> kg/m^3
titan_temps = titandata(1:end,3);    % temp (K)
titan_pressures = titandata(1:end,5);    % Pressure (Pa)
titan_mm = titandata(1:end,7);    % molar mass (kg/mol)

n = 120000;
altitudes = zeros(n,1);
temperatures = zeros(n,1);
densities = zeros(n,1);
for i=1:length(altitudes)
    altitudes(i) = (i-1)/100; %alt in km
end

for i=1:length(altitudes)
    alt = altitudes(i);
    rho= interp1(titan_altitude,titan_densities,alt);   %kg/m^3
    T_N2 = interp1(titan_altitude,titan_temps,alt,'makima','extrap');
    P_inf = interp1(titan_altitude,titan_pressures,alt,'makima','extrap');
    mm = interp1(titan_altitude,titan_mm,alt,'makima','extrap');
    
    densities(i) = rho;
    temperatures(i) = T_N2;
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