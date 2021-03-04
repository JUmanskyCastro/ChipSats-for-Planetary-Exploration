clear all; close all; clc;

n = 500000;
altitudes = zeros(n,1);
temperatures = zeros(n,1);
densities = zeros(n,1);
for i=1:length(altitudes)
    altitudes(i) = (i-1)/1000; %alt in km
end
R_E = 6371;
for i=1:length(altitudes)
    alt = altitudes(i);
    zeta = (alt-120)*(6356.766+120)/(6356.766+alt);
    if alt>750
        A= -3.701195E-12; B= -8.608611E-09; C= 5.118829E-05; D= -0.06600998; E= -6.137674;
        rho = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % kg/m^3
        T_inf = 1000-640*exp(-0.01875*zeta); %K
        A=2.813255E-11; B=-1.120689E-07; C=1.695568E-04; D=-0.1188941; E=14.56718;
        P_inf = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % Pa
    elseif alt>500
        A= 8.105631E-12; B= -2.358417E-09; C= -2.635110E-06; D= -0.01562608; E= -20.02246;
        rho = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % kg/m^3
        T_inf = 1000-640*exp(-0.01875*zeta);%K
        A=-7.835161E-11;B=1.964589E-07;C=-1.657213E-04;D=0.04305869;E=-14.77132;
        P_inf = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % Pa
    elseif alt>300
        A= 1.140564E-10; B= -2.130756E-07; C= 1.570762E-04; D= -0.07029296; E= -12.89844;
        rho = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % kg/m^3
        T_inf = 1000-640*exp(-0.01875*zeta);%K
        A=9.814674E-11; B=-1.654439E-07; C=1.148115E-04; D=-0.05431334; E=-2.011365;
        P_inf = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % Pa
    elseif alt>200
        A= 1.199282E-09; B= -1.451051E-06; C= 6.910474E-04; D= -0.1736220; E= -5.321644;
        rho = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % kg/m^3
        T_inf = 1000-640*exp(-0.01875*zeta);%K
        A=8.113942E-10; B=-9.822568E-07; C=4.687616E-04; D=-0.1231710; E=3.067409;
        P_inf = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % Pa
    elseif alt>150
        A= 1.906032E-08; B= -1.527799E-05; C= 0.004724294; D= -0.6992340; E= 20.50921;
        rho = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % kg/m^3
        T_inf = 1000-640*exp(-0.01875*zeta);%K
        A=1.209434E-08; B=-9.692458E-06; C=0.003002041; D=-0.4523015; E=19.19151;
        P_inf = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % Pa
    elseif alt>120
        A= 3.661771E-07; B= -2.154344E-04; C= 0.04809214; D= -4.884744; E= 172.3597;
        rho = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % kg/m^3
        T_inf = 1000-640*exp(-0.01875*zeta);%K
        A=2.283506E-07; B=-1.343221E-04;C=0.02999016; D=-3.055446; E=113.5764;
        P_inf = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % Pa
    elseif alt>110
        A= 0.00000; B= -8.854164E-05; C= 0.03373254; D= -4.390837; E= 176.5294;
        rho = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % kg/m^3
        T_inf = 240+12*(alt-110);%K
        A=0.000000; B=-6.539316E-05;C=0.02485568;D=-3.223620;E=135.9355;
        P_inf = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % Pa
    elseif alt>100
        A= -1.240774E-05; B= 0.005162063; C= -0.8048342; D= 55.55996; E= -1443.338;
        rho = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % kg/m^3
        T_inf = 263.1905-76.3232*sqrt(1-((alt-91)/-19.9429)^2);%K
        A=0.000000; B=6.693926e-05; C=-0.01945388; D=1.719080; E=-47.75030;
        P_inf = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % Pa
    elseif alt>91
        A=0.00000; B= 2.873405E-05; C= -0.008492037; D= 0.6541179; E= -23.62010;
        rho = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % kg/m^3
        T_inf = 263.1905-76.3232*sqrt(1-((alt-91)/-19.9429)^2);%K
        A=0.000000;B=3.304895e-05;C=-0.009062730;D=0.6516698;E=-11.03037;
        P_inf = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % Pa
    elseif alt>86
        A=0.00000; B= -3.322622E-06; C= 9.111460E-04; D= -0.2609971; E= 5.944694;
        rho = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % kg/m^3
        T_inf = 186.8673;%K
        A=0.000000; B=2.159582E-06;C=-4.836957e-04;D=-0.1425192;E=13.47530;
        P_inf = exp(A*alt^4 + B*alt^3 + C*alt^2 + D*alt + E); % Pa
    else
        gAlt=(R_E*1000/(R_E*1000+1000*alt))*1000*alt; %solves for Geopotential Height in m
        [T_inf, a_inf, P_inf, rho] = atmoscoesa(gAlt, 'none');
    end
    temperatures(i) = T_inf;
    densities(i) = rho;
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