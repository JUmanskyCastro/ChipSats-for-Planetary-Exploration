% Earth Thermal Single 
% using Hunter's orbital dynamics
% with Justin's thermal equations 
% last editted: 1/21/2021

%works. matches justin's ch4 results
%TODO test with ch 2 results
%TODO implement more ambitious heating model

%1/12 conclusions:
%it's fine as is for 1cmx1cm. all other paramaters seem pretty insignificant
%this is because peak temps happen while Kn>>10, so ST is just 1
%on initial testing, it seems that only gamma is important.

%for more realistic chipsats, Stanton number is in the transitional phase
%during temp spike
% molecular mass has a small impact. 
% but the Re2 number definitely matters.
% Need to calculate that properly :( - Calculated! Actually seems pretty
% consistent at ~2.72 for air. See if the same holds for the others
%it goes out of wack at slower velocities, but this is after the peak temp
% to be safe, the Re2 parmaeter ia now calculated in the initial
% conditions!
%
%for other planetary bodies:
% With ST, Kn, and Re2 calcs fixed, this can now be safely adapted to other
% planets. Pressure model not needed anymore, and mm seems pretty consistent. but can throw in laters.
%
% just make sure u have the right mm, gamma, and diam (unused now)
% Equilibrium temperatures: mars is 210K, titan is 85K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% parameters
h = 50/1000/1000;                     % spacecraft height, km, 35
w = 50/1000/1000;                     % spacecraft width, km
%d = 2/1000/1000;                      % spacecraft depth, km
p.l = h;
%p.Az = [h*d, h*d, h*h]*(10^-12);      % z-direction of coil surface, km^2
%p.A = h*h;                    % surface area, km^2
p.Cd_FM = 2.67;        % drag coefficient 2.67 face on, 1.43 averaged, in FM
p.Cd_SS = 1.28;  % drag coefficient 1.28 face on, 0.8155 averaged, in subsonic
p.omega = [0, 0, 7.292e-5]';  % angular velocity of Earth, rad/s 
p.mu_e = 398600.4418;         % Earth gravitational parameter, km^3/s^2
p.J2 = 1.7555e10;             % J2, km^5/s^2
p.R_E = 6371;                 % radius of Earth, km
G = 6.67e-11/(10^9);               % Earth's gravitational constant, km^3/(kg*s^2)
M_E = 5.972e24;             % mass of Earth, kg
initial_alt = 350;   % km, 374
v_orbital = sqrt(G*M_E/(p.R_E+initial_alt));    %km/s, initial orbital velocity
inclination = deg2rad(50);  % around ISS orbital inclination, radians

%ratio of specific heat
p.gammaAir = 1.4;
p.R = 8.314; %J/(mol*K) Universal gas constant
p.kb = 1.380649e-23; %boltzmann constant J/K

%viscosity model constants
%get temperature vs. altitude profile
%use alt in ode15s to find temp, use temp to find dyn visc.
p.mu0_air = 1.716e-5; %reference viscosity in centipoise at reference temperature To
p.S_air = 111; %Sutherland's constant
p.T0_air = 273; %reference temperature in degreees Kelvin


%Re2 parameter calculation
gammaAir = p.gammaAir;
mm = 28.97/1000; %/N_A; %kg/mol molar mass of air
a_inf = sqrt(gammaAir*p.R*255/mm); %m/s
M1 = 1000*v_orbital/a_inf;
Mp2 = sqrt(2/(gammaAir-1))*(M1^2 - 1)/(sqrt(1+(gammaAir-1)/2*M1^2)*sqrt(2*gammaAir/(gammaAir-1)*(M1^2)-1));
p.Re2_param = Mp2*sqrt(pi*gammaAir/2); 
%this gets used in the ODE function as Re2 = Re2_param/Kn2


% Justin's parameters
p.N_A = 6.022*10^23; %Avogadro's number
p.Cp = 1090; %J/(kg K) specific heat 678 silicon wafer, 1090 kapton, about same for FR4
p.epsilon = 0.85; %emissivity
p.sbsigma = 5.67e-8; %stephen-boltzman constant W/(m^2 K^4)
p.diam = 0.3e-9; %mean diameter of air molecule (m)
% Justin's initial conditions
Q_aero0 = 0;
Q_rad0 = 0;
T0 = 250;


% US standard atmosphere viscosity model
% https://www.engineeringtoolbox.com/standard-atmosphere-d_604.html
% p.heights = [-1000 0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 15000 ...
%     20000 25000 30000 40000 50000 60000 70000 80000 400000]';
% p.dyn_visc = [1.821 1.789 1.758 1.726 1.694 1.661 1.628 1.595 1.561 1.527 ...
%     1.493 1.458 1.422 1.422 1.448 1.475 1.601 1.704 1.584 1.438 1.321 0]';
% for i=1:length(p.dyn_visc)
%     p.dyn_visc(i)=p.dyn_visc(i)*10^-5;
%     p.heights(i)=p.heights(i)/1000;
% end

%p.visc = polyfit(height,dyn_visc,3);

% MSISE-90 Model of Earth's Upper Atmosphere
% http://www.braeunig.us/space/atmos.htm
p.height_mm = zeros(26,1); %km
for i=1:length(p.height_mm)
    p.height_mm(i)=(i-1)*20;
end
p.molmass = [28.9502,28.9502,28.9502,28.9502,29.0175,27.7137,25.8745,24.5349,23.4225,22.4106, ...
    21.4734, 20.6108, 19.8292, 19.1337, 18.5256, 18.0015, 17.5537, 17.1721, ...
    16.8449, 16.5597, 16.3044, 16.0669, 15.8360, 15.6008, 15.3508, 15.0760]'./1000;

%p.mm_f = polyfit(height_mm,molmass,3);

% % density model
% p.altitudes = 0:20:500;      %altitudes, km
% p.densities = [1.17, 9.49e-2, 4.07e-3, 3.31e-4, 1.68e-5, 5.08e-7, 1.80e-8,...
%     3.26e-9, 1.18e-9, 5.51e-10, 2.91e-10, 1.66e-10, 9.91e-11, 6.16e-11, ...
%     3.94e-11, 2.58e-11, 1.72e-11, 1.16e-11, 7.99e-12, 5.55e-12, 3.89e-12, ...
%     2.75e-12, 1.96e-12, 1.4e-12, 1.01e-12, 7.3e-13]*(1000^3);   %kg/km^3
% 
% alt = zeros(5001,1);
% rho = zeros(5001,1);
% 
% for i = 1:length(alt)
%     alt(i)=(i-1)*0.1;
%     rho(i) = interp1(p.altitudes,p.densities,alt(i));
% end

% % plot of altitudes and densities
% figure
% semilogy(alt,rho)
% xlabel('altitude (km)')
% ylabel('rho (kg/km^3)')
% title('Atmospheric Density with Altitude')
% 
% figure
% plot(altitudes,densities)
% xlabel('altitude (km)')
% ylabel('rho (kg/km^3)')
% title('Atmospheric Sample Data')

% initial conditions
totaltime = linspace(0, 10*24*3600, 10*24*3600);     % total time

accuracy = 1e-6;
options = odeset('AbsTol', accuracy, 'RelTol', accuracy,'Events', @stopEventsFxn);

n = 1;      %number of simulations,100
vd = 1/1000;     % deployment velocity, km/s
figure;       % create figure for batch plotting
r = zeros(length(totaltime),n);
x = zeros(length(totaltime),n); y = zeros(length(totaltime),n); z = zeros(length(totaltime),n);
xpos = zeros(length(totaltime),n); ypos = zeros(length(totaltime),n);
xdot = zeros(length(totaltime),n); ydot = zeros(length(totaltime),n); zdot = zeros(length(totaltime),n);
Q_aero = zeros(length(totaltime),n); Q_rad = zeros(length(totaltime),n); T = zeros(length(totaltime),n);
dvdot = zeros(length(totaltime),n); time = zeros(length(totaltime),n);
mass = zeros(length(n),1);
area = zeros(length(n),1);
Tmax = zeros(length(n),1);
TimeToGround = zeros(length(n),1);
lat = zeros(1,n); long = zeros(1,n);

% varying simulation parameters
for i = 1:n
    i           % tells us what simulation number we are on
    %area(i) = h*w*normrnd(0.6,0.01);  % km^2, tumbling average area (1 = worse case, face on)
    area(i) = h*w;
    p.Area_surf = 2*area(i);
    p.Area = area(i); %YOU NEED TO FIND A GOOD BACKING FOR THIS RATIO
    q = 0; %2*pi/n*i;
    %mass(i) = normrnd(3/1000,0.1/1000);                 % mass of spacecraft, kg 
    mass(i) = 3/1000; %kg, 5.8250e-06 for justin's 1cmx1cm chipsat
    p.m = mass(i);
    
    initial_position = [normrnd(initial_alt,0.1/1000), 0, 0];     % km, randomizing height of chipsat in deployer
    initial_position = (initial_position/norm(initial_position))*(p.R_E+norm(initial_position));
    %initial_velocity = [0, v_orbital*cos(inclination)+vd*cos(q), v_orbital*sin(inclination)+vd*sin(q)];   %km/s
    initial_velocity = [0, v_orbital, 0];   %km/s
    initial_dv = [0];                            %km/s^2
    %initial_A = A/(1000^2);                      % initial area, m^2
    w0 = [initial_position'; initial_velocity'; initial_dv; Q_aero0; Q_rad0; T0];

    f = @(t,w)  rhs(t,w,p);
    [tarray zarray] = ode15s(f, totaltime, w0, options);
    [~,U] = cellfun(@(t,w)  rhs(t,w.',p), num2cell(tarray), num2cell(zarray,2),'uni',0);
    U = cell2mat(U);
    
    for j = 1:length(tarray)
        x(j,i) = zarray(j,1);
        y(j,i) = zarray(j,2);
        z(j,i) = zarray(j,3);
        xdot(j,i) = zarray(j,4);
        ydot(j,i) = zarray(j,5);
        zdot(j,i) = zarray(j,6);
        dvdot(j,i) = zarray(j,7);
        Q_aero(j,i) = zarray(j,8); 
        Q_rad(j,i) = zarray(j,9); 
        T(j,i) = zarray(j,10);
        r(j,i) = sqrt(x(j,i)^2+y(j,i)^2+z(j,i)^2)-p.R_E;
        time(j,i) = tarray(j);
        
        if x(j,i)>=0
            xpos(j,i) = x(j,i)-(p.R_E*abs(x(j,i))/(r(j,i)+p.R_E));
        else
            xpos(j,i) = x(j,i)+(p.R_E*abs(x(j,i))/(r(j,i)+p.R_E));
        end
        if y(j,i)>=0
            ypos(j,i) = y(j,i)-(p.R_E*abs(y(j,i))/(r(j,i)+p.R_E));
        else
            ypos(j,i) = y(j,i)+(p.R_E*abs(y(j,i))/(r(j,i)+p.R_E));
        end
       

    end
    
    index = find(r(:,i)<0.1 & r(:,i)>0);     % find index when chip is at 0 km
    last = index(end);
    ECI = [x(last,i) y(last,i) z(last,i)]; % final 0 km altitude states
    ECIx = ECI(1); ECIy = ECI(2); ECIz = ECI(3);
    lat(i) = rad2deg(atan2(ECIz,sqrt(ECIx^2+ECIy^2)));  %deg
    long(i) = rad2deg(atan2(ECIy,ECIx));    %deg
    
    Tmax(i) = max(zarray(:,10))-273; %max temp, C
    TimeToGround(i) = length(tarray); %s
    
   % plotting Hunter & Justin
    yyaxis left
    plot(time(1:last,i)/(24*3600),r(1:last,i))
    ylabel('Altitude (km)')
    xlabel('Time (days)')
    hold on
    yyaxis right
    plot(time(1:last,i)/(24*3600),T(1:last,i)-273)
    ylabel('Temperature (C)')
    title('ChipSat Atmospheric Entry')
 
end

% figure
% plot(tarray/(24*3600),U)

% plot on world map 
figure;
worldmap('World')
load coastlines
plotm(coastlat,coastlon)
scatterm(lat,long,'r','filled')
title('ChipSat Landing Locations')

figure;
axis equal
hold on
for i=1:n
    plot(xpos(:,i),ypos(:,i));
end
title('ChipSat Deorbit Trajectories')

%     
% %r = sqrt(x^2+y^2+z^2)-p.R_E;        %altitude
% 

% figure(2)
% plot(tarray/86400,x)
% title('Orbital Decay from drag')
% xlabel('Time (Days)')
% ylabel('Height (km)')
% 
% figure
% semilogy(tarray,dvdot*1000)
% title('Delta V from drag')

figure;
yyaxis left
plot(tarray/(3600),U(:,1))
ylabel('Mach number')
hold on
yyaxis right
plot(time(1:last,i)/(3600),T(1:last,i)-273,'--')
ylabel('Temperature (C)')
xlabel('Time (hours)')
title('Mach Number vs. Time');

figure;
plot(tarray/(24*3600),U(:,2))
title('Stanton Number & Temp vs. Time');
hold on
yyaxis right
plot(time(1:last,i)/(24*3600),T(1:last,i)-273)
ylabel('Temperature (C)')

figure;
loglog(U(:,3),U(:,2))
title('Stanton Number vs. Knudsen Number');

figure;
plot(U(:,4),U(:,5))
title('T_inf vs. Altitude');
xlabel('Temperature (K)')
ylabel('Altitude (km)')

figure;
semilogy(tarray/3600,U(:,6)) %Re
hold on
semilogy(tarray/3600,U(:,3)) %Kn
title('Re & Kn vs. Time(hr)');

figure;
semilogy(tarray/3600,U(:,7)) %Re2
title('Re2 vs. Time(hr)')
hold on
yyaxis right
plot(time(1:last,i)/(3600),T(1:last,i)-273)
ylabel('Temperature (C)')

figure;
semilogy(tarray/3600,U(:,8)) %Re2 const.
title('Re2 parameter vs. Time(hr)')
hold on
yyaxis right
plot(time(1:last,i)/(3600),T(1:last,i)-273)
ylabel('Temperature (C)')

figure;
semilogy(tarray/3600,U(:,9)) %Cd
title('Drag Coeff. vs. Time(hr)')
hold on
yyaxis right
plot(time(1:last,i)/(3600),T(1:last,i)-273)
ylabel('Temperature (C)')

function [wdot,u] = rhs(t,w,p) %[zdot, u] if you want to export additional variables

    %unpack parameters
    names = fieldnames(p);
    for i=1:length(names)
      eval([names{i} '= p.' names{i} ';' ]); % e.g.  a = p.a;
    end
    
    %unpack w
    x = w(1); y = w(2); z = w(3); xdot = w(4); ydot = w(5); zdot = w(6); dv = w(7);
    Q_aero = w(8); Q_rad = w(9); T = w(10); 
    
    r = [x; y; z];  %km
    rdot = [xdot; ydot; zdot]; %km/s
    
    Vs = rdot - cross(omega,r); %km/s 
    
    %A_eff = Az*Vs/abs(Vs);  %km^2
    
%     rho= interp1(altitudes,densities,norm(r)-R_E);   %kg/km^3

   alt = norm(r)-p.R_E    %km
   zeta = (alt-120)*(6356.766+120)/(6356.766+alt);
% new density model
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
    
    rho_m = rho;    % kg/m^3
    rho = rho*10^9; %conversion to kg/km^3
     
    %CHECK UNITS!!!!!!
    %rho = 1;
   
    % Justin's thermal code
     %mu = polyval(visc,alt); % dynamic viscosity (N s/m^2)
     
     mu = mu0_air*((T_inf/T0_air)^3/2)*(T0_air+S_air)/(T_inf+S_air);
     
     %mu = interp1(heights,dyn_visc,alt,'makima','extrap'); %N s/m^2
     mm = interp1(height_mm,molmass,alt,'makima','extrap'); %kg/mol
    
     L = l*1000; %characteristic length, m
     Re = rho_m*1000*norm(Vs)*L/mu;
     %mm = 28.97/1000; %/N_A; %kg/mol molar mass of air
     %mm = polyval(mm_f,alt); %kg/mol
     a_inf = sqrt(gammaAir*R*T_inf/mm); %m/s
     Ma = 1000*norm(Vs)/a_inf;
     %zeta_a = R*255/(sqrt(2)*pi*(diam^2)*N_A*P_inf); %mean free path
     zeta_a = mu/rho_m*sqrt(pi*mm/(N_A*2*kb*T_inf));
     Kn = zeta_a/L;
     %Kn = sqrt(gammaAir*pi/2)*Ma/Re; %consistent with above, not below
     %Kn = kb*T_inf/(sqrt(2)*diam^2*P_inf*L);
     %Re2 = 3.33/Kn; %works for air. more generally, see below
    
     rho0_rho1 = (1+(gammaAir-1)/2*Ma^2)^(1/(gammaAir-1));
     rho1_rho2 = ((gammaAir-1)*(Ma^2)+2)/(gammaAir+1)/Ma^2;
     Kn2 = Kn*rho1_rho2*rho0_rho1;
     Re2 = Re2_param/Kn2;

%      Mp2 = sqrt(2/(gammaAir-1))*(M1^2 - 1)/(sqrt(1+(gammaAir-1)/2*M1^2)*sqrt(2*gammaAir/(gammaAir-1)*(M1^2)-1));
%      Re2 = abs(Mp2)/Kn2*sqrt(pi*gammaAir/2); 
     
     %Re2 = sqrt(pi*gammaAir/2)/Kn*rho2_rho1*rho1_rho0*Mp2;
     %Mp2*sqrt(pi*gammaAir/2)/rho1_rho2/rho0_rho1
     %next steps: find the more general form of the Re2 formula - done
     %if ambitious, add in hunter's full heat transfer equations
     %good news is your results are starting to match!
     
      ST_FM = 1;
      Cs = 1/sqrt(2);
      ST_C = Cs*2.1/sqrt(Re2);
      %C_FM = 2;
      %if Re>2300
      %    Cc = 0.0592*Re^(-0.2);
      %else
      %    Cc = 0.664*Re^(-0.5);
      %end
      %Cd = Cd_SS/sqrt(1+(Cd_SS/Cd_FM)^2);
      

      if Kn>10
          ST = ST_FM;
          Cd = Cd_FM;
          %Cskin = C_FM;
      elseif Kn>0.01
          ST = ST_C/sqrt(1+(ST_C/ST_FM)^2);
          %Cskin = C_FM/(1+(Cc/C_FM)^2);
          Cd = Cd_FM;
      else
          ST = ST_C;
          %Cskin = Cc;
          Cd = Cd_SS/sqrt(1+(Cd_SS/Cd_FM)^2);
          if Cd < Cd_SS
              Cd=Cd_SS;
          end
      end
    
      %ST=0.1;
    Qdot_aero = 1/2*ST*rho_m*(1000*norm(Vs))^3*(Area*1000^2); %might have to change sdot to velocity resultant vector magnitude
    Qdot_rad = sbsigma*epsilon*(T^4-255^4)*(Area_surf*1000^2); 
    Tdot = (Qdot_aero - Qdot_rad)/(m*Cp);
    
    drag =(1/2)*Cd*(Area)*rho*((Vs')*Vs)/m.*(-Vs)/norm(Vs);  % deceleration due to drag, km/s^2
    % Hunter's code did not include Cd
    
    dvdot = norm(drag); %km/s^2
    
%     if norm(rdot)>0.35
    xddot = -mu_e*x/((x^2 + y^2 + z^2)^(3/2))  + (5*J2*x*(-x^2-y^2+2*z^2))/(2*(x^2+y^2+z^2)^(7/2)) + J2*x/(x^2+y^2+z^2)^(5/2) + drag(1);
    yddot = -mu_e*y/((x^2 + y^2 + z^2)^(3/2))  + (5*J2*y*(-x^2-y^2+2*z^2))/(2*(x^2+y^2+z^2)^(7/2)) + J2*y/(x^2+y^2+z^2)^(5/2) + drag(2);
    zddot = -mu_e*z/((x^2 + y^2 + z^2)^(3/2))  + (5*J2*z*(-x^2-y^2+2*z^2))/(2*(x^2+y^2+z^2)^(7/2)) - 2*J2*z/(x^2+y^2+z^2)^(5/2) + drag(3);
     

    wdot = [xdot; ydot; zdot; xddot; yddot; zddot; dvdot; Qdot_aero; Qdot_rad; Tdot];
    u = [Ma, ST, Kn, T_inf,alt,Re,Re2,Re2*Kn2,Cd];
end 

% finalpos = zeros(n,3); 

function [position,isterminal,direction] = stopEventsFxn(t,w)
%Stops integration after reaching altitude of 0
x = w(1); y = w(2); z = w(3);
r = [x; y; z];  %km
altitude = norm(r)-6371;        % hard coded earth radius 
position = altitude<0; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end