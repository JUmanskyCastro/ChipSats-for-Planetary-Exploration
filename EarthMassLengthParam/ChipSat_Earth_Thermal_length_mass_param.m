% Hunter Earth Thermal Batch 
% using orbital dynamics
% Earth densities
% with Justin's thermal equations 
% last editted: 12/15/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% parameters
%h = 50/1000/1000;                     % spacecraft height, km, 35
%w = 50/1000/1000;                     % spacecraft width, km
%d = 2/1000/1000;                      % spacecraft depth, km
%p.Az = [h*d, h*d, h*h]*(10^-12);      % z-direction of coil surface, km^2
%p.A = h*h;                    % surface area, km^2
p.Cd_FM = 2.67;              % drag coefficient 2.67 face on, 1.43 averaged, in FM
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

% Justin's parameters
p.N_A = 6.022*10^23; %Avogadro's number
p.Cp = 1090; %J/(kg K) specific heat 678 silicon wafer, 1090 kapton, about same for FR4
p.epsilon = 0.85; %emissivity
p.sbsigma = 5.67e-8; %stephen-boltzman constant
p.diam = 0.3e-9; %mean diameter of air molecule (m)

% Justin's initial conditions
Q_aero0 = 0;
Q_rad0 = 0;
T0 = 250;

%viscosity model constants
%get temperature vs. altitude profile
%use alt in ode15s to find temp, use temp to find dyn visc.
p.mu0_air = 1.716e-5; %reference viscosity in centipoise at reference temperature To
p.S_air = 111; %Sutherland's constant
p.T0_air = 273; %reference temperature in degreees Kelvin


%Re2 parameter calculation
mm = 28.97/1000; %/N_A; %kg/mol molar mass of air
a_inf = sqrt(p.gammaAir*p.R*255/mm); %m/s
M1 = 1000*v_orbital/a_inf;
Mp2 = sqrt(2/(p.gammaAir-1))*(M1^2 - 1)/(sqrt(1+(p.gammaAir-1)/2*M1^2)*sqrt(2*p.gammaAir/(p.gammaAir-1)*(M1^2)-1));
p.Re2_param = Mp2*sqrt(pi*p.gammaAir/2); 
%this gets used in the ODE function as Re2 = Re2_param/Kn2

% MSISE-90 Model of Earth's Upper Atmosphere
% http://www.braeunig.us/space/atmos.htm
p.height_mm = zeros(26,1); %km
for i=1:length(p.height_mm)
    p.height_mm(i)=(i-1)*20;
end
p.molmass = [28.9502,28.9502,28.9502,28.9502,29.0175,27.7137,25.8745,24.5349,23.4225,22.4106, ...
    21.4734, 20.6108, 19.8292, 19.1337, 18.5256, 18.0015, 17.5537, 17.1721, ...
    16.8449, 16.5597, 16.3044, 16.0669, 15.8360, 15.6008, 15.3508, 15.0760]'./1000;



% initial conditions
totaltime = linspace(0, 6*30*24*3600, 6*30*24*3600);     % total time (6 months)

accuracy = 1e-6;
options = odeset('AbsTol', accuracy, 'RelTol', accuracy,'Events', @stopEventsFxn);

n = 20;      %number of simulations,100
vd = 0; %1/1000;     % deployment velocity, km/s
%figure;       % create figure for batch plotting
%r = zeros(n,n,length(totaltime));
%x = zeros(n,n,length(totaltime)); y = zeros(n,n,length(totaltime)); z = zeros(n,n,length(totaltime));
%xdot = zeros(n,n,length(totaltime)); ydot = zeros(length(totaltime),n); zdot = zeros(n,n,length(totaltime));
%Q_aero = zeros(n,n,length(totaltime)); Q_rad = zeros(length(totaltime),n); T = zeros(n,n,length(totaltime));
%dvdot = zeros(n,n,length(totaltime)); time = zeros(n,n,length(totaltime));
mass = zeros(n,1);
area = zeros(n,1);
Tmax = zeros(n,n);
lat = zeros(n,n); long = zeros(n,n);
B_coeff = zeros(n,n);
TimeToGround = zeros(n,n);

% varying simulation parameters
for i = 1:n
    for j = 1:n
        (i-1)*n+j           % tells us what simulation number we are on
        %area(i) = h*w*normrnd(0.6,0.01);  % km^2, tumbling average area (1 = worse case, face on)
        p.l = i/2/100/1000; %increments chipsat length by 0.1cm, then converts to km
        area(i) = (p.l)^2;
        p.Area = area(i);
        p.Area_surf = 2*area(i);
        q = 0; %2*pi/n*i;
        mass(j) = (j^3)/1000/1000; %increments chipsat mass by 0.1grams, then converts to kg 
        p.m = mass(j);
        B_coeff(i,j) = mass(j)/(p.Cd_FM*area(i)*10^6);

        initial_position = [initial_alt, 0, 0];     % km, randomizing height of chipsat in deployer
        initial_position = (initial_position/norm(initial_position))*(p.R_E+norm(initial_position));
        initial_velocity = [0, v_orbital*cos(inclination)+vd*cos(q), v_orbital*sin(inclination)+vd*sin(q)];   %km/s
        %initial_dv = [0];                            %km/s^2
        %initial_A = A/(1000^2);                      % initial area, m^2
        w0 = [initial_position'; initial_velocity'; Q_aero0; Q_rad0; T0];

        f = @(t,w)  rhs(t,w,p);
        [tarray zarray] = ode15s(f, totaltime, w0, options);
        %[~,U] = cellfun(@(t,w)  rhs(t,w.',p), num2cell(tarray), num2cell(zarray,2),'uni',0);
        %U = cell2mat(U);
% 
        %for k = 1:length(tarray)
%             x(i,j,k) = zarray(k,1);
%             y(i,j,k) = zarray(k,2);
%             z(i,j,k) = zarray(k,3);
%             xdot(i,j,k) = zarray(k,4);
%             ydot(i,j,k) = zarray(k,5);
%             zdot(i,j,k) = zarray(k,6);
%             
%             Q_aero(i,j,k) = zarray(k,7); 
%             Q_rad(i,j,k) = zarray(k,8); 
%             T(i,j,k) = zarray(k,9);
            %r(i,j,k) = sqrt(zarray(k,1)^2+zarray(k,2)^2+zarray(k,3)^2)-p.R_E;
%             time(i,j,k) = tarray(k);
        %end

        %index = find(r(i,j,:)<1 & r(i,j,:)>0);     % find index when chip is at 0 km
        last = length(tarray);
        ECI = [zarray(last,1) zarray(last,2) zarray(last,3)]; % final 0 km altitude states
        ECIx = ECI(1); ECIy = ECI(2); ECIz = ECI(3);
        lat(i,j) = rad2deg(atan2(ECIz,sqrt(ECIx^2+ECIy^2)));  %deg
        long(i,j) = rad2deg(atan2(ECIy,ECIx));    %deg
        Tmax(i,j) = max(zarray(:,9))-273; %max temp, C
        TimeToGround(i,j) = length(tarray);
       % plotting Hunter & Justin
       % yyaxis left
       % plot(time(1:last,i)/(24*3600),r(1:last,i))
       % ylabel('Altitude (km)')
       % xlabel('Time (days)')
       % hold on
       % yyaxis right
       % plot(time(1:last,i)/(24*3600),T(1:last,i)-273)
       % ylabel('Temperature (C)')
       % title('ChipSat Atmospheric Entry')
    end
end

% figure
% plot(tarray/(24*3600),U)

% plot on world map 
figure;
worldmap('World')
load coastlines
plotm(coastlat,coastlon)
hold on
for j=1:n
    scatterm(lat(:,j),long(:,j),'r','filled')
end
hold off

lengths = zeros(n,1);
masses = zeros(n,1);
for i=1:n
    lengths(i) = i/2; %cm
    masses(i) = (i^3)/1000; %grams
end

%3D parametric plot
figure;
surf(masses,lengths, Tmax)
% h1 = axes;
set(gca, 'Ydir', 'reverse')
xlabel('Mass (g)')
ylabel('Length (cm)')
zlabel('Max Temperature (C)')
title('ChipSat Aerothermal Heating Parametric Analysis')

Bcoeffs = zeros(n*n,1); Tmaxs = zeros(n*n,1); TimesToGround = zeros(n*n,1);
for i=1:n
    for j=1:n
        Bcoeffs(n*(i-1)+j)= B_coeff(i,j);
        Tmaxs(n*(i-1)+j)= Tmax(i,j);
        TimesToGround(n*(i-1)+j)= TimeToGround(i,j);
    end
end
Bcoeffs = sort(Bcoeffs); Tmaxs = sort(Tmaxs); TimesToGround = sort(TimesToGround);
powerfit_Tmax = fit(Bcoeffs,Tmaxs,'power2')
quadraticfit_Time = fit(Bcoeffs,TimesToGround/24/3600,'poly2')

figure;
plot(powerfit_Tmax,Bcoeffs,Tmaxs)
xlabel('Ballistic Coefficient (kg/m^2)')
ylabel('Max Temperature (C)')
title('Peak Temperature vs. Ballistic Coefficient')

figure;
plot(quadraticfit_Time,Bcoeffs,TimesToGround/24/3600)
xlabel('Ballistic Coefficient (kg/m^2)')
ylabel('Deorbit Time (days)')
title('Deorbit time (from 400km altitude) vs. Ballistic Coefficient')

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


function [wdot] = rhs(t,w,p) %[zdot, u] if you want to export additional variables

    %unpack parameters
    names = fieldnames(p);
    for i=1:length(names)
      eval([names{i} '= p.' names{i} ';' ]); % e.g.  a = p.a;
    end
    
    %unpack w
    x = w(1); y = w(2); z = w(3); xdot = w(4); ydot = w(5); zdot = w(6);
    Q_aero = w(7); Q_rad = w(8); T = w(9); 
    
    r = [x; y; z];  %km
    rdot = [xdot; ydot; zdot]; %km/s
    
    Vs = rdot - cross(omega,r); %km/s 
    
    %A_eff = Az*Vs/abs(Vs);  %km^2
    
%     rho= interp1(altitudes,densities,norm(r)-R_E);   %kg/km^3

   alt = norm(r)-p.R_E;    %km

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
    
%     if norm(rdot)>0.35
    xddot = -mu_e*x/((x^2 + y^2 + z^2)^(3/2))  + (5*J2*x*(-x^2-y^2+2*z^2))/(2*(x^2+y^2+z^2)^(7/2)) + J2*x/(x^2+y^2+z^2)^(5/2) + drag(1);
    yddot = -mu_e*y/((x^2 + y^2 + z^2)^(3/2))  + (5*J2*y*(-x^2-y^2+2*z^2))/(2*(x^2+y^2+z^2)^(7/2)) + J2*y/(x^2+y^2+z^2)^(5/2) + drag(2);
    zddot = -mu_e*z/((x^2 + y^2 + z^2)^(3/2))  + (5*J2*z*(-x^2-y^2+2*z^2))/(2*(x^2+y^2+z^2)^(7/2)) - 2*J2*z/(x^2+y^2+z^2)^(5/2) + drag(3);
     

    wdot = [xdot; ydot; zdot; xddot; yddot; zddot; Qdot_aero; Qdot_rad; Tdot];
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