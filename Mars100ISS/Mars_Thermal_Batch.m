% Mars Thermal Batch 
% using Hunter's orbital dynamics
% Mars densities
% with Justin's thermal equations 
% last editted: 12/16/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% import Mars data
% make sure excel file is in same folder as this matlab file 
marsdata = xlsread('Mars_data_all.xlsx');
p.mars_altitude = marsdata(1:end,2);       % km
p.mars_densities = marsdata(1:end,4);    % kg/m^3 
p.mars_temps = marsdata(1:end,5); % Kelvin
p.mars_pressures = marsdata(1:end,6); %N/m^2
%MARSGRAM data. tpdmsy11.txt & tpdloy11.txt Ls=180, Lat =7.5


% parameters
h = 50/1000/1000;                     % spacecraft height, km, 35
w = 50/1000/1000;                     % spacecraft width, km
%d = 2/1000/1000;                      % spacecraft depth, km
%p.Az = [h*d, h*d, h*h]*(10^-12);      % z-direction of coil surface, km^2
%p.A = h*h;                    % surface area, km^2
p.l = h;
p.Cd_FM = 2.67;        % drag coefficient 2.67 face on, 1.43 averaged, in FM
p.Cd_SS = 1.28;  % drag coefficient 1.28 face on, 0.8155 averaged, in subsonic
%p.m = 3/1000;                 % mass of spacecraft, kg 
p.omega = [0, 0, 7.0902e-05]';  % angular velocity of Mars, rad/s
p.R_M = 3389.5;                 % radius of Mars, km
G = 6.674*10^(-11)/(10^9);       % gravitational constant, km^3/(kg s^2)
M = 6.4171e23;            % mass of Mars, kg
p.mu_m = G*M;         % Mars gravitational parameter, km^3/s^2
p.J2 = (1960.45e-6)*(p.R_M^2)*(p.mu_m);             % J2, km^5/s^2
initial_alt = 200;   % km, 374
v_orbital = sqrt(p.mu_m/(p.R_M+initial_alt));    %km/s, initial orbital velocity
inclination = deg2rad(50);  % around ISS orbital inclination, radians

% Justin's parameters
p.N_A = 6.022*10^23; %Avogadro's number
p.Cp = 1090; %J/(kg K) specific heat 678 silicon wafer, 1090 kapton, about same for FR4
p.epsilon = 0.85; %emissivity
p.sbsigma = 5.67e-8; %stephen-boltzman constant
p.diam = 0.33e-9; %CO2 molecule diameter

%ratio of specific heat
p.gammaCO2 = 1.28; %carbon dioxide
p.R = 8.314; %J/(mol*K) Universal gas constant
p.kb = 1.380649e-23; %boltzmann constant J/K

%viscosity model constants
%get temperature vs. altitude profile
%use alt in ode15s to find temp, use temp to find dyn visc.
p.mu0_CO2 = 1.37e-5; %reference viscosity in centipoise at reference temperature To
p.S_CO2 = 222; %Sutherland's constant
p.T0_CO2 = 273; %reference temperature in degreees Kelvin

% Justin's initial conditions
Q_aero0 = 0;
Q_rad0 = 0;
T0 = 250;

%Re2 parameter calculation
gammaCO2 = p.gammaCO2;
mm = 28.97/1000; %/N_A; %kg/mol molar mass of air
a_inf = sqrt(gammaCO2*p.R*210/mm); %m/s
M1 = 1000*v_orbital/a_inf;
Mp2 = sqrt(2/(gammaCO2-1))*(M1^2 - 1)/(sqrt(1+(gammaCO2-1)/2*M1^2)*sqrt(2*gammaCO2/(gammaCO2-1)*(M1^2)-1));
p.Re2_param = Mp2*sqrt(pi*gammaCO2/2); 
%this gets used in the ODE function as Re2 = Re2_param/Kn2


% % Justin's viscosity model
% height = [-1000 0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 15000 ...
%     20000 25000 30000 40000 50000 60000 70000 80000 200000]';
% dyn_visc = [1.821 1.789 1.758 1.726 1.694 1.661 1.628 1.595 1.561 1.527 ...
%     1.493 1.458 1.422 1.422 1.448 1.475 1.601 1.704 1.584 1.438 1.321 0]';
% 
% for i=1:length(dyn_visc)
%     dyn_visc(i)=dyn_visc(i)*10^-5;
% end
% 
% p.visc = polyfit(height,dyn_visc,3);

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
totaltime = linspace(0, 5*24*3600, 5*24*3600);     % total time

accuracy = 1e-8;
options = odeset('AbsTol', accuracy, 'RelTol', accuracy,'Events', @stopEventsFxn);

n = 100;      %number of simulations,100
vd = 1/1000;     % deployment velocity, km/s
figure;       % create figure for batch plotting
r = zeros(length(totaltime),n);
x = zeros(length(totaltime),n); y = zeros(length(totaltime),n); z = zeros(length(totaltime),n);
xpos = zeros(length(totaltime),n); ypos = zeros(length(totaltime),n);
%xdot = zeros(length(totaltime),n); ydot = zeros(length(totaltime),n); zdot = zeros(length(totaltime),n);
%Q_aero = zeros(length(totaltime),n); Q_rad = zeros(length(totaltime),n); T = zeros(length(totaltime),n);
%dvdot = zeros(length(totaltime),n); time = zeros(length(totaltime),n);
mass = zeros(length(n),1);
area = zeros(length(n),1);
Tmax = zeros(length(n),1);
TimeToGround = zeros(length(n),1);
lat = zeros(1,n); long = zeros(1,n);

% varying simulation parameters
for i = 1:n
    i           % tells us what simulation number we are on
    area(i) = h*w*normrnd(1,0.01);  % km^2, tumbling average area (1 = worse case, face on)
    %area(i) = h*w;
    p.Area_surf = 2*area(i);
    p.Area = area(i);
    q = 2*pi/n*i;
    mass(i) = normrnd(3/1000,0.1/1000);                 % mass of spacecraft, kg 
    %mass(i) = 3/1000; %kg, 5.8250e-06 for justin's 1cmx1cm chipsat
    p.m = mass(i);
    
    %initial_position = [initial_alt, 0, 0];
    initial_position = [normrnd(initial_alt,0.1/1000), 0, 0];     % km, randomizing height of chipsat in deployer
    initial_position = (initial_position/norm(initial_position))*(p.R_M+norm(initial_position));
    initial_velocity = [0, v_orbital*cos(inclination)+vd*cos(q), v_orbital*sin(inclination)+vd*sin(q)];   %km/s
    %initial_velocity = [0, v_orbital, 0];   %km/s
    initial_dv = 0;                            %km/s^2
    %initial_A = A/(1000^2);                      % initial area, m^2
    w0 = [initial_position'; initial_velocity'; initial_dv; Q_aero0; Q_rad0; T0];

    f = @(t,w)  rhs(t,w,p);
    [tarray zarray] = ode15s(f, totaltime, w0, options);
    %[~,U] = cellfun(@(t,w)  rhs(t,w.',p), num2cell(tarray), num2cell(zarray,2),'uni',0);
    %U = cell2mat(U);
    
    for j = 1:length(tarray)
        x(j,i) = zarray(j,1);
        y(j,i) = zarray(j,2);
        z(j,i) = zarray(j,3);
        %xdot(j,i) = zarray(j,4);
        %ydot(j,i) = zarray(j,5);
        %zdot(j,i) = zarray(j,6);
        %dvdot(j,i) = zarray(j,7);
        %Q_aero(j,i) = zarray(j,8); 
        %Q_rad(j,i) = zarray(j,9); 
        T(j,i) = zarray(j,10);
        r(j,i) = sqrt(x(j,i)^2+y(j,i)^2+z(j,i)^2)-p.R_M;
        %time(j,i) = tarray(j);
        
        if x(j,i)>=0
            xpos(j,i) = x(j,i)-(p.R_M*abs(x(j,i))/(r(j,i)+p.R_M));
        else
            xpos(j,i) = x(j,i)+(p.R_M*abs(x(j,i))/(r(j,i)+p.R_M));
        end
        if y(j,i)>=0
            ypos(j,i) = y(j,i)-(p.R_M*abs(y(j,i))/(r(j,i)+p.R_M));
        else
            ypos(j,i) = y(j,i)+(p.R_M*abs(y(j,i))/(r(j,i)+p.R_M));
        end
       
    end
    
    index = find(r(:,i)<1 & r(:,i)>0);     % find index when chip is at 0 km
    last = index(end);
    ECI = [x(last,i) y(last,i) z(last,i)]; % final 0 km altitude states
    ECIx = ECI(1); ECIy = ECI(2); ECIz = ECI(3);
    lat(i) = rad2deg(atan2(ECIz,sqrt(ECIx^2+ECIy^2)));  %deg
    long(i) = rad2deg(atan2(ECIy,ECIx));    %deg
    
    Tmax(i) = max(zarray(:,10))-273; %max temp, C
    TimeToGround(i) = length(tarray); %s
    
end

% plot on world map 
figure;
worldmap('World')
scatterm(lat,long,'r','filled')
%title('ChipSat Landing Locations')

figure;
axis equal
hold on
for i=1:n
    plot(xpos(:,i),ypos(:,i));
end
%title('ChipSat Deorbit Trajectories')
xlabel('ECI X - Altitude (km)')
ylabel('ECI Y - Altitude (km)')

figure;
histogram(Tmax)
xlabel('Peak Temperature (C)')
ylabel('Number of ChipSats')

figure;
histogram(TimeToGround/3600)
xlabel('Deorbit Time (min)')
ylabel('Number of ChipSats')


function [wdot] = rhs(t,w,p) %[zdot, u] if you want to export additional variables

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

    alt = norm(r)-p.R_M;    %km

    % density model 
   if alt>220
       alt_m=alt*1000;  % convert km-->m
       rho = exp(2.65472e-11*(alt_m^5) - 2.45558e-08*(alt_m^4) + 6.31410e-06*(alt_m^3) + 4.73359e-04*(alt_m^2) - 0.443712*alt_m + 23.79408);
   else 
        rho = interp1(mars_altitude,mars_densities,alt,'makima','extrap'); 
   end
   
    rho_m = rho;    % kg/m^3
    rho = rho*10^9; %conversion to kg/km^3
     
    %rho = 1;
  
    T_CO2 = interp1(mars_altitude,mars_temps,alt,'makima','extrap');
    P_inf = interp1(mars_altitude,mars_pressures,alt,'makima','extrap');
    mu = mu0_CO2*((T_CO2/T0_CO2)^3/2)*(T0_CO2+S_CO2)/(T_CO2+S_CO2);
      
    % Justin's thermal code
%      mu = polyval(visc,z); % dynamic viscosity (N s/m^2)
      L = l*1000; %characteristic length, m
      Re = rho_m*1000*norm(Vs)*L/mu;
      mm = 44.01/1000; %kg/mol molecular mass of CO2
      
      a_inf = sqrt(gammaCO2*R*T_CO2/mm); %m/s
      Ma = 1000*norm(Vs)/a_inf;
      %zeta_a = R*255/(sqrt(2)*pi*(diam^2)*N_A*P_inf); %mean free path
      zeta_a = mu/rho_m*sqrt(pi*mm/(N_A*2*kb*T_CO2));
      Kn = zeta_a/L;
      %Kn = sqrt(gammaAir*pi/2)*Ma/Re; %consistent with above, not below
      %Kn = kb*T_inf/(sqrt(2)*diam^2*P_inf*L);
      
      
      M1 = Ma;
      rho0_rho1 = (1+(gammaCO2-1)/2*M1^2)^(1/(gammaCO2-1));
      rho1_rho2 = ((gammaCO2-1)*(M1^2)+2)/(gammaCO2+1)/M1^2;
      Kn2 = Kn*rho1_rho2*rho0_rho1;
      Re2 = Re2_param/Kn2;      

      ST_FM = 1;
      Cs = 1/sqrt(2);
      ST_C = Cs*2.1/sqrt(Re2);
      
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
        Cd = Cd_SS/sqrt(1+(Cd_SS/Cd_FM)^2);
        if Cd < Cd_SS
            Cd=Cd_SS;
        end        
        %Cskin = Cc;
    end    
      
     %ST=0.1;
    Qdot_aero = 1/2*ST*rho_m*(1000*norm(Vs))^3*(Area*1000^2); %might have to change sdot to velocity resultant vector magnitude
    Qdot_rad = sbsigma*epsilon*(T^4-210^4)*(Area_surf*1000^2); 
    Tdot = (Qdot_aero - Qdot_rad)/(m*Cp);

    drag =(1/2)*Cd*(Area)*rho*((Vs')*Vs)/m.*(-Vs)/norm(Vs);  % deceleration due to drag, km/s^2
    % Hunter's code did not include Cd
    
    dvdot = norm(drag); %km/s^2
    
    xddot = -mu_m*x/((x^2 + y^2 + z^2)^(3/2))  + (5*J2*x*(-x^2-y^2+2*z^2))/(2*(x^2+y^2+z^2)^(7/2)) + J2*x/(x^2+y^2+z^2)^(5/2) + drag(1);
    yddot = -mu_m*y/((x^2 + y^2 + z^2)^(3/2))  + (5*J2*y*(-x^2-y^2+2*z^2))/(2*(x^2+y^2+z^2)^(7/2)) + J2*y/(x^2+y^2+z^2)^(5/2) + drag(2);
    zddot = -mu_m*z/((x^2 + y^2 + z^2)^(3/2))  + (5*J2*z*(-x^2-y^2+2*z^2))/(2*(x^2+y^2+z^2)^(7/2)) - 2*J2*z/(x^2+y^2+z^2)^(5/2) + drag(3);
  
    
    wdot = [xdot; ydot; zdot; xddot; yddot; zddot; dvdot; Qdot_aero; Qdot_rad; Tdot];
end 

function [position,isterminal,direction] = stopEventsFxn(t,w)
%Stops integration after reaching altitude of 0
x = w(1); y = w(2); z = w(3);
r = [x; y; z];  %km
altitude = norm(r)-3389.5;        % hard coded Mars radius 
position = altitude<0; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end