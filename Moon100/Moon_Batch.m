% Moon Batch 
% using Hunter's orbital dynamics 
% last editted: 12/16/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% parameters
p.omega = [0, 0, 2.6638e-06]';  % angular velocity of moon, rad/s, locked with Saturn
p.R_M = 1737;                 % radius of moon, km
G = 6.674*10^(-11)/(10^9);       % gravitational constant, km^3/(kg s^2)
M = 7.34767e22;            % mass of moon, kg
p.mu_M = G*M;         % moon gravitational parameter, km^3/s^2
p.J2 = 0; %(2.7221e-5)*(p.R_M^2)*(p.mu_M); % J2, km^5/s^2
initial_alt = 100;   % km
v_orbital = sqrt(G*M/(p.R_M+initial_alt));    %km/s, initial orbital velocity of chipsat mothership
inclination = deg2rad(0);  % orbital inclination, radians

% initial conditions
totaltime = linspace(0, 1*24*3600, 1*24*3600);     % total time

accuracy = 1e-6;
options = odeset('AbsTol', accuracy, 'RelTol', accuracy,'Events', @stopEventsFxn);

n = 100;      %number of simulations,100

%orbital transfer calc
rp = p.R_M; %periapsis (km), we want elliptical orbit to meet lunar surface
ra = p.R_M+initial_alt; %apoapsis (km), current orbital height
a = (rp+ra)/2; %semi major axis for vis viva eqn
r = ra; %current position for vis viva eqn
v_desired = sqrt(p.mu_M*(2/r-1/a)); %km/s, velocity at apoapsis needed to reach surface at periapsis
v_desired = v_desired-5/1000; %km/s extra margin to account for additional dispersion
vd1 = v_desired - v_orbital; %deployment velocity, km/s negative value for retrograde deploy that slows velocity to desired value
vd2 = 5/1000; %km/s, individual sprite kick velocity, radially deployed

figure;       % create figure for batch plotting
r = zeros(length(totaltime),n);
x = zeros(length(totaltime),n); y = zeros(length(totaltime),n); z = zeros(length(totaltime),n);
xdot = zeros(length(totaltime),n); ydot = zeros(length(totaltime),n); zdot = zeros(length(totaltime),n);
xpos = zeros(length(totaltime),n); ypos = zeros(length(totaltime),n);
time = zeros(length(totaltime),n);
mass = zeros(length(n),1); area = zeros(length(n),1);
TimeToGround = zeros(length(n),1); V_impact = zeros(length(n),1);

lat = zeros(1,n); long = zeros(1,n);

% varying simulation parameters
for i = 1:n
    i           % tells us what simulation number we are on
    q = 2*pi/n*i;
    mass(i) = normrnd(3,0.1)/1000;  % mass of spacecraft, kg 5.8250e-06 for justin's 1cmx1cm chipsat
    p.m = mass(i);
    
    initial_position = [0,initial_alt, 0];     % km, randomizing height of chipsat in deployer
    initial_position = (initial_position/norm(initial_position))*(p.R_M+norm(initial_position));
    %initial_velocity = [(v_desired+vd2*sin(q))*cos(inclination), vd2*cos(q), (v_desired+vd2*sin(q))*sin(inclination)];   %km/s
    initial_velocity = [(v_desired)*cos(inclination),vd2*cos(q),(v_desired)*sin(inclination)+vd2*sin(q)];   %km/s

    w0 = [initial_position'; initial_velocity'];

    f = @(t,w)  rhs(t,w,p);
    [tarray zarray] = ode45(f, totaltime, w0, options);
    %[~,U] = cellfun(@(t,w)  rhs(t,w.',p), num2cell(tarray), num2cell(zarray,2),'uni',0);
    %U = cell2mat(U);
    
    for j = 1:length(tarray)
        x(j,i) = zarray(j,1);
        y(j,i) = zarray(j,2);
        z(j,i) = zarray(j,3);
        xdot(j,i) = zarray(j,4);
        ydot(j,i) = zarray(j,5);
        zdot(j,i) = zarray(j,6);
        r(j,i) = sqrt(x(j,i)^2+y(j,i)^2+z(j,i)^2)-p.R_M;
        time(j,i) = tarray(j);
        
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
    
    TimeToGround(i) = tarray(end); %s
    %V_impact(i) = U(length(tarray),2);
   
    temp = zeros(length(tarray),1);
    for k = 1:length(tarray)
        temp(k) = 5.5;
    end
    yyaxis left
    plot(time(1:length(tarray),i)/(60),r(1:length(tarray),i))
    ylabel('Altitude (km)')
    xlabel('Time (min)')
    hold on
    yyaxis right
    plot(time(1:length(tarray),i)/(60),temp(1:length(tarray)),'--')
    ylabel('Temperature (C)')
    legend('Altitude','Temperature','Location','NorthEast')
    %title('ChipSat Altitude vs. Time')
 
end

%figure;
%plot(tarray/(24*3600),U)


figure;
worldmap('World')
%axesm('MapProjection','hammer','Grid','on','Frame','on')
scatterm(lat,long,'r','filled')
%title('ChipSat Landing Locations')

figure;
%axis equal
hold on
for i=1:n
    plot(xpos(:,i),ypos(:,i));
end
%title('ChipSat Deorbit Trajectories')
xlabel('MCI X - Altitude (km)')
ylabel('MCI Y - Altitude (km)')

%error elipse
%[r_ellipse, X0, Y0] = error_ellipse((lat-mean(lat))*p.R_M,(long-mean(long))*p.R_M);


%scatter plot of landing locations
figure;
scatter(deg2rad(long-mean(long))*p.R_M,deg2rad(lat-mean(lat))*p.R_M)
%hold on
%plot((r_ellipse(:,1) + X0), (r_ellipse(:,2) + Y0), 'r-');
xlabel('x (km)'); ylabel('y (km)');
%axis equal

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
  
    r = [x; y; z];  %km
    rdot = [xdot; ydot; zdot]; %km/s
    
    Vs = rdot - cross(omega,r); %km/s 
    
    %A_eff = Az*Vs/abs(Vs);  %km^2
    
%     rho= interp1(altitudes,densities,norm(r)-R_E);   %kg/km^3

    alt = norm(r)-p.R_M    %km

    xddot = -mu_M*x/((x^2 + y^2 + z^2)^(3/2))  + (5*J2*x*(-x^2-y^2+2*z^2))/(2*(x^2+y^2+z^2)^(7/2)) + J2*x/(x^2+y^2+z^2)^(5/2);
    yddot = -mu_M*y/((x^2 + y^2 + z^2)^(3/2))  + (5*J2*y*(-x^2-y^2+2*z^2))/(2*(x^2+y^2+z^2)^(7/2)) + J2*y/(x^2+y^2+z^2)^(5/2);
    zddot = -mu_M*z/((x^2 + y^2 + z^2)^(3/2))  + (5*J2*z*(-x^2-y^2+2*z^2))/(2*(x^2+y^2+z^2)^(7/2)) - 2*J2*z/(x^2+y^2+z^2)^(5/2);

    
    wdot = [xdot; ydot; zdot; xddot; yddot; zddot];
    %u = [alt,norm(Vs)];
end 

function [position,isterminal,direction] = stopEventsFxn(t,w)
%Stops integration after reaching altitude of 0
x = w(1); y = w(2); z = w(3);
r = [x; y; z];  %km
altitude = norm(r)-1737;        % hard coded Titan radius 
position = altitude<0; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end