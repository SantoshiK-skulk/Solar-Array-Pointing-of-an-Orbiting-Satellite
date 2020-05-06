% Run through / animate horizons data 
% Hamil Shah

clc; 
clear; 
close all;

%% Read Data

% Format
% Julian Day,   Calendar Day (+ Time),  X,  Y,  Z,  Vx,     Vy,     Vz

% Coordinate System: Ecliptic and Mean Equinox of Reference Epoch 
% Data in 1 minute increments
% Positions:    km
% Velocity:     km/s

fid = fopen('iss_data_2hz.txt');
iss = textscan(fid,'%f,%[^,],%f,%f,%f,%f,%f,%.15f,',14400);
fclose(fid);

fid = fopen('sun_data_2hz.txt');
sun = textscan(fid,'%f,%[^,],%f,%f,%f,%f,%f,%.15f,',14400);
fclose(fid);

p_iss = [iss{3},iss{4},iss{5}];
v_iss = [iss{6},iss{7},iss{8}];

p_sun = [sun{3},sun{4},sun{5}];
v_sun = [sun{6},sun{7},sun{8}];
%% Animate ISS 2 hr. orbit data (minute by minute)
clearvars all -except p_iss v_iss p_sun v_sun
% 
% 
% figure(1); clf
% hold on; 
% grid on;
% title('Position of ISS relative to Earth for 2 hours');
% xlabel('x'); ylabel('y'); zlabel('z');
% 
% % Earth (assumed perfect model for visual plot)
Re = 6371;  % Radius of earth
% [x,y,z] = sphere(30);
% E = surf(x*Re, y*Re, z*Re);
% colormap summer
% shading interp
% 
% axis equal
% axis manual
% for i2 = 1:40:length(p_iss(:,1))
%     a = plot3(p_iss(i2,1), p_iss(i2,2), p_iss(i2,3),'O');
%     b = line(p_iss(1:i2,1), p_iss(1:i2,2), p_iss(1:i2,3));
%     axis equal
%     drawnow
%     
%     if i2 < length(p_iss(:,1))-40
%         delete(a);
%         delete(b);
% %         delete(c);
%     end
% end

%% Plot useful information?
iss_alt = sum(p_iss.^2,2).^.5-Re;
figure(2); clf;
plot(0:(length(iss_alt)-1),iss_alt)
title("ISS Altitude over 2 hrs (km)")
ylabel("Altitude (km)");
xlabel("Time (minutes)");


%% Calculate desired angles

% Convert ECI (Eart centered) frame ISS-Sun vector to LVLH frame
giss = p_iss;
viss = v_iss;

gsun = p_sun;
vsun = v_sun;

zliss = -giss;
yliss = cross(zliss,viss);
xliss = cross(yliss,zliss);

xliss_unit = xliss./(sum(xliss.^2,2.).^.5);
yliss_unit = yliss./(sum(yliss.^2,2).^.5);
zliss_unit = zliss./(sum(zliss.^2,2).^.5);

u_vec = p_sun - p_iss;

u_xform = zeros(size(u_vec));

X = [1 0 0]';
Y = [0 1 0]';
Z = [0 0 1]';
for i = 1:size(u_vec,1)
    u = u_vec(i,:)';
    x = xliss_unit(i,:)';
    y = yliss_unit(i,:)';
    z = zliss_unit(i,:)';

    Q = [...
    cos(getAngle(X,x)),cos(getAngle(X,y)),cos(getAngle(X,z));...
    cos(getAngle(Y,x)),cos(getAngle(Y,y)),cos(getAngle(Y,z));...
    cos(getAngle(Z,x)),cos(getAngle(Z,y)),cos(getAngle(Z,z))];

    u_prime = Q'*u;
    u_xform(i,:) = u_prime';
    
    % Calculated desired angle in controllable directions
    alpha_des(i,:) = atan2d(u_prime(3),u_prime(1)); % main / SARJ motors
    beta_des(i,:) = atan2d(u_prime(2),sqrt(sum(u_prime([1,3]).^2))); % small / individual motors
%     yaw_des(i,:) = atan2d(u_prime(2),u_prime(1));
end

% Orbit is approx 1150 indicies long
% Put eclipse at ~31.5 minutes in for 30 min
t_sample = 0.5;
e_start = 1 + 31.5*60*(1/t_sample);
e_end = e_start + 30*60*(1/t_sample);

idx = find(abs(diff(alpha_des)) > 100); idx = idx(1);
alpha_des((idx+1):end) = alpha_des((idx+1):end) - (alpha_des(idx+1)-alpha_des(idx))-(360/93/60*0.5);

% alpha_des(e_start:e_end) = alpha_des(e_start-1);
beta_des(e_start:e_end) = beta_des(e_start-1);

alpha_vel_des = ones(size(alpha_des)).*(-360/(93*60)); 
alpha_vel_des(end+1) = alpha_vel_des(end);
beta_vel_des = diff(beta_des)./0.5;
beta_vel_des(end+1) = beta_vel_des(end);

alpha_des = [(1:length(alpha_des))'.*0.5-0.5 alpha_des];
beta_des = [(1:length(beta_des))'.*0.5-0.5 beta_des];
save('angle_desired','alpha_des','alpha_vel_des','beta_des','beta_vel_des');
figure(3); clf;
plot(beta_des(:,2),'LineWidth',2)
title("Desired BGA Angle");
ylabel("Angle (deg.)");
xlabel("Index");

figure(4); clf;
plot(alpha_des(:,2),'LineWidth',2);
title("Desired SARJ angle");
ylabel("Angle (deg.)");
xlabel("Index");
