clc;clear; close all;
%motor parameters
Vdc = 120;              %DC voltage 120V
a_max = 1;              %max current 1A
theta_dot_max = 0.5;    %max angular speed 0.5 deg/s
time_max = 30;          %time to reach max angular speed
t_sample = 0.5;

%derived parameters
tau = 0.6321*time_max;  %time constant: reach 63% of max speed
rpm = theta_dot_max/6;  %rpm of max speed

%plant estimation - this is the best I can do without knowing 
%things like inductance, resistance, etc
%http://ctms.engin.umich.edu/CTMS/index.php?aux=Activities_DCmotorB
%Plant = rpm/(tau*s + 1)
[A,B,C,D] = tf2ss([0 rpm*6], [tau 1]);%convert to continuous ss
plant_ss = ss(A,B,C,D);                    %create sysc object for ss
plant_ssd = c2d(plant_ss,t_sample,'zoh')            %discretize

%augment matrices to get states for alpha, beta angles
%x vector is [theta_a; theta_B; theta_a_dot; theta_B_dot]
%u vector is [u_a; u_B]
A2 = [1 0 t_sample 0; 0 1 0 t_sample;...
      0 0 plant_ssd.A 0; 0 0 0 plant_ssd.A];        
B2 = [0 0; 0 0; plant_ssd.B 0; 0 plant_ssd.b];
C2 = [1 0 0 0;0 1 0 0];
D2=[0 0;0 0];

Plant_ssd_augmented = ss(A2,B2,C2,D)

%In controllable canonical form

A1=[1 t_sample;0 plant_ssd.A]; B1=[0;plant_ssd.B]; C1=[1 0]; D1=[0];
[n,d]=ss2tf(A1,B1,C1,D1)
p=poly(A1)
roots(p)
p1=poly([0.8257+0.1865i,0.8257-0.1865i])
%1.0000   -1.6514    0.7166
%CEd=A1*A1-1.778644*A1+0.818731*eye(2);
CEd=A1*A1*1.6514*A1+0.7166*eye(2);
%2x2
%K=[0 1]*inv(ctrb(A1,B1))*CEd
K=[0.2641    0.6537]
A1dash=(A1-B1*K)
[n1dash,d1dash]=ss2tf(A1dash,B1,C1,D1)
%K0=0.1625
K0=0.2643
B1dash=B1*K0
closedLoop=ss(A1dash,B1dash,C1,D1,t_sample)
figure(1)
step(closedLoop)

controlLoop1=ss(A1dash,B1dash,-K,[K0],t_sample)
figure(2)
step(controlLoop1)


%4x4
K2=[K(1) 0 K(2) 0;0 K(1) 0 K(2)]
A2dash=(A2-B2*K2)

B2dash=B2*K0
closedLoop2=ss(A2dash,B2dash,C2,D2,t_sample)
figure(3)
step(closedLoop2)

controlLoop2=ss(A2dash,B2dash,-K2,[K0 0;0 K0],t_sample)
figure(4)
step(controlLoop2)




