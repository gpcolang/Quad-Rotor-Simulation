%% HW4 Gabriel Colangelo Quadrotor Simulation
clear
close
clc

data = importdata('quadData.mat');
time = data(:,1); % seconds
global t0 tf
t0 = time(1,1);
tf = time(end);

% Quad Rotor Data
x = data(:,2); % Actual position in meters
y = data(:,3);
z = data(:,4);
u = data(:,5); % meters per second
v = data(:,6); 
w = data(:,7);
phi = data(:,8); % Actual euler angle in radians
theta = data(:,9);
psi = data(:,10);
p = data(:,11); % radians per second
q = data(:,12);
r = data(:,13);

% Part A
input = [u v w p q r]; 
init = [x(1,1) y(1,1) z(1,1) phi(1,1) theta(1,1) psi(1,1)]'; % Initial conditions
options = odeset('AbsTol',1e-8,'RelTol',1e-8);
[T,XX] = ode45(@diffeq,time,init,options,input); % Linear Simulation Part A

figure
plot(time,phi,'.r',time,XX(:,4),'-.k',time,theta,'.b',time,XX(:,5),'-.y',time,psi,'.c',time,XX(:,6),'-.m')
xlim([t0 tf])
legend('Actual \phi','Simulated \phi','Actual \theta','Simulated \theta','Actual \psi','Simulated \psi')
title('Part A Vehicle Attitude in Euler Angles vs. Time')
xlabel('Time [s]')
ylabel('Attitude [rad]')

figure
plot(time,x,'.r',time,XX(:,1),'-.k',time,y,'.b',time,XX(:,2),'-.y',time,z,'.c',time,XX(:,3),'-.m')
xlim([t0 tf])
legend('Actual X','Simulated X','Actual Y','Simulated Y','Actual Z','Simulated Z','Location','SouthWest')
title('Part A Position of the Center of Mass')
xlabel('Time [s]')
ylabel('Position [m]')

% Part B

% True Quaternion Values
q0 = (cos(psi/2).*cos(theta/2).*cos(phi/2))+(sin(psi/2).*sin(theta/2).*sin(phi/2));
q1 = (cos(psi/2).*cos(theta/2).*sin(phi/2))-(sin(psi/2).*sin(theta/2).*cos(phi/2));
q2 = (cos(psi/2).*sin(theta/2).*cos(phi/2))+(sin(psi/2).*cos(theta/2).*sin(phi/2));
q3 = (sin(psi/2).*cos(theta/2).*cos(phi/2))-(cos(psi/2).*sin(theta/2).*sin(phi/2));

input1 = [p q r u v w]; 
init1 =  [q0(1,1) q1(1,1) q2(1,1) q3(1,1) x(1,1) y(1,1) z(1,1)];
options = odeset('AbsTol',1e-8,'RelTol',1e-8);
[T1,X1] = ode45(@diffeq1,time,init1,options,input1); % Simulation part B

figure
plot(time,q0,'.r',time,X1(:,1),'-.k',time,q1,'.b',time,X1(:,2),'-.y',time,q2,'.c',time,X1(:,3),'-.m',time,q3,'.g',time,X1(:,4),'-.r')
xlim([t0 tf])
legend('Actual q0','Simulated q0','Actual q1','Simulated q1','Actual q2','Simulated q2','Actual q3','Simulated q3')
title('Part B Quaternion Component Values vs Time')
xlabel('Time [s]')
ylabel('Quaternion Component Value')

figure
plot(time,x,'.r',time,X1(:,5),'-.k',time,y,'.b',time,X1(:,6),'-.y',time,z,'.c',time,X1(:,7),'-.m')
xlim([t0 tf])
legend('Actual X','Simulated X','Actual Y','Simulated Y','Actual Z','Simulated Z','Location','SouthWest')
title('Part B Position of the Center of Mass')
xlabel('Time [s]')
ylabel('Position [m]')

function xdot = diffeq1(t,x,input1)

global t0 tf

q0 = x(1,1); % q0
q1 = x(2,1); % q1
q2 = x(3,1); % q2
q3 = x(4,1); % q3
X = x(5,1); % x
y = x(6,1); % y
z = x(7,1); % z

input1 = interp1(linspace(t0,tf,length(input1)),input1,t)';

p = input1(1,1); % p
q = input1(2,1); % q
r = input1(3,1); % r
u = input1(4,1); % u
v = input1(5,1); % v
w = input1(6,1); % w

% Quaternion Rotational Rate kinematics
xdot(1,1) = ((-p*q1) + (-q*q2) + (-r*q3))/2;
xdot(2,1) = ((p*q0) + (r*q2) + (-q*q3))/2;
xdot(3,1) = ((q*q0)+ (-r*q1) + (p*q3))/2;
xdot(4,1) = ((r*q0) + (q*q1) + (-p*q2))/2;

% Tranlational Kinematics
xdot(5,1) = u*(q0^2 + q1^2 - q2^2 - q3^2) - v*(2*q0*q3 - 2*q1*q2) + w*(2*q0*q2 + 2*q1*q3);
xdot(6,1) = v*(q0^2 - q1^2 + q2^2 - q3^2) + u*(2*q0*q3 + 2*q1*q2) - w*(2*q0*q1 - 2*q2*q3);
xdot(7,1) = w*(q0^2 - q1^2 - q2^2 + q3^2) - u*(2*q0*q2 - 2*q1*q3) + v*(2*q0*q1 + 2*q2*q3);
end

function xdot = diffeq(t,x,input)

global t0 tf

s = @(y)sin(y);
c = @(y)cos(y);
T = @(y)tan(y);

X = x(1,1); % x
y = x(2,1); % y
z = x(3,1); % z
phi = x(4,1); % phi
theta = x(5,1); % theta
psi = x(6,1); % psi

input = interp1(linspace(t0,tf,length(input)),input,t)';

u = input(1,1); % u
v = input(2,1); % v
w = input(3,1); % w
p = input(4,1); % p
q = input(5,1); % q
r = input(6,1); % r

% Navigation Equations and 321 Euler Rate Kinematics
xdot(1,1) = (c(theta)*c(psi)*u) + (s(phi)*s(theta)*c(psi)-c(phi)*s(psi))*v + (c(phi)*s(theta)*c(psi)+s(phi)*s(psi))*w;
xdot(2,1) = (c(theta)*s(psi)*u) + (s(phi)*s(theta)*s(psi)+c(phi)*c(psi))*v + (c(phi)*s(theta)*s(psi)-s(phi)*c(psi))*w;
xdot(3,1) = -s(theta)*u + (s(phi)*c(theta))*v + (c(phi)*c(theta))*w;
xdot(4,1) = p + q*s(phi)*T(theta) + r*c(phi)*T(theta);
xdot(5,1) = q*c(phi) - r*s(phi);
xdot(6,1) = (q*s(phi) + r*c(phi))/(c(theta));
end



