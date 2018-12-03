%Pendulum_shm
%modify original file "pendulum.m" to simulate SHM:


%relative change in the total energy during one cycle and is the function
%"delta" uniformly small during the cycle? choose theta0 = 0.25 and its
%derivative = 0 and omega0^2=9


function [period,sol,kenergy,penergy,totenergy] = Pendulum_shm(omega0,theta0,thetad0,grph) 
% Finds the period of a SHM given the angular frequency
% arm and initial conditions. All angles in radians.

%Setting initial conditions
if nargin==0
    error('Must input length and initial conditions')
end
if nargin==1
   theta0 = pi/2;
   thetad0=0;
   grph=0;
end
if nargin==2
    thetad0 = 0;
    grph=1;
end
if nargin==3
    grph=1;
end
m=1;
g=9.81;
R=g/omega0^2; %omega0 is an input here
%omega = sqrt(g/R);
T= 2*pi/omega0;
% number of oscillations to graph
N = 10; %number of oscillations


tspan = [0 N*T];
%opts = odeset('events',@events,'refine',6); %Here for future event finder
opts = odeset('refine',6);
r0 = [theta0 thetad0];
[t,w] = ode45(@proj,tspan,r0,opts,g,R,m);
sol = [t,w];
ind= find(w(:,2).*circshift(w(:,2), [-1 0]) <= 0);
ind = chop(ind,4);
period= 2*mean(diff(t(ind)));


% Small-angle approximation solution
delta = atan(theta0/(omega0*thetad0));
y = theta0*sin(omega0*t+delta);

grph=0;
if grph % Plot Solutions of exact and small angle
    subplot(2,1,1)
    plot(t,w(:,1),'k-',t,y,'b--')
    legend('Exact','Small Angle')
    title('Exact vs Approximate Solutions')
    xlabel('t')
    ylabel('\phi')
    subplot(2,1,2)
    plot(t,w(:,1)-y,'k-')
    title('Difference between Exact and Approximate')
    xlabel('t')
    ylabel( '\Delta\phi')
end


time=floor(5*T); %floor (t): rounds each element of the duration array t to the nearest number of seconds less than or equal to that element.
index=find(t<=time);
index=max(index);
sol=sol(1:index,1:3);

kenergy = 0.5*m.*sol(:,3).*sol(:,3); %sol(:,3) is the velocity: v^2
penergy = 0.5*m*omega0^2.*sol(:,2).*sol(:,2); %sol(:,2) is the position: x^2
totenergy = kenergy + penergy; 
%delta_n = (totenergy(:) - totenergy(1)) / totenergy(1); %Total energy E(:) reshapes all elements of E into a single column vector. 
    
    %This has no effect if E is already a column vector.
    %MATLAB starts indexing at 1 unlike python


end
%-------------------------------------------
%
function rdot = proj(t,r,g,R,m)
    rdot = [r(2); -g/R*r(1)];
end