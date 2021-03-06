%modify original file "pendulum.m" to simulate SHM. Use it to calculate: 
%relative change in the total energy during one cycle and is the function
%"delta" uniformly small during the cycle? choose theta0 = 0.25 and its
%derivative = 0 and omega0^2=9
function [period,sol]= Pendulummodified(R,theta0,thetad0,grph) 
%kenergy,penergy,totenergy]

%given values omega0 = 3, theta0 = 0.25 and dertheta = 0
% Finds the period of a nonlinear pendulum given the length of the pendulum
% arm and initial conditions. All angles in radians.

%Setting initial conditions
if nargin==0 %Number of function input arguments
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

    %kenergy = 1/2*1.*w(:,2).*w(:,2); %setting m = 1
    %penergy = 1*g*R*(1-cos(w(:,1)));
    %totenergy = kenergy + penergy; 
    %delta_n = (totenergy(:) - totenergy(1)) / totenergy(1); %Total energy E(:) reshapes all elements of E into a single column vector. 
    
    %This has no effect if E is already a column vector.
    %MATLAB starts indexing at 1 unlike python

%if nargin==4
    %grph=0;
%end
m = 1;
g = 9.81;
omega = sqrt(g/R); 
%R = (g/omega)^2;

T= 2*pi/omega;
% number of oscillations to graph
N = 10;


tspan = [0 N*T]; %tspan = [t0 tf], integrates the system of differential equations w'=f(t,w) from 0 to N*T with initial conditions r0.
%opts = odeset('events',@events,'refine',6); %Here for future event finder
opts = odeset('refine',6);
r0 = [theta0 thetad0];
[t,w] = ode45(@proj,tspan,r0,opts,g,R); %Solve nonstiff differential equations ? medium order method

sol = [t,w]; %sol = ode45(___) returns a structure that you can use with deval to evaluate the solution at any point on the interval [0 N*T].
ind= find(w(:,2).*circshift(w(:,2), [-1 0]) <= 0); %take the entire 2nd column and shift to the left since [-1 -1] means shift ypward to the left 
%and [1 1] is shift downward to the right 
%the <= 0 is for the find function: true?

ind = chop(ind,4);
period= 2*mean(diff(t(ind)));

kenergy = 0.5 * m *(R*w(ind(1):ind(3),2)).^2;
penergy = m * g * R *(1-cos(w(ind(1):ind(3),1)));
totenergy = kenergy + penergy;
energy0 = 0.5*(R*thetad0).^2 + g*R*(1-cos(theta0));
denergy = (totenergy(:) - energy0)./energy0;
delta = atan(theta0/(omega*thetad0)); %atan: Inverse tangent in radians
y = theta0*sin(omega*t+delta); 


% Small-angle approximation solution

%delta = atan(theta0/(omega0*thetad0)); 
%y = theta0*sin(omega0*t+delta);
    
    
     
     
    %Total energy E(:) reshapes all elements of E into a single column vector.  
    %This has no effect if E is already a column vector.
    %MATLAB starts indexing at 1 unlike python
   
    
if grph
    figure
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
    
    figure
    subplot(2,1,1)
    plot(t(ind(1):ind(3)),denergy,'c*')
    title('Total Energy during one cycle')
    xlabel('t')
    ylabel('\Deltan')
    
    figure
    subplot(2,1,1)
    plot(t,w(:,1),'m')
    title('Position VS Time')
    xlabel('t')
    ylabel('\theta')
    
    figure
    subplot(2,1,1)
    plot(t,w(:,2),'m')
    title('Velocity VS Time')
    xlabel('t')
    ylabel('d\theta /dt')
    
    figure 
    subplot(2,1,1)
    plot(w(:,1),w(:,2),'g')
    title('Phase Space')
    xlabel('theta')
    ylabel('d\theta /dt')
end

end

%-------------------------------------------
%
function rdot = proj(t,r,g,R)
    rdot = [r(2); -g/R*sin(r(1))];
end
% if you add in the command window: the code ALWAYS returns the first line
%gives back one answer
% only but if we do: [T,Solution] = pendulum(0.4)
%you need to give it something to assign to 