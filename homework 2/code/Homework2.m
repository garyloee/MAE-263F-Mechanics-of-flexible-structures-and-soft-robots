close all;
clear all;
clc

global Fg M dt
global kappaBar E I G J voronoiL
global A refL

%% Setup simulation
% number DOF
N=50; % number of nodes
ne=N-1; % number of edges
ndof=3*N+ne; %location of each nodes plus twisting of each section

% time frame
dt=0.01; % time step size
totalTime=5; % simulation time
Nsteps=round(totalTime/dt); % total amount of steps


%% Set up physical parameters
% grometry of the rod
l=0.2; % rodLength of 20cm
Rn=0.02; % Natural curvature bend radius
r0=0.001; % Cross section radius of 1mm

% Material Parameters
rho=1000; % density of the rod
E=10e6; % Young's Modulus of 10 MPA
nu=0.5; % poisson's ratio
G=E/(2*(1+nu)); % shear modulus

% stiffness notions
I=pi*r0^4/4; % moment of inertia I
J=pi*r0^4/2; % moment of inertia J
A=pi*r0^2; % cross section area A

% gravity-External acceleration
g=[0;0;-9.81];

%% set up each node
% mass matrix
totalMass=A*l*rho; % total mass in KG
dMass=totalMass/ne; % mass of each section

% asign mass to each dof
massVector=zeros(ndof,1); % initializing mass for all dof of each node
for i=1:N % setup mass of midle nodes
    if i==1 || i==N
        massVector([4*i-3;4*i-2;4*i-1])=dMass/2; % half the mass for first and last node
    else
        massVector([4*i-3;4*i-2;4*i-1])=dMass;
    end
end

for i=1:ne
    massVector(4*i)=1/2*dMass*r0^2; % rotaional Inertia of each section = 1/2 mr^2
end

M=diag(massVector); % preparing proper Mass matrix

%% Initialize nodes vector
nodes=zeros(N,3);
dTheta=l/Rn*(1/(N-1)); % delta theta
for i=1:N
    nodes(i,1)=Rn*cos((i-1)*dTheta); % initializing xk
    nodes(i,2)=Rn*sin((i-1)*dTheta);
    nodes(i,3)=0;
end

% initial condition
q0=zeros(ndof,1);
for i=1:N
    % putting vector into each dof, with forth index being theta =0
    q0([4*i-3;4*i-2;4*i-1])=nodes(i,:); 
end

% velocity = 0
u=zeros(ndof,1);

%% defining reference length of each section for stretching force
refL=zeros(ne,1);
for i=1:ne
    refL(i)=norm(nodes(i+1,:)-nodes(i,:)); % formulating refL: length of each section
end

% defining "Voronoi Length", used for bending and twisting
voronoiL=zeros(N,1);
for i=(1:N)
    if  i==1
        voronoiL(i)=1/2*refL(i); % first node half length back
    elseif i==N
        voronoiL(i)=1/2*refL(i-1); % last node half length front
    else
        voronoiL(i)=1/2*refL(i-1)+1/2*refL(i); % middle nodes whole length front and back
    end
end

%% Initializing reference frame
a1=zeros(ne,3); % initialize first reference director
a2=zeros(ne,3); % second reference director
tangent=getTangent(q0); % impliment function to get tangent of all direction for each node

% compute the reference director
t0=tangent(1,:); % tangent of the first edge
t1=[0;0;-1]; % creating an arbitrary vector
a1Tmp=cross(t0,t1); % finding a vector at right angle with t0
if abs(a1Tmp)<1e-6 % if somehow right angle vector is not found
    t1=[0,1,0]; % choose another arbitrary vector
    a1Tmp=cross(t0,t1);
end

% lebaling first section
a1(1,:)=a1Tmp/norm(a1Tmp); % getting unit vector thats perpendicular to tangent
a2(1,:)=cross(tangent(1,:),a1(1,:));

% lebal the rest
for i=2:ne
    t0=tangent(i-1,:); % last tangent
    t1=tangent(i,:); % current tangent
    a10=a1(i-1,:); % get privious a1
    a11=parallelTransport(a10,t0,t1); % move a10 from previous tangent to current tangent

    a1(i,:)=a11/norm(a11); % saving current section
    a2(i,:)=cross(t1,a1(i,:)); % getting current a2
end

%% building Material Frame
theta=q0(4:4:end); % initializing bending angle theta
[m1,m2]=computeMaterialDirectors(a1,a2,theta);

% initializing twist
refTwist=zeros(N,1);

% natural curve
kappaBar=getKappa(q0,m1,m2);

%% building simulation environment
Fg=zeros(ndof,1); % initializing gravitational force
for i=1:N % incerting g force to all node first three dof
    Fg([4*i-3,4*i-2,4*i-1])=massVector([4*i-3,4*i-2,4*i-1]).*g; 
end

% tolerence for solution finding
tol=E*I/l^2*(1e-6); % guess magnitude of EI/l

% DOF constraints
fixInd=1:7; % first two nodes and the section inbetween fixed
freeInd=8:ndof; % rest of the dof free

%% Running Simulation
currentTime=0;
endZ=zeros(Nsteps,1); % initializing space to record final length of each instant

for step=1:Nsteps
    fprintf('Current time = %f\n', currentTime);
    [q,u,a1,a2]=objfun(q0,u,a1,a2,freeInd,tol,refTwist);
    
    currentTime=currentTime+dt; % update time

    q0=q; % pass on new q as old q

    endZ(step)=q(end); % storing final length for plotting

    % plotting
    if mod(step,100)==0 % get the remainder = 0, plot 1 sec, 2 sec, 3 sec only to save time and resource
        theta=q(4:4:end); % get theta of all nodes
        [m1,m2]=computeMaterialDirectors(a1,a2,theta);  % get the director of m1 m2
        plotrod(q,a1,a2,m1,m2,currentTime); % visualize the rod of each second with info above
    end
end

% Visualization
figure(2);
timeArray = (1:1:Nsteps) * dt; 
plot( timeArray, endZ, 'ro-');
xlabel('Time, t [sec]');
ylabel('z-coordinate of last node, \delta_z [m]');



