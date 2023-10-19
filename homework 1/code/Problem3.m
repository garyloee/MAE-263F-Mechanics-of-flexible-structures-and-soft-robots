clear variables
close all

% Declare the amount of Spheres (N)
N=50;
ndof=N*2; % number of degree of freedom

dt=0.01; % Time step (second)

rodLength=1; % Rod Length (meter)
dLength=rodLength/(N-1); % distance between spheres

% Radius of each sphere in meters
Ro=0.013; % outer radius
ri=0.011; % inner radius

% moment of Inertia
I=pi/4*(Ro^4-ri^4);


% insert force
P=2000; % Newton force
insertF=zeros(2*N,1); 
d=0.75; % meters from x=0
nodeP=round((N-1)*d)+1;
insertF(2*nodeP)=-P;

% Density of different material in kg/m^3
rho=2700;

g=9.8; 

visc=0; % viscosity µ

totalTime=1; %total simulation time in seconds

% Utility Parameters
YM=7e10; %GPa
ne=N-1; % number of edges, between all neighboring. spheres
EI=YM*I; % bending stiffness= Young's Modulus * Area moment of Inertia
EA=YM*pi*(Ro^2-ri^2); % Axial Rigidity, YM * area

% Initial Geometry Configuration, o-O-o
nodes=zeros(N,2); % N row 2 col zeros matrix
for n=1:N % entering location of each sphere
    nodes(n,1)=(n-1)*dLength; % x cordinate, (0,dL,2dL,...)
%    nodes(n,2)=0; % y cordinate (0,0,0,...)
end

% Mass matrix M, (x1,0;0,y1;...)
M=zeros(ndof,ndof);

for m =1:N
    M(2*m-1,2*m-1)=pi*(Ro^2-ri^2)*rodLength*rho/(N-1); % Mass of first sphere(volume * density)
    M(2*m,2*m)=pi*(Ro^2-ri^2)*rodLength*rho/(N-1); % same mass different axis
end



% Viscous damping matrix C
%{
C=zeros(ndof,ndof);
for m =1:N
    C(2*m-1,2*m-1)=6*pi*visc*R(m);
    C(2*m,2*m)=6*pi*visc*R(m);
end
%}

% no weight?
% Weight vector W (single column matrix pointing downward)

W=zeros(ndof,1);
for m =1:N
    W(2*m)=-pi*(Ro^2-ri^2)*rodLength*rho/(N-1)*g;
end




% Initial condition
q0=zeros(ndof,1); %location [x1;y1;x2;y2;...]
for n=1:N
    q0(2*n-1)=nodes(n,1); %x1,x2,x3
    q0(2*n)=nodes(n,2);%y1,y2,y3
end

q=q0; % start from last instant
u=(q-q0)/dt;

tol=EI/rodLength^2*1e-3; % small tolerance acounting for our magnitude, closing to zero

% Time stepping
steps=round(totalTime/dt)+1;

all_mid_y=zeros(steps,1);
all_mid_v=zeros(steps,1);

all_mid_y(1)=q(4);
all_mid_v(1)=u(4);

% create ymax matrix to store displacement over time
ymax=zeros(steps,1);

% setting up plot 1
plotPoints=[0, totalTime]; 
figure(1);
clf
fig=gcf;
fig.Position=[800 500 500 400];
grid on
ylim([-0.25 0.05])
xlabel('x [meter]');
ylabel('y [meter]');
title('Shape of Structure at Different Time')

for n=1:steps
    fprintf('Time = %f\n',(n-1)*dt);

    q=q0; %start point
    qfree=q(3:2*N-1);

    % finding solution to force equillibrium
    err=tol*10;
    while (err>tol && n>1) 
        f=M/dt*((q-q0)/dt-u)-insertF; % resulting force
        J=M/dt^2;
        
        % Elastic Forces
        for m=1:N-1 % linear or streching
            xk=q(2*m-1);
            yk=q(2*m);
            xkp1=q(2*m+1);
            ykp1=q(2*m+2);
            l_k=dLength;
            dF=gradEs(xk,yk,xkp1,ykp1,l_k,EA);
            dJ=hessEs(xk,yk,xkp1,ykp1,l_k,EA);
            f(2*m-1:2*m+2)=f(2*m-1:2*m+2)+dF;
            J(2*m-1:2*m+2,2*m-1:2*m+2)=J(2*m-1:2*m+2,2*m-1:2*m+2)+dJ;
        end
        
        for m=1:N-2 % Bending
            xkm1=q(2*m-1);
            ykm1=q(2*m);
            xk=q(2*m+1);
            yk=q(2*m+2);
            xkp1=q(2*m+3);
            ykp1=q(2*m+4);
            curvature0=0;
            l_k=dLength;
            dF=gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            dJ=hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);

            f(2*m-1:2*m+4)=f(2*m-1:2*m+4)+dF;
            J(2*m-1:2*m+4,2*m-1:2*m+4)=J(2*m-1:2*m+4,2*m-1:2*m+4)+dJ;

        end
        

        % Viscous Force
        %{
        f=f+C*(q-q0)/dt;
        J=J+C/dt;
        %}

        % Weight;
        f=f-W;
        
        ffree=f(3:2*N-1);
        Jfree=J(3:2*N-1,3:2*N-1);
        % seeking next q
        qfree=qfree-Jfree\ffree;
        q(3:2*N-1)=qfree;
        
        err=sum(abs(ffree)); % reach zero to exit
    end

    % passing current q u to next instant as 'old'
    u=(q-q0)/dt;
    q0=q;  
    % saving informaiton
    all_mid_y(n)=q(N);
    all_mid_v(n)=u(N);

    % finding and saving ymax
    ymax(n)=min(q);

    % Ploting
    
    
    if ismember((n-1)*dt,plotPoints)
        hold on
        leg=sprintf('Time = %2.2f',(n-1)*dt);
        plot(q(1:2:end),q(2:2:end),'.-','DisplayName',leg);
        
    end

    % saving angle
    %{
    v1=[q(3)-q(1),q(4)-q(2),0];
    v2=[q(5)-q(3),q(6)-q(4),0];
    
    angles(n)=atan2(norm(cross(v1,v2)),dot(v1,v2));
    %}


end

legend
hold off % finish plot 1

figure(2);
fig=gcf;
fig.Position=[800 100 1000 400];
subplot(1,2,1)
timeArray=(0:steps-1)*dt;
plot(timeArray,all_mid_y,'m-');
xlabel('Time, t [sec]');
ylabel('Position of mid-node, v [meter/sec]');
title('Mid-node y axis Position')

subplot(1,2,2)
plot(timeArray,all_mid_v,'-');
xlabel('Time, t [sec]');
ylabel('Velocity of mid-node, v [meter/sec]');
title('Mid-node y axis Velocity')

% diplaying terminal y velocity
m=(N+1)/2;
fprintf('Terminal y axis Velocity of Node %d : %f meter/sec\n',m, u(2*m))
    
figure(3)
fig=gcf;
fig.Position=[1300 500 500 400];
clf
plot(timeArray,ymax,'.-');
xlabel('Time, t [sec]');
ylabel('Maximum Vertical Displacement, y_m_a_x [meter]');
title('Maximum Vertical Displacement over Time')
hold on
plot(timeArray(n), ymax(n), 'r.', 'MarkerFaceColor','r')
hold off
text(timeArray(n), ymax(n), sprintf('(Steady State Displacement = %f m)',ymax(n)), 'Horiz','right', 'Vert','bottom')

% diplaying turning angle
%{
figure(3)
fig=gcf;
fig.Position=[1300 500 500 400];
plot(timeArray,angles,'r-');
xlabel('Time, t [sec]');
ylabel('Turning Angle, ° [degree]');
title('Mid-node Turning Angle')
%}





problem3Euler()








