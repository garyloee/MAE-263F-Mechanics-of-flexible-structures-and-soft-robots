close all

clear variables
% Declare the amount of Spheres (N)
N=50;
ndof=N*2; % number of degree of freedom

dt=0.01; % Time step (second)

rodLength=0.1; % Rod Length (meter)
dLength=rodLength/(N-1); % distance between spheres

% Profile of the beam
b=0.01;
h=0.002;

% moment of Inertia
I=b*h^3/12;

% Density of different material in kg/m^3
rho=2700;

g=9.8; 

% Viscous damping matrix C
visc=0; % viscosity µ
% C=zeros(ndof,ndof);
% for m =1:N
%     C(2*m-1,2*m-1)=6*pi*visc*R(m);
%     C(2*m,2*m)=6*pi*visc*R(m);
% end



% Weight vector W (single column matrix pointing downward)
W=zeros(ndof,1);
for m =1:N
    W(2*m)=-b*h*rodLength*rho/(N-1)*g;
end

% Utility Parameters
YM=200e9; %GPa
ne=N-1; % number of edges, between all neighboring. spheres
EI=YM*I; % bending stiffness= Young's Modulus * Area moment of Inertia
EA=YM*b*h; % Axial Rigidity, YM * area

% Mass matrix M, (x1,0;0,y1;...)
M=zeros(ndof,ndof);

for m =1:N
    M(2*m-1,2*m-1)=b*h*rodLength*rho/(N-1); % Mass of first sphere(volume * density)
    M(2*m,2*m)=b*h*rodLength*rho/(N-1); % same mass different axis
end

% initialize force insertion matrix
insertF=zeros(2*N,1); 

% initialize inserted P magnitude
Psteps=[1,5:5:200];

% create ymax matrix to store displacement over varying P
ymax=zeros(length(Psteps),1);
for p=1:length(Psteps)
  
    % insert force
    P=Psteps(p); % Newton force
    insertF(2*N)=-P;
        
    totalTime=1; %total simulation time in seconds
        
    % Initial Geometry Configuration, o-O-o
    nodes=zeros(N,2); % N row 2 col zeros matrix
    for n=1:N % entering location of each sphere
        nodes(n,1)=(n-1)*dLength; % x cordinate, (0,dL,2dL,...)
    
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
    
    qmax=zeros(steps,1);
    
    for n=1:steps
        %fprintf('Time = %f\n',(n-1)*dt);
    
        q=q0; %start point
        qfree=q(5:2*N);
    
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
            
            ffree=f(5:2*N);
            Jfree=J(5:2*N,5:2*N);
            % seeking next q
            qfree=qfree-Jfree\ffree;
            q(5:2*N)=qfree;
            
            err=sum(abs(ffree)); % reach zero to exit
        end
    
        % passing current q u to next instant as 'old'
        u=(q-q0)/dt;
        q0=q;  
        % saving informaiton
        qmax(n)=q(2*N);
        
    end
    if Psteps(p)==10
        figure(3);
        fig=gcf;
        fig.Position=[800 100 500 400];
        
        timeArray=(0:steps-1)*dt;
        plot(timeArray,qmax,'m-');
        xlabel('Time, t [sec]');
        ylabel('Displacement, x [meter]');
        title('Displacement over time at P=10')
    end

    % finding and saving ymax
    ymax(p)=min(q);
    
    fprintf('Max Displacement when P = %d : %f meters\n',Psteps(p), ymax(p))
end

% Calculating euler
euymax=zeros(length(Psteps),1);
bk=0;
cc=rodLength;
for z=1:length(Psteps)
    pp=Psteps(z);
    euymax(z)=-pp*rodLength^3/(3*EI);% minus sign to acount for axis
    dif=(euymax(z)-ymax(z))/ymax(z);
    
    if dif>0.1 && bk==0
        diverg=[pp,dif];
        bk=1;
    end
end
if bk==1
    fprintf('DER and Euler deviates %.1f %% at P = %d\n',diverg(2)*100,diverg(1));
else
    fprintf('Did not deviate more than 10 %%\n');
end
% plotting ymax vs P 
figure(2)
clf
kk1=plot(Psteps,ymax,'.-');
hold on
kk2=plot(Psteps,euymax,'r.-');
xlabel('External Force Magnitude, P [N]');
ylabel('Maximum Vertical Displacement, y_m_a_x [meter]');
ylim([min(euymax) 0])
title('Maximum Vertical Displacement under Different Load')
plot(Psteps(p), ymax(p), 'b.', 'MarkerFaceColor','r')
plot(Psteps(p), euymax(p), 'm.', 'MarkerFaceColor','r')
legend([kk1,kk2],'Discrete simulation','Euler Beam Theory')
hold off
text(Psteps(p), ymax(p), sprintf('(Simulation Displacement = %f m)',ymax(p)), 'Horiz','right', 'Vert','bottom')
text(Psteps(p), euymax(p), sprintf('(Euler Theory Displacement = %f m)',euymax(p)), 'Horiz','right', 'Vert','bottom')
















