clear variables
close all
global dt q N q0 tol M u dLength EA EI Ro HRo visc mpNodes ndof C

%% Initialize environment and parameters
visc = 0.1; % Viscosity [Pa-s]

% Defining object parameters
N = 50; % Number of nodes
ndof = N * 2; % Number of degree of freedom
mpNodes=15; % Number of midpiece nodes
rodLength = 50e-2; % Rod length [m]
dLength = rodLength / (N - 1); % Nodal distance between spheres
Ro = 0.5e-2; % outer radius [m]
HRo=Ro*3; % head radius [m]
flgRo=Ro/2; % flagellum radius [m]
I = pi / 2 * (Ro^4); % Moment of inertia [kg/m^2]
rho = 1000; % Density [kg/m^3]

% Utility Parameters
YM = 5e6; % Young's Modulus [Pa]
ne = N - 1; % number of edges, between all neighboring. spheres
EI = YM * I; % bending stiffness = Young's Modulus * Area moment of Inertia
EA = YM * pi * Ro^2; % Axial Rigidity, YM * area

% Initial Geometry Configuration
nodes = zeros(N, 2); % N row 2 col zeros matrix
for n = 1:N % entering location of each sphere
    nodes(n, 1) = (n - 1) * dLength; % x cordinate, (0,dL,2dL,...)
end

% Mass matrix M, (x1,0;0,y1;...)
M = zeros(ndof, ndof);
for m = 1:mpNodes
    M(2 * m - 1, 2 * m - 1) = pi * Ro^2 * dLength * rho; % Mass of first sphere(volume * density)
    M(2 * m, 2 * m) = pi * Ro^2 * dLength * rho; % same mass different axis
end
for m = mpNodes+1:N
    M(2 * m - 1, 2 * m - 1) = pi * flgRo^2 * dLength * rho; % Mass of first sphere(volume * density)
    M(2 * m, 2 * m) = pi * flgRo^2 * dLength * rho; % same mass different axis
end

M(1,1) = M(1,1)+3/4*pi * HRo^3 * rho / (N - 1); 
M(2,2) = M(2,2)+3/4*pi * HRo^3 * rho / (N - 1); 

% Viscous damping matrix C
C = zeros(ndof, ndof);
C(1,1) = 6 * pi * visc * HRo^3;
C(2,2) = 6 * pi * visc * HRo^3;

% Initial conditions
q0 = zeros(ndof, 1); %location [x1;y1;x2;y2;...]
for n = 1:N
    q0(2 * n - 1) = nodes(n, 1); %x1,x2,x3
    q0(2 * n) = nodes(n, 2); %y1,y2,y3
end

q = q0; % start from last instant
u = zeros(size(q));

% simulation variable
tol = EI / rodLength^2 * 1e-3; % small tolerance accounting for our magnitude, closing to zero

% Time stepping paramters
totalTime = 15; % Total simulation time [s]
dt = 0.005; % Time step [s]
steps = round(totalTime / dt) + 1; % total time*100

%% Running and designing simulation
% simulating path control uncomment line 74-82 and comment 84
% (the velocity and displacement graph would not funtion correctly)

%{
[qdata,vdata]=swim(1, 500, 0.2, 0.1, 0);
[qdata,vdata]=swim(501, 590, 0.2, 0.1, 0.5);
[qdata,vdata]=swim(591, 1090, 0.2, 0.1,0);
[qdata,vdata]=swim(1091, 1180, 0.2, 0.1,0.5);
[qdata,vdata]=swim(1181, 1680, 0.2, 0.1,0);
[qdata,vdata]=swim(1681, 1770, 0.2, 0.1,0.5);
[qdata,vdata]=swim(1771, 2270, 0.2, 0.1,0);
[qdata,vdata]=swim(2271, 2360, 0.2, 0.1,0.5);
[qdata,vdata]=swim(2361, 2861, 0.2, 0.1,0);
%}
[qdata,vdata]=swim(1, steps, 0.2, 0.1,0);
Qdata=zeros(1,steps);
Vdata=zeros(1,steps);
timestep=linspace(0,totalTime,steps);
for i=1:steps
    Qdata(i)=qdata(2*i-1);
    Vdata(i)=vdata(2*i-1);
end
figure(2)

clf
fig = gcf;
fig.Position = [400 200 600 600];
plot(timestep,Qdata, 'r.-');
xlabel('time [s]');
ylabel('displacement [m]')
title('Displacement vs Time')
figure(3)

clf
fig = gcf;
fig.Position = [800 200 600 600];
plot(timestep,Vdata, 'r.-');
xlabel('time [s]');
ylabel('velocity [m/s]')
title('Velocity vs Time')
localMax=max(Vdata(end-50:end));
localmin=min(Vdata(end-50:end));

fprintf('terminal velocity = %d', (localMax+localmin)/2);
%% defining swimming function to make coding easier
function [qdata,vdata]=swim(nstart, nend, F, Period, turn)
    global dt q N q0 tol M u dLength EA EI visc mpNodes C
    qdata=zeros(2*(nend-nstart),1);

    for n = nstart:nend
        if rem(n-1,100)==0
            fprintf('Time = %f\n', (n - 1) * dt);
        end

        q = q0; % Start point
        actuate = F * sin(2 * pi * n * dt / Period)+turn*F; % Change from time to frequency
        qfree=q(5:2*N);

        % Finding solution to force equilibrium
        err = tol * 10;
        while (err > tol && n > 1)
            f = M / dt * ((q - q0) / dt - u);
            J = M / dt^2;

            % Elastic Forces
            for m = 1:mpNodes % linear or stretching
                xk = q(2 * m - 1);
                yk = q(2 * m);
                xkp1 = q(2 * m + 1);
                ykp1 = q(2 * m + 2);
                l_k = dLength;
                dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
                dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
                f(2 * m - 1:2 * m + 2) = f(2 * m - 1:2 * m + 2) + dF;
                J(2 * m - 1:2 * m + 2, 2 * m - 1:2 * m + 2) = J(2 * m - 1:2 * m + 2, 2 * m - 1:2 * m + 2) + dJ;
            end

            for m = mpNodes+1:N - 1 % linear or stretching
                xk = q(2 * m - 1);
                yk = q(2 * m);
                xkp1 = q(2 * m + 1);
                ykp1 = q(2 * m + 2);
                l_k = dLength;
                dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA/4);
                dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA/4);
                f(2 * m - 1:2 * m + 2) = f(2 * m - 1:2 * m + 2) + dF;
                J(2 * m - 1:2 * m + 2, 2 * m - 1:2 * m + 2) = J(2 * m - 1:2 * m + 2, 2 * m - 1:2 * m + 2) + dJ;
            end

            for m = 1:mpNodes % Bending
                xkm1 = q(2 * m - 1);
                ykm1 = q(2 * m);
                xk = q(2 * m + 1);
                yk = q(2 * m + 2);
                xkp1 = q(2 * m + 3);
                ykp1 = q(2 * m + 4);
                curvature0 = actuate*(m/mpNodes); % actuating by changing natural curvature
                l_k = dLength;
                dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
                dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);

                f(2 * m - 1:2 * m + 4) = f(2 * m - 1:2 * m + 4) + dF;
                J(2 * m - 1:2 * m + 4, 2 * m - 1:2 * m + 4) = J(2 * m - 1:2 * m + 4, 2 * m - 1:2 * m + 4) + dJ;

            end

            for m = mpNodes+1:N - 2 % Bending
                xkm1 = q(2 * m - 1);
                ykm1 = q(2 * m);
                xk = q(2 * m + 1);
                yk = q(2 * m + 2);
                xkp1 = q(2 * m + 3);
                ykp1 = q(2 * m + 4);
                curvature0 = 0; % actuating by changing natural curvature
                l_k = dLength;
                dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI/16);
                dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI/16);

                f(2 * m - 1:2 * m + 4) = f(2 * m - 1:2 * m + 4) + dF;
                J(2 * m - 1:2 * m + 4, 2 * m - 1:2 * m + 4) = J(2 * m - 1:2 * m + 4, 2 * m - 1:2 * m + 4) + dJ;


            end

            % Parameters for Gray and Hancock model
            a = 0.25e-6; % Characteristic radius of the sperm tail [m]
            lamda = 10e-6; % Wavelength of the flagellar beat [m]

            % Calculate drag forces
            for m = 1:N

                % Calculate drag coefficients from Gray and Hancock model
                Ct = (2 * pi * visc) / (log(2 * lamda / a) - 0.5);
                Cn = (4 * pi * visc) / (log(2 * lamda / a) + 0.5);

                % calculate node tangent
                if m == 1
                    tangent=[q(2 * m + 1) - q(2 * m-1),q(2*m+2)-q(2*m)];
                    tang=(tangent./norm(tangent))';

                elseif m == N
                    tangent=[q(2 * m-1) - q(2 * m - 3),q(2 * m )-q(2 * m - 2)];
                    tang=(tangent./norm(tangent))';
                    
                else
                    tangent=[q(2 * m + 1) - q(2 * m - 3),q(2 * m + 2)-q(2 * m - 2)];
                    tang=(tangent./norm(tangent))';
                end

                if m==mpNodes
                    a = 0.125e-6; % Characteristic radius of the sperm tail [m]
                    lamda = 10e-6; % Wavelength of the flagellar beat [m]

                end

                % Tangential drag force
                f(2 * m - 1:2*m) = f(2 * m - 1:2*m) + dLength * Ct * tang * (tang') * (q(2*m-1:2*m) - q0(2*m-1:2*m)) / dt;
                
                J(2 * m - 1:2*m,2 * m - 1:2*m) = J(2 * m - 1:2*m,2 * m - 1:2*m) + dLength * Ct * tang * (tang') / dt;
                
                % Normal drag force
                f(2 * m - 1:2*m) = f(2 * m - 1:2*m) + dLength * Cn * ( eye(2) - tang * (tang')) * (q(2*m-1:2*m) - q0(2*m-1:2*m)) / dt;
                
                J(2 * m - 1:2*m,2 * m - 1:2*m) = J(2 * m - 1:2*m,2 * m - 1:2*m) + dLength * Cn * ( eye(2) - tang * (tang')) / dt;
                
            end

            % Viscous Force
            f = f + C * (q - q0) / dt;
            
            J = J + C / dt;

            % Seeking next q
            q = q - J \ f;

            err = sum(abs(f)); % reach zero to exit
        end
        
        
        % passing current q u to next instant as 'old'
        u = (q - q0) / dt;
        q0 = q;
        
        qdata(2*n-1:2*n)=q(1:2);
        vdata(2*n-1:2*n)=u(1:2);
        
        % Plotting
        figure(1);
        clf
        fig = gcf;
        fig.Position = [200 200 600 600];
        plot(q(1:2:end), q(2:2:end), 'r.-');
        hold on
        plot(q(1),q(2),'o', 'MarkerFaceColor','r');
        xlim([-0.6 0.5])
        ylim([-0.2 0.9])
        grid on
        drawnow
    end
end
