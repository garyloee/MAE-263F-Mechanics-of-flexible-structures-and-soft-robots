function [q,u,a1,a2]=objfun(q0,u,a1,a2,freeInd,tol,refTwist)

% get same global variables just to make it easier
global Fg M dt

q=q0; % initial approximation
iter=1; % initialize iteration

% setup simulation environment
err=10*tol;
while err>tol % below tolerence to locate solution
    % compute reference frame
    [a1iter,a2iter]=computeTimeParallel(a1,q0,q);

    % compute reference Twist
    tangent=getTangent(q);
    refTwistIter=computeRefTwist(a1iter,tangent,refTwist);

    % compute material frame
    theta=q(4:4:end); % getting bending angle of all sections
    [m1,m2]=computeMaterialDirectors(a1iter,a2iter,theta);

    % compute f and J
    [Fb,Jb]=getFb(q,m1,m2); % bending force
    [Ft,Jt]=getFt(q,refTwistIter);% get twisting force
    [Fs,Js]=getFs(q); % stretching force

    % Equation of motion
    f=M/dt*((q-q0)/dt-u)-Fb-Ft-Fs-Fg;
    J=M/dt^2-Jb-Jt-Js;

    ffree=f(freeInd); % only doing Newton meathod on free nodes
    Jfree=J(freeInd,freeInd); % essentially fixing the first two nodes

    % updating q
    qfree=Jfree \ ffree;
    q(freeInd)=q(freeInd)-qfree; % similar to HW1

    % converging
    err=sum(abs(ffree));

    fprintf('iter = %d, error = %f\n', iter, err);
    iter=iter+1;
end

% getting new u velocity
u=(q-q0)/dt;
a1=a1iter; % passing on informaiton 
a2=a2iter;

end
