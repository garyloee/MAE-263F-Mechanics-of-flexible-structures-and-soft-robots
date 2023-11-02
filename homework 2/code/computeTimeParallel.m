function [a1,a2]=computeTimeParallel(a10,q0,q)


ne=(length(q)+1)/4-1;
tangent0=getTangent(q0); % previous tangent
tangent=getTangent(q); % current tangent

a1=zeros(ne,3); % initializing space for new a1 a2
a2=zeros(ne,3);

for i=1:ne % go though all sections

    t0=tangent0(i,:); % get previous tangent
    t=tangent(i,:); % get current tangent
 
    a1Previous=a10(i,:); % get previous a1 at this node
    a1Current=parallelTransport(a1Previous,t0,t); % move old to new

    % make sure dot product = 0
    a1Current=a1Current-dot(a1Current,t)*t;
    a1Current=a1Current/norm(a1Current);
 
    a1(i,:)=a1Current; % save current data
    a2(i,:)=cross(t,a1Current);
end
end