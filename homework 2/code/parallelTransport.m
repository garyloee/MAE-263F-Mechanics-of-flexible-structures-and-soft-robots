function d=parallelTransport(u,t1,t2)
% transporting u from t1 to t2

% finding mutual normal vector
rightVector=cross(t1,t2);

if norm(rightVector)==0 % t1 and t2 are on the same vector direction
    d=u;
else
    % find unit right angle vector
    ru=rightVector/norm(rightVector);
    % cautious validation
    ru=ru-dot(ru,t1)*t1; % dot product of normal vectors should be 0
    ru=ru/norm(ru);
    ru=ru-dot(ru,t2)*t2; % dot product of normal vectors should be 0
    ru=ru/norm(ru);

    n1=cross(t1,ru);
    n2=cross(t2,ru);
    d=dot(u,t1)*t2+dot(u,n1)*n2+dot(u,ru)*ru;
end
end


    