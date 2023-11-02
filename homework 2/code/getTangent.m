function tangent=getTangent(q)
N=(length(q)+1)/4; % get N from length of q
tangent=zeros(N-1,3); % initialize tangent with size of(number of edges, 3 DOF)
for i=1:N-1
    gradient=q([4*i+1,4*i+2,4*i+3])-q([4*i-3,4*i-2,4*i-1]); % get the change in each dimension
    tangent(i,:)=gradient/norm(gradient); % finding unit tangent vector
end
end

