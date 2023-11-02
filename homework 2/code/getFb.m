function [Fb,Jb]=getFb(q,m1,m2)
global kappaBar E I voronoiL
N=(length(q)+1)/4;
Fb=zeros(size(q));
Jb=zeros(length(q),length(q));
for i=2:N-1 % Compute bending force at each internel node
 
    node0=transpose(q(4*i-7:4*i-5)); % get previous node
    node1=transpose(q(4*i-3:4*i-1)); % get current node
    node2=transpose(q(4*i+1:4*i+3)); % get next node
 
    m1p=m1(i-1,:); % m1 of previous edge 
    m1c=m1(i,:); % m1 of current edge
    m2p=m2(i-1,:); % m2 of previous edge
    m2c=m2(i,:); % m2 of current edge
 
    [dF,dJ]=gradEb_hessEb(node0,node1,node2,m1p,m2p,m1c,m2c,...
        kappaBar(i,:),voronoiL(i),E*I);
 
    ind=4*i-7:4*i+3; % the three nodes
    Fb(ind)=Fb(ind)-dF;
    Jb(ind,ind)=Jb(ind,ind)-dJ;
 
end
end