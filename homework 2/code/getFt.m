function [Ft,Jt]=getFt(q, refTwist)
global G J voronoiL

N=(length(q)+1)/4;
Ft=zeros(size(q));
Jt=zeros(length(q),length(q));

for i=2:N-1 % Compute bending force at each internel node
 
    node0=transpose(q(4*i-7:4*i-5));
    node1=transpose(q(4*i-3:4*i-1));
    node2=transpose(q(4*i+1:4*i+3));
    thetaE=q(4*i-4);
    thetaF=q(4*i);
 
    [dF,dJ]=gradEt_hessEt(node0,node1,node2,...
    thetaE,thetaF,refTwist(i),voronoiL(i),G*J);
 
    ind=4*i-7:4*i+3; % previous current next node
    Ft(ind)=Ft(ind)-dF;
    Jt(ind,ind)=Jt(ind,ind)-dJ;
 
end
end