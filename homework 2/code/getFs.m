function [Fs,Js]=getFs(q)
global E A refL
N=(length(q)+1)/4;
ne=N-1;
Fs=zeros(size(q));
Js=zeros(length(q),length(q));
for i=1:ne % run through all sections

    node1=transpose(q(4*i-3:4*i-1)); % current node
    node2=transpose(q(4*i+1:4*i+3)); % next node
 
    % get newton method index
    [dF,dJ]=gradEs_hessEs(node1,node2,refL(i),E*A);
 
    ind=[4*i-3,4*i-2,4*i-1,4*i+1,4*i+2,4*i+3]; % two nodes current and next
    Fs(ind)=Fs(ind)-dF;
    Js(ind,ind)=Js(ind,ind)-dJ;
 
end
end