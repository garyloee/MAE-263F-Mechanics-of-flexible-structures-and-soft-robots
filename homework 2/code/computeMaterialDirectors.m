function [m1,m2]=computeMaterialDirectors(a1,a2,theta)
% initialize m1 m2
m1=zeros(size(a1));
m2=zeros(size(a1));

for i=1:length(theta)
    costi=cos(theta(i)); % getting cosine and sine of gending angle
    sinti=sin(theta(i));
    m1(i,:)=costi*a1(i,:)+sinti*a2(i,:); % getting the material directors of two axis
    m2(i,:)=costi*a2(i,:)-sinti*a1(i,:);
end
end