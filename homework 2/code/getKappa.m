function kappa=getKappa(q,m1,m2)
% calculating kappa 1 and kappa 2
kappa=zeros(length(m1)+1, 2); % size of (N, 2 kappa)

for i=2:length(m1) % avoiding first and last
    nim1=q([4*i-7,4*i-6,4*i-5]); % finding nodes before and after current node
    ni=q([4*i-3,4*i-2,4*i-1]);
    nip1=q([4*i+1,4*i+2,4*i+3]);

    m10=m1(i-1,:); % getting m1 m2 vector of previous edge
    m11=m1(i,:);

    m20=m2(i-1,:);
    m21=m2(i,:);

    currentKappa=computekappa(nim1,ni,nip1,m10,m20,m11,m21); % return [current kappa1, current kappa2]

    kappa(i,1)=currentKappa(1); % saving data
    kappa(i,2)=currentKappa(2);
end
end