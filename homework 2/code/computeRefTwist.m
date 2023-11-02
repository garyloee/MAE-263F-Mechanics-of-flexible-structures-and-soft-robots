function refTwist=computeRefTwist(a1,tangent,refTwist)

% setup
% get ne
[ne, ~] = size(a1);
%ne=length(a1);
for i=2:ne % run through all middle nodes
    
    u0=a1(i-1,:); % get previous a1
    u1=a1(i, :); % get current a1
    t0=tangent(i-1,:); % get previous tangent
    t1=tangent(i,:); % get current tangent
    ut=parallelTransport(u0, t0, t1); % purform parallel transport
    % refTwist(i) = signedAngle(ut, u1, t1); 
    % accounting for u1 = ut
    % error when rotate >360 degree
    ut=rotateAxisAngle(ut,t1,refTwist(i));
    refTwist(i)=refTwist(i)+signedAngle(ut,u1,t1); 
    
end
end
