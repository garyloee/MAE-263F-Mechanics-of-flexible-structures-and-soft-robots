function plotrod(q,a1,a2,m1,m2,ctime)
% Function to plot a rod structure in 3D space

% Calculate the number of vertices
nv=(length(q)+1)/4;

% Extract x, y, z coordinates from the generalized coordinates q
x1=q(1:4:end);
x2=q(2:4:end);
x3=q(3:4:end);

% Calculate the total length of the rod
L=sum(sqrt((x1(2:end)-x1(1:end-1)).^2+...
    (x2(2:end)-x2(1:end-1)).^2+...
    (x3(2:end)-x3(1:end-1)).^2)) / 3;

% Scale the parameters based on the length of the rod
a1=0.1*L*a1;
a2=0.1*L*a2;
m1=0.1*L*m1;
m2=0.1*L*m2;

% Create a 3D plot
h1=figure(1);
clf() % Clear the current figure

% Plot the rod vertices
plot3(x1,x2,x3,'ko-');
hold on

% Highlight the initial position with a red triangle
plot3(x1(1),x2(1),x3(1),'r^');

% Loop through each vertex and plot the associated vectors
for i=1:nv-1
    xa=q(4*i-3:4*i-1);
    xb=q(4*i+1:4*i+3);
    xp=(xa+xb)/2;
    
    % Plot vector a1 in blue
    p1=plot3([xp(1),xp(1)+a1(i,1)],[xp(2),xp(2)+a1(i,2)],...
        [xp(3),xp(3)+a1(i,3)],'b--','LineWidth',2);
    
    % Plot vector a2 in cyan
    p2=plot3([xp(1),xp(1)+a2(i,1)],[xp(2),xp(2)+a2(i,2)],...
        [xp(3),xp(3)+a2(i,3)],'c--','LineWidth',2);
    
    % Plot vector m1 in red
    p3=plot3([xp(1),xp(1)+m1(i,1)],[xp(2),xp(2)+m1(i,2)],...
        [xp(3),xp(3)+m1(i,3)],'r-');
    
    % Plot vector m2 in green
    p4=plot3([xp(1),xp(1)+m2(i,1)],[xp(2),xp(2)+m2(i,2)],...
        [xp(3),xp(3)+m2(i,3)],'g-');
end

% Release the hold on the plot
hold off

% Add legend for vectors
legend([p1,p2,p3,p4],'a_1','a_2','m_1','m_2');

% Add title with the current time
title(num2str(ctime,'t=%f'));

% Set axis properties
axis equal
view(3);
xlabel('x');
ylabel('y');
zlabel('z');

% Display the plot
drawnow
end
