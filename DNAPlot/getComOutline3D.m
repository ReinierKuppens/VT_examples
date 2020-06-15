


function [X,Y,Z,C]=getComOutline3D(P,r) 

% P = [1,1,1]
% r = 5; 

[X1,Y1,Z1] = sphere(40);

X = r.*X1 + P(1);
Y = r.*Y1 + P(2); 
Z = r.*Z1 + P(3);

C1 = [ones(10,10),zeros(10,10),ones(10,10),zeros(10,10)];
C2 = [zeros(10,10),ones(10,10),zeros(10,10),ones(10,10)];

C(:,:,1)  = [C1;C1;C2;C2];
C(:,:,2)  = [C1;C1;C2;C2];
C(:,:,3)  = [C1;C1;C2;C2];
% figure;surface(X,Y,Z,C);view(3);axis equal 
% colormap gray 
% xlabel('x')
% ylabel('y')
% zlabel('z')