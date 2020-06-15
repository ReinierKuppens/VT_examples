

clear all 
close all 
clc 

alpha = linspace(0,2*pi,400)


t1 = 0; 
t2 = 1; 
t3 = 2;

rho_1 = 1; 
rho_2 = 1; 
rho_3 = 1; 

T = @(a,tx)[cos(a),-sin(a),tx;
            sin(a), cos(a),0 ;
            0     , 0     ,1 ];
        
       
DNA           = initializeDNA;
DNA.incstr    = [1 3 6 ];
DNA.edgelabel = [1 1 1 ];

m1 = rand; 
m2 = rand;
m3 = rand;


com1 = rand(2,1);
com2 = rand(2,1);
com3 = rand(2,1);


DNA.Hpar = [[t1;0],[t2+t1;0],[t3+t2+t1;0]];



for k = 1:numel(alpha) 

    C1_prime{k} = T(rho_1*alpha(k),t1)*[com1;1];
    C2_prime{k} = T(rho_1*alpha(k),t1)*T(rho_2*alpha(k),t2)*[com2;1];
    C3_prime{k} = T(rho_1*alpha(k),t1)*T(rho_2*alpha(k),t2)*T(rho_3*alpha(k),t3)*[com3;1];
   
end

DNA.Mpar(2,:) = [C1_prime{1}(1:2,:).' m1]
DNA.Mpar(3,:) = [C2_prime{1}(1:2,:).' m2]
DNA.Mpar(4,:) = [C3_prime{1}(1:2,:).' m3] 


for k = 1:numel(alpha) 
    
state(k,:) = [  C1_prime{k}(1:2,:).', (rho_1)*alpha(k),...
                C2_prime{k}(1:2,:).', (rho_2+rho_1)*alpha(k),...
                C3_prime{k}(1:2,:).', (rho_3+rho_2+rho_1)*alpha(k)];

end


t = alpha;



animateDNA(DNA,t,state) 

















