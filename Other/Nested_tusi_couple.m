

clear all
close all
clc



alpha = 0:0.5:2*pi;
beta  = 0:0.5:2*pi;
gamma = 0:0.5:2*pi;
zeta  = gamma*0.1;

r1    = 3;
r2    = 2;
r3    = 1;

T01 = @(a,b,r) [cos(b) -sin(b) r*cos(a);
    sin(b)  cos(b) r*sin(a);
    0          0         1];

O3  = [ r2;0];
O4  = [-r2;0];

p1 = [-r3; 0];
p2 = [ r3; 0];
p3 = [-r3; 0];
p4 = [ r3; 0];

count = 1;

for k = 1:numel(alpha)
    for ii = 1:numel(beta)
        for jj = 1:numel(gamma)
            
            O2_prime(:,count) = T01(alpha(k),0,r1)*[0;0;1];
            O3_prime(:,count) = T01(alpha(k),beta(ii),r1)*[O3;1];
            O4_prime(:,count) = T01(alpha(k),beta(ii),r1)*[O4;1];
            
            p1_prime(:,count) = T01(0,gamma(jj),0)*[p1;1] + T01(alpha(k),beta(ii),r1)*[O3;1];
            p2_prime(:,count) = T01(0,gamma(jj),0)*[p2;1] + T01(alpha(k),beta(ii),r1)*[O3;1];
            
            p3_prime(:,count) = T01(0,zeta(jj),0)*[p3;1] + T01(alpha(k),beta(ii),r1)*[O4;1];
            p4_prime(:,count) = T01(0,zeta(jj),0)*[p4;1] + T01(alpha(k),beta(ii),r1)*[O4;1];
            
                
            a(count,1) = alpha(k);
            b(count,1) = beta(ii);
            c(count,1) = gamma(jj);
            d(count,1) = zeta(jj);
                        count = count + 1;

        end
    end
end

% keyboard 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                     Make DNA                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DNA           = initializeDNA;

DNA.incstr(1) = inc2str([1;1])
DNA.incstr(2) = inc2str([0;1;1])
DNA.incstr(3) = inc2str([0;0;1;1])
DNA.incstr(4) = inc2str([0;0;1;0;1]);

DNA.incstr(5) = inc2str([1;0;0;1])
DNA.incstr(6) = inc2str([1;0;0;1])
DNA.incstr(7) = inc2str([1;0;0;0;1])
DNA.incstr(8) = inc2str([1;0;0;0;1])


DNA.edgelabel = [1 1 1 1 2 2 2 2]
% DNA.edgelabel = [1 1 1 1 ];

DNA.Mpar(1,:)   = NaN;
DNA.Mpar(2:5,:) = ones(4,3);

DNA.Hpar  = [0  O2_prime(1,1) O3_prime(1,1) O4_prime(1,1)
    0  O2_prime(2,1) O3_prime(2,1) O4_prime(2,1)];

DNA.Spar = [zeros(2,4)
    p1_prime(1:2,1) p2_prime(1:2,1) p3_prime(1:2,1) p4_prime(1:2,1)
    zeros(1,4)
    ones(1,4)];

t = 1:numel(O2_prime(1,:));

[~,comPoint] = getMassLocations(DNA);

m1 = comPoint{1}.';
m2 = comPoint{2}.';
m3 = comPoint{3}.';
m4 = comPoint{4}.';

count = 1; 
for k = 1:numel(alpha)
    for ii = 1:numel(beta)
        for jj = 1:numel(gamma)
            
            M1_prime(:,count) = T01(alpha(k),0,r1)*[m1;1];
            
            O2_prime(:,count) = T01(alpha(k),0,r1)*[m1;1];
            O3_prime(:,count) = T01(alpha(k),beta(ii),r1)*[O3;1];
            O4_prime(:,count) = T01(alpha(k),beta(ii),r1)*[O4;1];
            
            p1_prime(:,count) = T01(0,gamma(jj),0)*[p1;1] + T01(alpha(k),beta(ii),r1)*[O3;1];
            p2_prime(:,count) = T01(0,gamma(jj),0)*[p2;1] + T01(alpha(k),beta(ii),r1)*[O3;1];
            
            p3_prime(:,count) = T01(0,zeta(jj),0)*[p3;1] + T01(alpha(k),beta(ii),r1)*[O4;1];
            p4_prime(:,count) = T01(0,zeta(jj),0)*[p4;1] + T01(alpha(k),beta(ii),r1)*[O4;1];
                                    count = count + 1;

            
        end
    end
end

% keyboard 

state = [M1_prime(1:2,:).', a,p2_prime(1:2,:).', b,p3_prime(1:2,:).', c,p4_prime(1:2,:).', d];

animateDNA(DNA,t,state) 










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p1Squared = p1_prime(1,:).^2 + p1_prime(2,:).^2;
p2Squared = p2_prime(1,:).^2 + p2_prime(2,:).^2;
p3Squared = p3_prime(1,:).^2 + p3_prime(2,:).^2;
p4Squared = p4_prime(1,:).^2 + p4_prime(2,:).^2;


figure;hold on
plot(p1Squared);
plot(p2Squared);
plot(p3Squared);
plot(p4Squared);

plot(p1Squared+p2Squared+p3Squared+p4Squared)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;axis([-6,6,-6,6]);axis square
hold on;
h1 = plot([O2_prime(1,1),O3_prime(1,1)],[O2_prime(2,1),O3_prime(2,1)]);
h2 = plot([O2_prime(1,1),O4_prime(1,1)],[O2_prime(2,1),O4_prime(2,1)]);

h3 = plot([O3_prime(1,1),p1_prime(1,1)],[O3_prime(2,1),p1_prime(2,1)]);
h4 = plot([O3_prime(1,1),p2_prime(1,1)],[O3_prime(2,1),p2_prime(2,1)]);
h5 = plot([O4_prime(1,1),p3_prime(1,1)],[O4_prime(2,1),p3_prime(2,1)]);
h6 = plot([O4_prime(1,1),p4_prime(1,1)],[O4_prime(2,1),p4_prime(2,1)]);

plot(r1.*cos(alpha),r1.*sin(alpha))


for k = 1:10:numel(O3_prime(1,:))
    set(h1,'Xdata',[O2_prime(1,k),O3_prime(1,k)],'Ydata',[O2_prime(2,k),O3_prime(2,k)])
    set(h2,'Xdata',[O2_prime(1,k),O4_prime(1,k)],'Ydata',[O2_prime(2,k),O4_prime(2,k)])
    
    set(h3,'Xdata',[O3_prime(1,k),p1_prime(1,k)],'Ydata',[O3_prime(2,k),p1_prime(2,k)])
    set(h4,'Xdata',[O3_prime(1,k),p2_prime(1,k)],'Ydata',[O3_prime(2,k),p2_prime(2,k)])
    set(h5,'Xdata',[O4_prime(1,k),p3_prime(1,k)],'Ydata',[O4_prime(2,k),p3_prime(2,k)])
    set(h6,'Xdata',[O4_prime(1,k),p4_prime(1,k)],'Ydata',[O4_prime(2,k),p4_prime(2,k)])
    
    
    pause(0.1)
end











