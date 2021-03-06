function Es1 = MGspringenergy(in1,in2,L0,K,in5,in6)
%MGSPRINGENERGY
%    ES1 = MGSPRINGENERGY(IN1,IN2,L0,K,IN5,IN6)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    25-Dec-2017 16:52:50

com11 = in5(1,:);
com12 = in5(2,:);
ps11 = in1(1,:);
ps12 = in1(2,:);
ps21 = in2(1,:);
ps22 = in2(2,:);
theta1 = in6(3,:);
x1 = in6(1,:);
y1 = in6(2,:);
t3 = cos(theta1);
t4 = com12-ps22;
t5 = sin(theta1);
t6 = com11-ps21;
t2 = ps11-x1+t3.*t6-t4.*t5;
t7 = ps12-y1+t3.*t4+t5.*t6;
t8 = L0-sqrt(t2.^2+t7.^2);
Es1 = K.*t8.^2.*(1.0./2.0);
