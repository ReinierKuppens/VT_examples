







function [x1,y1,x2,y2,x3,y3]=getComOutline(p,r)
hold on

t1 = 0:0.01:pi/2;
t2 = pi:0.01:3/2*pi;
t3 = 0:0.01:2*pi;

x1 = r*cos(t1)+p(1);
y1 = r*sin(t1)+p(2);
x2 = r*cos(t2)+p(1);
y2 = r*sin(t2)+p(2);
x3 = r*cos(t3)+p(1);
y3 = r*sin(t3)+p(2);

x1 = [p(1) x1];
y1 = [p(2) y1];

x2 = [p(1) x2];
y2 = [p(2) y2];


% fill(x3,y3,'w')
% plot(x3,y3,'k')
% fill(x1,y1,'k')
% fill(x2,y2,'k')

end%Function

