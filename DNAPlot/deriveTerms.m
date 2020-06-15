%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 29-06-2015                                             %
%   Last modified: 18-10-2016                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Derive the constraint equations
%

clear all; close all; clc; 

% % cf = 'D:\pkuppens\Box Sync\Reinier PhD\3_Matlab\Mechanical_EA_Toolbox_PRKuppens_19-12-2017\Essentials\simulation\deriveTerms';
% % 
% % cd(cf)

% if ~strcmp(cf,pwd)
%     cd('../Mechanical_engine/deriveTerms')
% end

%Declare all symbolic variables

x  = sym('x%d',[2,1]);   y = sym('y%d',[2,1]);   theta = sym('theta%d',[2,1]);   
xd = sym('xd%d',[2,1]); yd = sym('yd%d',[2,1]); thetad = sym('thetad%d',[2,1]); %thetad1,...,thetadn

com1    = sym('com1%d',[2,1]); com2     = sym('com2%d',[2,1]);
ps1     = sym('ps1%d',[2,1]);  ps2      = sym('ps2%d',[2,1]);
H       = sym('H1%d',[2,1]);   

syms K L0 theta0 t alpha a w

%Rotation matrix
R = @(theta) [cos(theta),-sin(theta);sin(theta),cos(theta)];

%State vectors
q1  = [x(1); y(1); theta(1)];
qd1 = [xd(1); yd(1); thetad(1)];

q2 = [x(1); y(1); theta(1); x(2); y(2); theta(2)];
qd2 = [xd(1); yd(1); thetad(1); xd(2); yd(2); thetad(2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       
%                       Scleronomic constraint equations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Making symbolic equations')

%Scleronomic rotational hinge constraint Ground-Mass
C1      = [x(1);y(1)]-(R(theta(1))*(com1-H)+H);       % zeroth derivative
C1_q    = jacobian(C1,q1);                            % First derviatives (d)
C1_qq   = jacobian(C1_q*qd1,q1)*qd1;                  % Second derivatives (dd)

%Scleronomic rotational hinge constraint Mass-Mass
C2      = (R(theta(1))*(H-com1)+[x(1);y(1)]) - (R(theta(2))*(H-com2)+[x(2);y(2)]);     %zeroth derivative
C2_q    = jacobian(C2,q2);                            % First derviatives (d)
C2_qq   = jacobian(C2_q*qd2,q2)*qd2;                  % Second derivative (dd)

%Scleronomic prismatic hinge constraint Ground-Mass
C3pend   = [x(1);y(1)]-(R(theta(1))*(com1-H)+H);       % Constraints of a pendulum
C3      = C3pend(1)*sin(alpha)-C3pend(2)*cos(alpha);    % Constraint of a pendulum on a track  
C3_q    = jacobian(C3,q1);                            % First derviatives (d)
C3_qq   = jacobian(C3_q*qd1,q1)*qd1;                  % Second derivatives (dd)
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       
%                       Rheonomic constraint equations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Rheonomic prismatic hinge with asin(wt) motor constraint

motorFunction = a*sin(w*t);
% motorFunction = a-exp(-w*t);

C4(1)       = C3;
C4(2)       = C3pend(1) + C3pend(2) - motorFunction;




C4_q        = jacobian(C4,q1);                            % First derviatives (d)
C4_qq       = jacobian(C4_q*qd1,q1)*qd1;                  % Second derivatives (dd)

C4_t        = jacobian(C4,t);
C4_qt       = jacobian(C4_t,q1)*qd1;
C4_tt       = jacobian(C4_t,t);

 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       
%                       Force equations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Spring Ground-Mass
Ls1 = ps1-(R(theta(1))*(ps2-com1)+[x(1);y(1)]);
Es1 = 0.5*K*(sqrt(Ls1.'*Ls1)-L0)^2;
Fs1 =  - jacobian(Es1,q1).';



%Spring Mass-Mass
Ls2 = (R(theta(1))*(ps1-com1)+[x(1);y(1)]) - (R(theta(2))*(ps2-com2)+[x(2);y(2)]);
Es2 = 0.5*K*(sqrt(Ls2.'*Ls2)-L0)^2;
Fs2 = - jacobian(Es2,q2).';

keyboard 

% keyboard

%Torsion Ground-Mass
Et1 = 0.5*K*(theta(1)-theta0)^2;
Ft1 = -jacobian(Et1,q1).';

%Torsion Mass-Mass
% Et2 = 0.5*K*(theta(2)-theta(1))^2;
Et2 = 0.5*K*((theta(2)-theta(1))-theta0)^2;
Ft2 = -jacobian(Et2,q2).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       
%                       Build matlabfunctions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Building matlab functions')

matlabFunction(C1,      'file','MGH','vars',{H,com1,q1});    
matlabFunction(C1_q,    'file','MGH_q','vars',{H,com1,q1});  
matlabFunction(C1_qq,   'file','MGH_qq','vars',{H,com1,q1,qd1});  

matlabFunction(C2,      'file','MMH','vars',{H,com1,com2,q2});  
matlabFunction(C2_q,    'file','MMH_q','vars',{H,com1,com2,q2});  
matlabFunction(C2_qq,   'file','MMH_qq','vars',{H,com1,com2,q2,qd2});  

matlabFunction(C3,      'file','MGP','vars',{H,alpha,com1,q1});  
matlabFunction(C3_q,    'file','MGP_q','vars',{H,alpha,com1,q1});  
matlabFunction(C3_qq,   'file','MGP_qq','vars',{H,alpha,com1,q1,qd1}); 

matlabFunction(C4,      'file','MGPm','vars',{H,alpha,a,w,com1,t,q1});
matlabFunction(C4_q,    'file','MGPm_q','vars',{H,alpha,a,w,com1,t,q1});
matlabFunction(C4_qq,   'file','MGPm_qq','vars',{H,alpha,a,w,com1,t,q1,qd1});

matlabFunction(C4_t,    'file','MGPm_t','vars',{H,alpha,a,w,com1,t,q1});
matlabFunction(C4_qt,   'file','MGPm_qt','vars',{H,alpha,a,w,com1,t,q1,qd1});
matlabFunction(C4_tt,   'file','MGPm_tt','vars',{H,alpha,a,w,com1,t,q1});

matlabFunction(Ls1,     'file','MGspringlength','vars',{ps1,ps2,L0,K,com1,q1});  
matlabFunction(Ls2,     'file','MMspringlength','vars',{ps1,ps2,L0,K,com1,com2,q2});  

matlabFunction(Fs1,     'file','MGspringforce','vars',{ps1,ps2,L0,K,com1,q1});  
matlabFunction(Fs2,     'file','MMspringforce','vars',{ps1,ps2,L0,K,com1,com2,q2});  

matlabFunction(Ft1,     'file','MGtorsion','vars',{K,theta0,q1});
matlabFunction(Ft2,     'file','MMtorsion','vars',{K,theta0,q2});

matlabFunction(Es1,     'file','MGspringenergy','vars',{ps1,ps2,L0,K,com1,q1});
matlabFunction(Es2,     'file','MMspringenergy','vars',{ps1,ps2,L0,K,com1,com2,q2});

matlabFunction(Et1,     'file','MGtorsionenergy','vars',{K,theta0,q1});
matlabFunction(Et2,     'file','MMtorsionenergy','vars',{K,theta0,q2});



disp('Constraints have been updated')











