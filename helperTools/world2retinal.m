function theta = world2retinal(vel,pos,eyePos)

% Transform world speed to retinal speed

if isrow(vel)
    vel = vel';
end

xL = eyePos(1);
xR = eyePos(2);

x0 = pos(1);
y0 = pos(2);
z0 = pos(3);

hL = sqrt(z0^2 + (x0 - xL)^2);
hR = sqrt(z0^2 + (x0 - xR)^2);

a  = xR - xL;

tXL = [-z0/(hL^2) 0 (x0-xL)/(hL^2)];
tXR = [-z0/(hR^2) 0 (x0-xR)/(hR^2)];
tYL = [-(y0*(x0-xL))/(hL*(hL^2 + y0^2)) hL/(hL^2 + y0^2) -(y0*z0)/(hL*(hL^2 + y0^2))];
tYR = [-(y0*(x0-xR))/(hR*(hR^2 + y0^2)) hR/(hR^2 + y0^2) -(y0*z0)/(hR*(hR^2 + y0^2))];

theta = [tXL;tXR;tYL;tYR]*vel;

% ATrue = pinv([tXL;tXR;tYL;tYR]);

% a1 = [(y0^2)*hL*z0*((x0-xR)*hL - (x0-xL)*(hL^2 + y0^2))/(a*hR*(hR^2 + y0^2)*(hL^2 + y0^2)^2),...
%       (y0^2)*hR*z0*((x0-xR)*hL - (x0-xL)*(hL^2 + y0^2))/(a*hL*(hR^2 + y0^2)*(hL^2 + y0^2)^2),...
%       hL*(hL^2 + y0^2)*(x0-xR)/((x0-xR)*(hL^2) - (x0-xL)*(hR^2)),...
%       -hR*(hR^2 + y0^2)*(x0-xL)/((x0-xR)*(hL^2) - (x0-xL)*(hR^2))];
% 
% A = [a1;...
%      ((x0-xR)*hL^2)/(z0*a)              -((x0-xL)*hR^2)/(z0*a)                    0                0;...
%      (hL^2)/a                           -(hR^2)/a                                 0                0];

% A = [((x0-xR)*hL^2)/(z0*a)              -((x0-xL)*hR^2)/(z0*a)                    0                0;...
%      y0*(z0^2 + (x0-xL)*(x0-xR))/(a*z0) -(y0*hR^2)*(z0^2 - (x0-xL)^2)/(a*z0*hL^2) (hL^2 + y0^2)/hL 0;...
%      (hL^2)/a                           -(hR^2)/a                                 0                0];
% 
% v1 = ATrue*theta;
% v2 = A*theta;
% 
% keyboard

end