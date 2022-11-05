function sigmaX = getVarX(pos,eyePos,sigN)

% Calculate variance in world coordinates (y dimension) after propogating retinal noise 

xL = eyePos(1);
xR = eyePos(2);

x0 = pos(1);
y0 = pos(2);
z0 = pos(3);

hL = sqrt(z0^2 + (x0 - xL)^2);
hR = sqrt(z0^2 + (x0 - xR)^2);

a  = xR - xL;

% A = [((x0-xR)*hL^2)/(z0*a)              -((x0-xL)*hR^2)/(z0*a)                    0                0;...
%      y0*(z0^2 + (x0-xL)*(x0-xR))/(a*z0) -(y0*hR^2)*(z0^2 - (x0-xL)^2)/(a*z0*hL^2) (hL^2 + y0^2)/hL 0;...
%      (hL^2)/a                           -(hR^2)/a                                 0                0];

a1 = [(y0^2)*hL*z0*((x0-xR)*hL - (x0-xL)*(hL^2 + y0^2))/(a*hR*(hR^2 + y0^2)*(hL^2 + y0^2)^2),...
      (y0^2)*hR*z0*((x0-xR)*hL - (x0-xL)*(hL^2 + y0^2))/(a*hL*(hR^2 + y0^2)*(hL^2 + y0^2)^2),...
      hL*(hL^2 + y0^2)*(x0-xR)/((x0-xR)*(hL^2) - (x0-xL)*(hR^2)),...
      -hR*(hR^2 + y0^2)*(x0-xL)/((x0-xR)*(hL^2) - (x0-xL)*(hR^2))];

A = [((x0-xR)*hL^2)/(z0*a)              -((x0-xL)*hR^2)/(z0*a)                    0                0;...
     a1;...
     (hL^2)/a                           -(hR^2)/a                                 0                0];

M = A*A'*sigN;

sigmaX = M(1,1) - M(1,2:3)/M(2:3,2:3)*M(2:3,1);

end