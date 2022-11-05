function vel = retinal2world(theta,pos,eyePos)

% Transform retinal speed to world speed

tXL = theta(1);
tXR = theta(2);
tYL = theta(3);
tYR = theta(4);

xL = eyePos(1);
xR = eyePos(2);

x0 = pos(1);
y0 = pos(2);
z0 = pos(3);

hL = sqrt(z0^2 + (x0 - xL)^2);
hR = sqrt(z0^2 + (x0 - xR)^2);

a  = xR - xL;

dx = 1/(z0*a) * ((x0-xR)*tXL*hL^2 - (x0-xL)*tXR*hR^2);
dy = tYL*(hL^2 +y0^2)/hL + ...
    (tXL*hL^2 - tXR*hR^2)*(y0*z0)/(a*hL^2) + ...
    ((x0-xR)*tXL*hL^2 - (x0-xL)*tXR*hR^2)*(y0*(x0-xL))/(z0*a*hL^2);
dz = 1/a * (tXL*hL^2 - tXR*hR^2);

vel = [dx;dy;dz];

end