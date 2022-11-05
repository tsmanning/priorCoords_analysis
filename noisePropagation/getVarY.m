function sigmaY = getVarY(pos,eyePos,sigN)

% Calculate variance in world coordinates (y dimension) after propogating retinal noise 

% Define eye positions in x
xL = eyePos(1);
xR = eyePos(2);

sigmaY = nan(size(pos,1),1);

% Loop over input stimulus positions (probably could vectorize this with some effort)
for ii = 1:size(pos,1)

    % Stimulus position
    x0 = pos(ii,1);
    y0 = pos(ii,2);
    z0 = pos(ii,3);

    % Distance of stimulus from each eye
    hL = sqrt(z0^2 + (x0 - xL)^2);
    hR = sqrt(z0^2 + (x0 - xR)^2);

    % IPD
    a  = xR - xL;

    % World motion to retinal motion transform
    tXL = [-z0/(hL^2) 0 (x0-xL)/(hL^2)];
    tXR = [-z0/(hR^2) 0 (x0-xR)/(hR^2)];
    tYL = [-(y0*(x0-xL))/(hL*(hL^2 + y0^2)) hL/(hL^2 + y0^2) -(y0*z0)/(hL*(hL^2 + y0^2))];
    tYR = [-(y0*(x0-xR))/(hR*(hR^2 + y0^2)) hR/(hR^2 + y0^2) -(y0*z0)/(hR*(hR^2 + y0^2))];

    % Get retinal to world transform (Cheat with pseudoinverse/numerical methods)
    A = pinv([tXL;tXR;tYL;tYR]);

    % (Don't cheat this time and find it arithmetically)
    % 1st derivation
    % A = [y0*(z0^2 + (x0-xL)*(x0-xR))/(a*z0) -(y0*hR^2)*(z0^2 - (x0-xL)^2)/(a*z0*hL^2) (hL^2 + y0^2)/hL 0;...
    %      ((x0-xR)*hL^2)/(z0*a)              -((x0-xL)*hR^2)/(z0*a)                    0                0;...
    %      (hL^2)/a                           -(hR^2)/a                                 0                0];

    % 2nd derivation (wrong)
    % a1 = [(y0^2)*hL*z0*((x0-xR)*hL - (x0-xL)*(hL^2 + y0^2))/(a*hR*(hR^2 + y0^2)*(hL^2 + y0^2)^2),...
    %       (y0^2)*hR*z0*((x0-xR)*hL - (x0-xL)*(hL^2 + y0^2))/(a*hL*(hR^2 + y0^2)*(hL^2 + y0^2)^2),...
    %       hL*(hL^2 + y0^2)*(x0-xR)/((x0-xR)*(hL^2) - (x0-xL)*(hR^2)),...
    %       -hR*(hR^2 + y0^2)*(x0-xL)/((x0-xR)*(hL^2) - (x0-xL)*(hR^2))];
    %
    % A = [a1;...
    %      ((x0-xR)*hL^2)/(z0*a)              -((x0-xL)*hR^2)/(z0*a)                    0                0;...
    %      (hL^2)/a                           -(hR^2)/a                                 0                0];

    % Find the covariance matrix
    M = A*A'*sigN;

    % Isolate y component of noise after propagating retinal noise into
    % world space
    sigmaY(ii) = M(1,1) - M(1,2:3)/M(2:3,2:3)*M(2:3,1);

end

end