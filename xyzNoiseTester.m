clear all
close all

% This is how we calc screen vel in m/s
vDeg    = 4;
hMet    = 0.5967;

dva     = 2*atand(hMet);
thisVel = vDeg*(hMet/dva);

velY = [0 thisVel 0]';
velX = [thisVel 0 0]';

%%%%%%%%%%%%%
posNear = [-0.134 0.0526 0.5];
posFar  = [ 0.268 0.1051 1.0];

posNearSimp = [0 0 0.5];
posFarSimp  = [0 0 1.0];

eyePos  = [-0.03 0.03];

thetaNearX = world2retinal(velX,posNear,eyePos);
thetaFarX  = world2retinal(velX,posFar,eyePos);
thetaNearY = world2retinal(velY,posNear,eyePos);
thetaFarY  = world2retinal(velY,posFar,eyePos);
% theta     = [0 0 4 4];

sigN      = 2;

%%%%%
sigmaYNear = getVarY(posNear,eyePos,sigN);
sigmaYFar  = getVarY(posFar,eyePos,sigN);

sigmaYNearSimp = getVarY(posNearSimp,eyePos,sigN);
sigmaYFarSimp  = getVarY(posFarSimp,eyePos,sigN);

ratY     = sigmaYFar/sigmaYNear;
ratSimpY = sigmaYFarSimp/sigmaYNearSimp;

%%%%%
sigmaXNear = getVarX(posNear,eyePos,sigN);
sigmaXFar  = getVarX(posFar,eyePos,sigN);

sigmaXNearSimp = getVarX(posNearSimp,eyePos,sigN);
sigmaXFarSimp  = getVarX(posFarSimp,eyePos,sigN);

ratX     = sigmaXFar/sigmaXNear;
ratSimpX = sigmaXFarSimp/sigmaXNearSimp;

figure;
hold on;
bar([ratY ratSimpY ratX ratSimpX]);
plot([2.5 2.5],[0 5],'k');
set(gca,'fontsize',15,'xtick',1:4,'xticklabel',{'y','y (Simp.)','x','x (Simp.)'});
ylabel('ratio \sigma^{2} far/near');