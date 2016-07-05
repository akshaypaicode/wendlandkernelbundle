function computeAffineMatricesandDerivatives(direc)

syms rotX rotY rotZ sX sY sZ shX shY shZ;

rX = [1 0 0; 0 cos(rotX) -sin(rotX); 0 sin(rotX) cos(rotX)];
rY = [cos(rotY) 0 sin(rotY); 0 1 0; -sin(rotY) 0 cos(rotY)];
rZ = [cos(rotZ) -sin(rotZ) 0; sin(rotZ) cos(rotZ) 0; 0 0 1];

r = rX*rY*rZ;

s = diag([1+sX,1+sY,1+sZ]);

shxy = [1 0 shX; 0 1 shY; 0 0 1];
shyx = [1 shX 0; 0 1 0; 0 shZ 1];
shyz = [1 0 0; shY 1 0; shZ 0 1];


affR = s*r*shxy*shyx*shyz;

drotX = diff(affR, rotX);
drotY = diff(affR, rotY);
drotZ = diff(affR, rotZ);

dsX = diff(affR, sX);
dsY = diff(affR, sY);
dsZ = diff(affR, sZ);

dshX = diff(affR, shX);
dshY = diff(affR, shY);
dshZ = diff(affR, shZ);


invaffR = inv(affR);

dirotX = diff(invaffR, rotX);
dirotY = diff(invaffR, rotY);
dirotZ = diff(invaffR, rotZ);

disX = diff(invaffR, sX);
disY = diff(invaffR, sY);
disZ = diff(invaffR, sZ);

dishX = diff(invaffR, shX);
dishY = diff(invaffR, shY);
dishZ = diff(invaffR, shZ);





save(strcat(direc,'/affine.mat'));
