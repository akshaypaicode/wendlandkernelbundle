function [f,df,T] = rigidCostFunction_SingleScale(p, I, J, points,centerI, centerJ, Itrival,Jtrival,scaleI,scaleJ,det,onlyRigid,direc,similarity,bins)

sizePoints = size(points,1);

df=zeros(numel(p),1);

rotx=p(1); roty=p(2);rotz=p(3);
t1=zeros(3,1);t1(1)=p(10);t1(2)=p(11);t1(3)=p(12);


%Rotation matrices
Rx=[1 0 0; 0 cos(rotx) -sin(rotx); 0 sin(rotx) cos(rotx)];
Ry=[cos(roty) 0 sin(roty); 0 1 0; -sin(roty) 0 cos(roty)];
Rz=[cos(rotz) -sin(rotz) 0; sin(rotz) cos(rotz) 0; 0 0 1];
R=Rz*Ry*Rx; 


%derived rotation matrices
dRx=Rz*Ry*[0 0 0; 0 -sin(rotx) -cos(rotx); 0 cos(rotx) -sin(rotx)];
dRy=Rz*[-sin(roty) 0 cos(roty); 0 0 0; -cos(roty) 0 -sin(roty)]*Rx;
dRz=[-sin(rotz) -cos(rotz) 0; cos(rotz) -sin(rotz) 0; 0 0 0]*Ry*Rx;

fwpoints=(points-repmat(centerI,sizePoints,1))*R+repmat(t1',[sizePoints 1])+repmat(centerI,sizePoints,1);


switch similarity
    case 'nmi'
        [res,gradSim(:,1),gradSim(:,2),gradSim(:,3)]=NMI3D_DET_PW3(fwpoints,Jtrival+2,I+2,[0  bins-1 0  bins-1 ], [bins bins], [0 0 0],[scaleI(1) scaleI(2) scaleI(3)],det);
    case 'pnorm'
        [res,gradSim(:,1),gradSim(:,2),gradSim(:,3)]=PNorm_det(fwpoints,Jtrival,I,[0  bins-1 0  bins-1 ],[140 140],[0 0 0],scaleI, 1,double(det));
end



gradSim(:,1)=gradSim(:,1)./scaleI(1);
gradSim(:,2)=gradSim(:,2)./scaleI(2);
gradSim(:,3)=gradSim(:,3)./scaleI(3);


df(10)=sum(gradSim(:,1));
df(11)=sum(gradSim(:,2));
df(12)=sum(gradSim(:,3));
%df(4)=sum(sum((points-repmat(centerI,sizePoints,1))*R.*gradSim,2))*~onlyRigid;
df(1)=sum(sum((points-repmat(centerI,sizePoints,1))*dRx.*gradSim,2));
df(2)=sum(sum((points-repmat(centerI,sizePoints,1))*dRy.*gradSim,2));
df(3)=sum(sum((points-repmat(centerI,sizePoints,1))*dRz.*gradSim,2));



switch similarity
    
    case 'nmi'
        f=6-res;
        
    case 'pnorm'
        f=res;
end

T=R;

end
