function [f,df] = rigidCostFunction_Sym(p, I, J, points,centerI, centerJ, Itrival,Jtrival,scaleI,scaleJ,det,onlyTranslate,direc,similarity,bins)

sizePoints = size(points,1);

df = zeros(length(p),1);

[fwpoints,dR]=applyRigidSym(p,points,centerI,0);
[invpoints,~,idR,invR]=applyRigidSym(p,points,centerJ,1);


t=zeros(3,1);t(1)=p(10);t(2)=p(11);t(3)=p(12);
centeredPtsI = points - repmat(centerI,sizePoints,1);
centeredPtsJ = points - repmat(centerJ,sizePoints,1) - repmat(t',[sizePoints 1]);


switch similarity
    case 'nmi'
        [res,gradSim(:,1),gradSim(:,2),gradSim(:,3)]=NMI3D_DET_PW3(fwpoints,Jtrival+2,I+2,[0  bins-1 0  bins-1 ], [bins bins], [0 0 0],[scaleI(1) scaleI(2) scaleI(3)],det);
        [invres,invd(:,1),invd(:,2),invd(:,3)]=NMI3D_DET_PW3(invpoints,Itrival+2,J+2,[0  bins-1 0  bins-1 ], [bins bins], [0 0 0],[scaleJ(1) scaleJ(2) scaleJ(3)],det);
    case 'pnorm'
        
        disp('pnorm');
        [res,gradSim(:,1),gradSim(:,2),gradSim(:,3)]=PNorm(fwpoints,Jtrival+2,I+2,[0  bins-1 0  bins-1 ],[bins bins],[0 0 0],scaleI, 2,double(det),1);
        [invres,invd(:,1),invd(:,2),invd(:,3)]=PNorm(invpoints,Itrival+2,J+2,[0  bins-1 0  bins-1 ],[bins bins],[0 0 0],scaleJ, 2,double(det),1);
end



gradSim(:,1)=gradSim(:,1)./scaleI(1);
gradSim(:,2)=gradSim(:,2)./scaleI(2);
gradSim(:,3)=gradSim(:,3)./scaleI(3);


invd(:,1)=invd(:,1)./scaleJ(1);
invd(:,2)=invd(:,2)./scaleJ(2);
invd(:,3)=invd(:,3)./scaleJ(3);


    
df(1) = sum(sum(centeredPtsI*dR{1}.*gradSim,2) + sum(centeredPtsJ*idR{1}.*invd,2));
df(2) = sum(sum(centeredPtsI*dR{2}.*gradSim,2) + sum(centeredPtsJ*idR{2}.*invd,2));
df(3) = sum(sum(centeredPtsI*dR{3}.*gradSim,2) + sum(centeredPtsJ*idR{3}.*invd,2));
df(4) = sum(sum(centeredPtsI*dR{4}.*gradSim,2) + sum(centeredPtsJ*idR{4}.*invd,2));
df(5) = sum(sum(centeredPtsI*dR{5}.*gradSim,2) + sum(centeredPtsJ*idR{5}.*invd,2));
df(6) = sum(sum(centeredPtsI*dR{6}.*gradSim,2) + sum(centeredPtsJ*idR{6}.*invd,2));
df(7) = sum(sum(centeredPtsI*dR{7}.*gradSim,2) + sum(centeredPtsJ*idR{7}.*invd,2));
df(8) = sum(sum(centeredPtsI*dR{8}.*gradSim,2) + sum(centeredPtsJ*idR{8}.*invd,2));
df(9) = sum(sum(centeredPtsI*dR{9}.*gradSim,2) + sum(centeredPtsJ*idR{9}.*invd,2));


id = invd*invR;

df(10)=sum(gradSim(:,1))-sum(id(:,1));
df(11)=sum(gradSim(:,2))-sum(id(:,2));
df(12)=sum(gradSim(:,3))-sum(id(:,3));





switch similarity
    
    case 'nmi'
        f=6-res+6-invres;
        
    case 'pnorm'
        f=res+invres;
end

end
