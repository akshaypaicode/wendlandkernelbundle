function [f,df] = rigidCostFunction_Sym_uni(p, I, J, points,centerI, centerJ, Itrival,Jtrival,scaleI,scaleJ,det,onlyTranslate,direc,similarity,bins)

sizePoints = size(points,1);

df = zeros(length(p),1);

[fwpoints,dR]=applyRigidSym(p,points,centerI,0);


t=zeros(3,1);t(1)=p(10);t(2)=p(11);t(3)=p(12);
centeredPtsI = points - repmat(centerI,sizePoints,1);


switch similarity
    case 'nmi'
        [res,gradSim(:,1),gradSim(:,2),gradSim(:,3)]=NMI3D_DET_PW3(fwpoints,Jtrival+2,I+2,[0  bins-1 0  bins-1 ], [bins bins], [0 0 0],[scaleI(1) scaleI(2) scaleI(3)],det);
    case 'pnorm'
        [res,gradSim(:,1),gradSim(:,2),gradSim(:,3)]=PNorm_det(fwpoints,Jtrival,I,[0 130 0 130],[140 140],[0 0 0],scaleI, 1,double(det));
end



gradSim(:,1)=gradSim(:,1)./scaleI(1);
gradSim(:,2)=gradSim(:,2)./scaleI(2);
gradSim(:,3)=gradSim(:,3)./scaleI(3);

    
df(1) = sum(sum(centeredPtsI*dR{1}.*gradSim,2));
df(2) = sum(sum(centeredPtsI*dR{2}.*gradSim,2));
df(3) = sum(sum(centeredPtsI*dR{3}.*gradSim,2));
df(4) = sum(sum(centeredPtsI*dR{4}.*gradSim,2));
df(5) = sum(sum(centeredPtsI*dR{5}.*gradSim,2));
df(6) = sum(sum(centeredPtsI*dR{6}.*gradSim,2));
df(7) = sum(sum(centeredPtsI*dR{7}.*gradSim,2));
df(8) = sum(sum(centeredPtsI*dR{8}.*gradSim,2));
df(9) = sum(sum(centeredPtsI*dR{9}.*gradSim,2));

df(10)=sum(gradSim(:,1));
df(11)=sum(gradSim(:,2));
df(12)=sum(gradSim(:,3));





switch similarity
    
    case 'nmi'
        f=6-res;
        
    case 'pnorm'
        f=res;
end

end
