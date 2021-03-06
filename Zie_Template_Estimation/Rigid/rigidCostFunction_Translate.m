function [f,df] = rigidCostFunction_Translate(p, I, J, points,centerI, centerJ, Itrival,Jtrival,scaleI,scaleJ,det,onlyRigid,direc,similarity,bins)

sizePoints = size(points,1);

df=zeros(numel(p),1);

t1=zeros(3,1);t1(1)=p(10);t1(2)=p(11);t1(3)=p(12);

fwpoints=points+repmat(t1',[sizePoints 1]);


switch similarity
    case 'nmi'
        [res,gradSim(:,1),gradSim(:,2),gradSim(:,3)]=NMI3D_DET_PW3(fwpoints,Jtrival+2,I+2,[0  bins-1 0  bins-1 ], [bins bins], [0 0 0],[scaleI(1) scaleI(2) scaleI(3)],det);
    case 'pnorm'
        [res,gradSim(:,1),gradSim(:,2),gradSim(:,3)]=PNorm(fwpoints,Jtrival+2,I+2,[0  bins-1 0  bins-1 ],[bins bins],[0 0 0],scaleI, 2,double(det),1);
        res = res/length(fwpoints);
        gradSim = gradSim/length(fwpoints);
end



gradSim(:,1)=gradSim(:,1)./scaleI(1);
gradSim(:,2)=gradSim(:,2)./scaleI(2);
gradSim(:,3)=gradSim(:,3)./scaleI(3);


df(10)=sum(gradSim(:,1));
df(11)=sum(gradSim(:,2));
df(12)=sum(gradSim(:,3));



switch similarity
    
    case 'nmi'
        f=6-res;
        
    case 'pnorm'
        f=res;
end

%T=(S*R)';

end
