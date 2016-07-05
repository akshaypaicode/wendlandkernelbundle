function [lp] = logPosteriorLambda(p, pts, image, templateImage, ffSpacing,ffSupportInterp, ffSupportCov, iter, sigma, voxelSize, dimP, freeformParams,numSamples, ssd)

lambda = p

numP = dimP(1)*dimP(2)*dimP(3)*dimP(4);

%lambda = lambda/numP;

ffSupportCov = ffSupportCov*2;

ssd = 0;

lp = [];

szP = [];
for k=1:iter
    sZ=size(freeformParams{1}{1}{k});
    sZ=sZ(1:3);
    szP=vertcat(szP(:),sZ(:));
    numelP(k)=sZ(1)*sZ(2)*sZ(3);
end

if lambda<=0
    
    lp = NaN;
    return
end


        

for s=1:numSamples
    
    for j=1:length(image)
        
        [fwpoints]=wendkernbuncompositioninterpolation(pts,freeformParams{s}{j},ffSpacing,ffSupportInterp,5,szP,numelP,iter,iter,voxelSize,0);
        
        [val]=imagegradBSP(fwpoints, templateImage, [0 0 0], voxelSize);%dtheta(v)

        ssd = (val - image{j}(:)).^2;

        
        [regLambda,~,C] = wendsquarednorm(freeformParams{s}{j}{end},ffSpacing,ffSupportCov, 1);

        
        C = lambda*C;
        
        C = sparse(blkdiag(C,C,C));
        
        logdetC = logdet(C,'chol');
        
        
        
        a = -1*logdetC;
        
        
        lp{s}(j) = (0.5)*a - (0.5*regLambda) - (numel(image{j})*log(sigma)) - ((0.5*sigma^-2)*ssd);
        
    end
    
end

lambda

lp = cell2mat(lp);

lp = sum(lp(:));

lp = 10^10-(lp/numSamples)
