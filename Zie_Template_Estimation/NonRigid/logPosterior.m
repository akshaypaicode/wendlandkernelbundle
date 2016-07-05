function [lp] = logPosterior(p, pts, image, templateImage, ffSpacing,ffSupportInterp, ffSupportCov, iter, sigma, lambda, voxelSize, dimP, freeformParams)

numP = dimP(1)*dimP(2)*dimP(3)*dimP(4);

%lambda = lambda/numP;

ffSupportCov = ffSupportCov*2;

dlpLamda = 0;

szP = [];
for k=1:iter
    sZ=size(freeformParams{1}{k});
    sZ=sZ(1:3);
    szP=vertcat(szP(:),sZ(:));
    numelP(k)=sZ(1)*sZ(2)*sZ(3);
end

lp = zeros(length(image),1);

parfor j=1:length(image)
    
    st = (numP*(j-1))+1;
    en = ((numP*(j-1))+1)+(numP-1);
    
    freeformParams{j}{end} = reshape(p(st:en),dimP);

    [fwpoints]=wendkernbuncompositioninterpolation(pts,freeformParams{j},ffSpacing,ffSupportInterp,5,szP,numelP,iter,iter,voxelSize,0);
    
    [val]=imagegradBSP(fwpoints, templateImage, [0 0 0], voxelSize);%dtheta(v)
    
    ssd = sum((val(:) - image{j}(:)).^2);
    
    [~,~,C] = wendsquarednorm(freeformParams{j}{end},ffSpacing,ffSupportCov, 1);
    
    C = sparse(blkdiag(C,C,C));
    
    C = lambda*C;
    
    logdetC = logdet(C,'chol');
    
    [regLambda] = freeformParams{j}{end}(:)'*C*freeformParams{j}{end}(:);
    
    a = -1*logdetC;
    
    
    lp(j,1)=(0.5)*a - (0.5*regLambda) - (numel(image{j})*log(sigma)) - ((0.5*sigma^-2)*ssd);
    
    %dlpLamda = dlpLamda + (0.5*(-1/detC)*(lambda^(length(C)-1))*detC) - (freeformParams{j}{end}(:)')*((C/lambda)*freeformParams{j}{end}(:));
    
    
    
end

lp = -1*sum(lp)