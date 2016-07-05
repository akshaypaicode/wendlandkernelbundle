function [dlp] = DlogPosterior(p, pts, image, templateImage, ffSpacing,ffSupportInterp, ffSupportCov, iter, sigma, lambda, voxelSize, dimP, freeformParams)

numP = dimP(1)*dimP(2)*dimP(3)*dimP(4);

%lambda = lambda/numP;

ffSupportCov = round(ffSupportCov)*2;

N = numP/3;


dlp = [];

ssd = 0;

szP = []; 
for k=1:iter
    sZ=size(freeformParams{1}{k});
    sZ=sZ(1:3);
    szP=vertcat(szP(:),sZ(:));
    numelP(k)=sZ(1)*sZ(2)*sZ(3);
end


parfor j=1:length(image)
    
    st = (numP*(j-1))+1;
    en = ((numP*(j-1))+1)+(numP-1);
    
    freeformParams{j}{end} = reshape(p(st:en),dimP);
    
    tic
    %for m = 1:length(pts)
    
    
    [fwpoints,idx,tt]=wendkernbuncompositionparallel_tbb(pts,freeformParams{j},ffSpacing,ffSupportInterp,5,szP,numelP,iter,iter,voxelSize,0);
  
    [val,gradSimx,gradSimy,gradSimz]=imagegradBSP(fwpoints, templateImage, [0 0 0], voxelSize);%dtheta(v)
    
    
    d = (val(:)-image{j}(:));
    
    
    %end
    toc
    
    temp = tt./5;
    
    d = repmat(d,[1 3]);
    vZ = repmat(voxelSize,[length(pts) 1]);
    d = d./vZ;
    
    
    [~,~,C] = wendsquarednorm(freeformParams{j}{end},ffSpacing,ffSupportCov, 1);%for covariance, 1 is dummy
    
    C = sparse(blkdiag(C,C,C));
    
    C = lambda*C;
    
    %logdetC = logdet(C,'chol');
    
    %regLambda = (freeformParams{j}{end}(:)')*(C*freeformParams{j}{end}(:));
    
    dReg = (freeformParams{j}{end}(:)'*(C+C'));
    
    
    df = dDdPFuncJac([gradSimx gradSimy gradSimz].*d,squeeze(temp(1:(ffSupportInterp(end)^3),1,:)),squeeze(temp(1:(ffSupportInterp(end)^3),2,:)),squeeze(temp(1:(ffSupportInterp(end)^3),3,:)),...
        squeeze(temp((ffSupportInterp(end)^3+1):(ffSupportInterp(end)^3*2),1,:)),squeeze(temp((ffSupportInterp(end)^3+1):(ffSupportInterp(end)^3*2),2,:)),squeeze(temp((ffSupportInterp(end)^3+1):(ffSupportInterp(end)^3*2),3,:)),...
        squeeze(temp((ffSupportInterp(end)^3*2)+1:(ffSupportInterp(end)^3*3),1,:)),squeeze(temp((ffSupportInterp(end)^3*2)+1:(ffSupportInterp(end)^3*3),2,:)),squeeze(temp((ffSupportInterp(end)^3*2)+1:(ffSupportInterp(end)^3*3),3,:)), idx, N,ffSupportInterp(end));
    
    
     %a = -1*logdetC;
    
    
    %lp =  (0.5)*a - (0.5*regLambda) - (numel(image{j})*log(sigma)) - ((0.5*sigma^-2)*ssd);
    
    dlp = vertcat(dlp,((-0.5*dReg(:)) - (sigma^-2)*df(:))); %dlp = dlp*10;

end

%(0.5)*log(1./detC) - (0.5*lambda*reg) - (numel(image)*log(sigma)) - ((0.5*sigma^-2)*ssd);

dlp = -1*dlp';
