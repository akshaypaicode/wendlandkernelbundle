function [image, jacobDet] = updateTempJacobian(freeformParams, ffSpacing, ffSupportInterp, iter, image, voxelSize)

pts = generateEvaluationPoints(voxelSize, size(image).*voxelSize, [1 1 1].*voxelSize);



szP = [];
for k=1:iter
    freeformParamsNew{k}=-1.*freeformParams{k};
    sZ=size(freeformParamsNew{k});
    sZ=sZ(1:3);
    szP=vertcat(szP(:),sZ(:));
    numelP(k)=sZ(1)*sZ(2)*sZ(3);
end

clear freeformParams

jacobDet = zeros(length(pts),1);

imVec = zeros(length(pts),1);

parfor i=1:length(pts)
    
    [ptsNew,d1,d2,d3]=wendkernbuncompositionJacobian_tbb(pts(i,:),freeformParamsNew,ffSpacing,ffSupportInterp,5,szP,numelP,iter,iter,voxelSize,0);
    
    jac = jacobiandettransform(d1,d2,d3);
    
    jacobDet(i,1) = jac;
    
    [val(i,1)]=imagegradBSP(ptsNew, image, [0 0 0], voxelSize);%dtheta(v)
    
    
    
end

image = reshape(val.*jacobDet, size(image));

jacobDet = reshape(jacobDet, size(image));




end