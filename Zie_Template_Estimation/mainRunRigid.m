function [rigidReg,maskReg] = mainRunRigid(target, source,targetVoxelSize,sourceVoxelSize,mask)

target = (target-min(target(:)));
target = (target./max(target(:)))*130;



%[rigidParams] = rigidRegister_Bidir([],target,[],source,[], targetVoxelSize, sourceVoxelSize,[2], [3,3,3],zeros(12,1),1,'nmi');
%[rigidParams] = rigidRegister_Bidir([],target,[],source,[],  targetVoxelSize, sourceVoxelSize,[1], [2,2,2],rigidParams,0,'nmi');
%[rigidParams] = rigidRegister_Bidir([],target,[],source,[], targetVoxelSize, sourceVoxelSize,[0.2], [2,2,2],rigidParams,0,'nmi');


rigidParams = zeros(12,1);

corCor1=generateEvaluationPoints(sourceVoxelSize, size(source).*sourceVoxelSize, [1 1 1].*sourceVoxelSize);
[corCor]=double(applyRigidSym(rigidParams(:),corCor1,(size(source)/2).*sourceVoxelSize,0));
rigidReg = interpn(reshape(corCor1(:,1),size(source)),reshape(corCor1(:,2),size(source)),reshape(corCor1(:,3),size(source)),source,reshape(corCor(:,1),size(source)),reshape(corCor(:,2),size(source)),reshape(corCor(:,3),size(source)),'cubic',0);
rigidReg = reshape(rigidReg,size(source));


maskReg = interpn(reshape(corCor1(:,1),size(source)),reshape(corCor1(:,2),size(source)),reshape(corCor1(:,3),size(source)),mask,reshape(corCor(:,1),size(source)),reshape(corCor(:,2),size(source)),reshape(corCor(:,3),size(source)),'nearest',0);
maskReg = reshape(maskReg,size(source));

end
