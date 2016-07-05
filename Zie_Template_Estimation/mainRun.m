clc
clear all

dataset='MGH10';
seq='g';

sourcedir=strcat('/home/pai/akshay/Sune_Data/');

%optionsText = '/home/pai/akshay/Source_Wend_KernBun_Template/';
addpath_recurse('/home/pai/Zie_Template_Estimation');
%addpath_recurse('/Users/akshay/Desktop/mountpoint/ReadData3D_version1k/');
%addpath_recurse('/Users/akshay/Desktop/mountpoint/Zie_Template_Estimation/freesurfer_matlab/');
p=pwd;
cd('/home/pai/Zie_Template_Estimation/C++/');
mex -O dDdPFuncJac.cpp
mex -O imagegradBSP.cpp
mex -O SplineInterpolation.cpp
mex -O wendkernbuncompositionparallel.cpp
mex -O wendkernbuncompositionJacobian_tbb.cpp -ltbb -ltbbmalloc
mex -O wendkernbuncompositionparallel_tbb.cpp -ltbb -ltbbmalloc
mex -O wendsquarednorm.cpp
mex -O wendkernbuncompositioninterpolation.cpp
mex -O jacobiandettransform.cpp
mex -O NMI3D_DET_PW3.cpp
mex -O wendkernbuncompositionJacobian.cpp
cd(p)


voxelSize = [1.00 1.33 1.00];

for i=1:10
    
    imt=strcat(sourcedir,'/EBrains/',seq,num2str(i),'.nii');
    at=strcat(sourcedir,'/Atlases/',seq,num2str(i),'lr.nii.gz');
    [~,~,ext] = fileparts(imt);
    if strcmp(ext,'.gz') ext = '.nii.gz'; end;
    [image{i}] = loadImage(imt,ext);
    [atlas{i}] = loadImage(at,ext);
    imageDimensions = size(image);
    
end

image{1} = image{1}-min(image{1}(:));
image{1} = (image{1}./max(image{1}(:)))*1;



meanImage = image{1};
atlases = atlas{1};
%image{1} = fitImage(image{1},voxelSize);
source = image{1};

parfor i = 2:length(image)
    
    
    image{i} = image{i}-min(image{i}(:));
    image{i} = (image{i}./max(image{i}(:)))*1;
    %image{i} = fitImage(image{i},voxelSize);
    target = image{i};
    [imReg, maskReg] = mainRunRigid(source,target,voxelSize,voxelSize,atlas{i});
    image{i} = imReg;
    atlas{i} = maskReg;
    image{i}(isnan(image{i}(:)))=0;
    meanImage = meanImage + imReg;
    atlases = atlas{i} + atlases;
    
end

clear imReg
meanImage = meanImage/length(image);

points = generateEvaluationPoints(voxelSize, size(meanImage), [10 10 10].*voxelSize);

points=reduceeval(points,atlases~=0,6,voxelSize);


ssd = 0;
for i=1:length(image)
    images{i}=SplineInterpolation(double(points), double(image{i}), [0 0 0], double(voxelSize)); %pts is not transformed to begin with
    temp=SplineInterpolation(double(points), double(meanImage), [0 0 0], double(voxelSize)); %pts is not transformed to begin wit
    
    ssd = ssd + sum((temp(:) - images{i}(:)).^2);
end

ssd = ssd./(length(images{1}(:))*length(images));


ffSpacing = [20 10 6];
ffSupportInterp = [2 2 2];
ffSupportCov = [3];
ffSpacingI=[];
imm = image{1};

for i=1:length(ffSpacing)
    ffTemp = [ffSpacing(i) ffSpacing(i) ffSpacing(i)].*voxelSize;
    ffSpacingI = vertcat(ffSpacingI(:),ffTemp(:));
end

numSamples=15;


for i=1:length(ffSpacing)
    
    %Estimating the covariance matrix
    
    ffSpacingInput = ffSpacingI(1:(3*(i-1))+3);
    
    
    for k=1:length(image)
        

        for l=1:numSamples
            [x,~,~]=ndgrid(voxelSize(1):ffSpacingInput(end-2):size(imm,1)*voxelSize(1),voxelSize(2):ffSpacingInput(end-1):size(imm,2)*voxelSize(2),voxelSize(3):ffSpacingInput(end):size(imm,3)*voxelSize(3));
            ffNumControlPoints=size(x)+2;clear x y z;
            freeformParams{l}{k}{i,1} = zeros([ffNumControlPoints 3]);
            size(freeformParams{l}{k}{i,1})
            
        end
    end
%     
    poolobj = gcp('nocreate');
    delete(poolobj);
    myCluster=parcluster('local'); myCluster.NumWorkers=22; parpool(myCluster,22);
    
    delete('log.txt');
    diary('log.txt');
    
    [meanImage,freeformParams] = EMTemplate(images, points, meanImage, freeformParams, ffSpacingInput,ffSupportInterp(1:i), ffSupportCov, i, sqrt(ssd), 1, voxelSize, image,numSamples);
    
    diary('off');
end







