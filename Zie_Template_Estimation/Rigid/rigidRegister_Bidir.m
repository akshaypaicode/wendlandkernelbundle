function [rigidParams, matchPenalty] = rigidRegister_Bidir(direc,target,targetAseg, source, sourceAseg,targetVoxelSize, sourceVoxelSize, smoothing,samplingDensity, rigidParams,onlyTranslate,similarity)
global debugGraph

% function rigidParams = rigidRegister(target, source, targetVoxelSize, sourceVoxelSize, optO)
%
% Performs multiscale rigid image registration of the input images.
%
% INPUT:
% - target         : Target series.
% - source         : Source series.
% - targetVoxelSize: The voxel size of the target ([x,y,z]).
% - sourceVoxelSize: The voxel size of the source ([x,y,z]).
% - optO           : Options objects.
%
% OUTPUT:
% - rigidParams:  The rigid transformation parameters (applied to source they transform into target).
% - matchPenalty: Final NMI match penalty

% Get parameters from the options object
maxIter            = 100;
tolFunExp          = -8;
tolXExp            = -8;
bins               = 130;

sourceVoxelSize = double(sourceVoxelSize(1:3));
targetVoxelSize = double(targetVoxelSize(1:3));


% Set optimization options
warning off MATLAB:nearlysingularmatrix
options.method          = 'lbfgs';
options.maxIter         = maxIter;
options.maxFunEvals     = 150;
options.tolFun          = 10^(tolFunExp);
options.tolX            = 10^(tolXExp);
options.derivativeCheck = 'off';
%  options.display = 'off';
options.LS_interp = 1; % Simpler line-search interpolation than default cubic, but cross-platform-robust

% Options names have changed on new minFunc 2012
options.optTol = options.tolX;
options.progTol = options.tolFun;

% Initialize rigid parameters (if not supplied in function call)


% Get rescaled dimensions of target and source
targetDim = size(target);
sourceDim = size(source);

% computeAffineMatricesandDerivatives(direc);

ptsInitial = generateEvaluationPoints(targetVoxelSize, targetDim.*targetVoxelSize, samplingDensity.*targetVoxelSize); % original

for i=1:length(smoothing)
    
    
    
    % Perform smoothing and rescale intensities according to desired bin size.
    s = smoothing(i)./double(targetVoxelSize);
    
    
    if s>0
        
        targetPreP = real(ifftn(scalen(fftn(target), s))); targetPreP = targetPreP - min(targetPreP(:));
        sourcePreP = real(ifftn(scalen(fftn(source), s))); sourcePreP = sourcePreP - min(sourcePreP(:));
        
    else
        
        targetPreP = target;
        sourcePreP = source;
        
    end
    
    targetPreP = targetPreP / max(targetPreP(:));
    sourcePreP = sourcePreP / max(sourcePreP(:));
    
    if( strcmp(similarity, 'nmi') )
        targetPreP = targetPreP * bins;
        sourcePreP = sourcePreP * bins;
    end
    
    Jtrivial = SplineInterpolation(double(ptsInitial), targetPreP, [0 0 0], double(targetVoxelSize));%interpolationWrapper('cubic', ptsInitial, double(targetPreP), double(targetVoxelSize), 10000);
    Itrivial = SplineInterpolation(double(ptsInitial), sourcePreP, [0 0 0], double(sourceVoxelSize));%interpolationWrapper('cubic', ptsInitial, double(sourcePreP), double(sourceVoxelSize), 10000);
    
    p = double([rigidParams(:)]);
    
    centerI = double([sourceDim(1)/2*sourceVoxelSize(1) sourceDim(2)/2*sourceVoxelSize(2) sourceDim(3)/2*sourceVoxelSize(3)]);
    centerJ = double([targetDim(1)/2*targetVoxelSize(1) targetDim(2)/2*targetVoxelSize(2) targetDim(3)/2*targetVoxelSize(3)]);
    
    %only translate
    
    if onlyTranslate
        
        disp('only translating')
        
        [p]=minFunc(@rigidCostFunction_Translate,p,options, double(sourcePreP), double(targetPreP), double(ptsInitial),centerI,centerJ, double(Itrivial),double(Jtrivial),double(sourceVoxelSize),double(targetVoxelSize),double(ones(size(ptsInitial,1),1)),0,direc,similarity,bins+5);
    else
        
        [p]=minFunc(@rigidCostFunction_Sym,p,options, double(sourcePreP), double(targetPreP), double(ptsInitial),centerI,centerJ, double(Itrivial),double(Jtrivial),double(sourceVoxelSize),double(targetVoxelSize),double(ones(size(ptsInitial,1),1)),0,direc,similarity,bins+5);
        
        
    end
    rigidParams = p
end

end
