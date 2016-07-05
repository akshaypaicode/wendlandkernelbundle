function p = runEMTemplate(image, pts, templateImage, freeformParams, ffSpacing,ffSupportInterp, ffSupportCov, iter, sigma, lambda, voxelSize, origImage, numSamples)

warning off MATLAB:nearlysingularmatrix % minFunc shouts uninterestingly
optionsPred.method          = 'sd';
optionsPred.tolFun          = 10^(-8);
optionsPred.tolX            = 10^(-8);
optionsPred.derivativeCheck = 'off';
optionsPred.LS_interp       = 1;  % Simpler line-search interpolation than default cubic, but cross-platform-robust
optionsPred.optTol          = optionsPred.tolX;
optionsPred.progTol         = optionsPred.tolFun;
optionsPred.corrections     = 30; % Limit number of gradients stored in minFunc
optionsPred.MaxIter         = 12;
optionsPred.MaxFunEvals     = 12;
optionsPred.numDiff         = 1;



p = lambda;

p = minFunc(@EMTemplate, p, optionsPred, image, pts, templateImage, freeformParams, ffSpacing,ffSupportInterp, ffSupportCov, iter, sigma, lambda, voxelSize, origImage, numSamples);




end




