function [templateImage,freeformParams] = EMTemplate(image, pts, templateImage, freeformParams, ffSpacing,ffSupportInterp, ffSupportCov, iter, sigma, lambda, voxelSize, origImage, numSamples)

p=[];


ffSupportInterp=ffSupportInterp(1:iter);

for i=1:length(image)
    
    p = vertcat(p,freeformParams{1}{i}{end}(:));
    
end

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



for j=1:10
    
    fprintf('Iteration...%d\n',j);
    
    sZ = size(freeformParams{1}{1}{end});
    
    numP = sZ(1)*sZ(2)*sZ(3)*sZ(4);
    
    options = foptions;     % Default options vector.
    options(1) = 1;		% Switch on diagnostics.
    options(5) = 0;		% Use persistence
    options(7) = 15;	% Number of steps in trajectory.
    options(14) = numSamples;	% Number of Monte Carlo samples returned.
    options(15) = 5;	% Number of samples omitted at start of chain.
    options(18) = 0.005;
    %options(9) = 1;
    

    [samples] = hmc(@logPosterior, p, options, @DlogPosterior, pts, image, templateImage, ffSpacing,ffSupportInterp, ffSupportCov, iter, sigma, lambda, voxelSize, sZ, freeformParams{1});
    
    templateImageNew = 0;
    jac = 0;
    
    ssd = 0;
    
    for k=1:numSamples
        for l=1:length(image)
            pSamp = samples(k,:);
            st = (numP*(l-1))+1;
            en = ((numP*(l-1))+1)+(numP-1);
            
            freeformParams{k}{l}{end} = reshape(pSamp(st:en),sZ);
        end
    end
    clear samples
    ffSupportCovNew=0;
    
    for k=1:numSamples
        
        fprintf('Max magnitude of deformation in sample...%d\t%f\n',k,max(freeformParams{k}{1}{end}(:)));
        
        for l=1:length(image)
            
            fprintf('Updating template for sample...%d\n',k);
    
            [temp, jacNew] = updateTempJacobian(freeformParams{k}{l}, ffSpacing, ffSupportInterp, iter, origImage{l}, voxelSize);
            
            
            mSave(:,:,:,l,k) = single(temp);
            
            ssd = ssd + sum((templateImage(:) - origImage{l}(:)).^2);
            
            templateImageNew = templateImageNew+temp; %clear temp
            jac = jac+jacNew; %clear jacNew

            
        end
        
        
        
    end
    
        
    pParam = [lambda];
    
    %update lambda
    pParam = minFunc(@logPosteriorLambda,pParam, optionsPred, pts, image, templateImage, ffSpacing,ffSupportInterp, ffSupportCov, iter, sigma, voxelSize, sZ, freeformParams, numSamples, ssd);
    
    
    %update I and sigma
    templateImage = templateImageNew./jac;
    sigma = sqrt(ssd/(numSamples*numel(templateImageNew)*length(image)));
    
    lambda = pParam;
    %ffSupportCov = ffSupportCovNew/(numSamples*length(image));
   
    
    tSave = single(templateImage);
    
    sum(isnan(tSave(:)))
    
     save(strcat('/home/pai/Zie_Template_Estimation/Out/meanImage_EM_KB_level',num2str(iter),'_iter',num2str(j),'.mat'),'tSave'); clear tSave
     save(strcat('/home/pai/Zie_Template_Estimation/Out/sampleImage_EM_KB_level',num2str(iter),'_iter',num2str(j),'.mat'),'mSave'); clear mSave
    
    %update the parameters
%     for m=1:length(image)
%         
%         temp=0;
%         
%         for n=1:numSamples
%         
%             temp = temp+freeformParams{n}{m}{end};
%             
%         end
%         
%         temp = temp./numSamples;
%         freeformParams{1}{m}{end} = temp;
%         
%         
%     end
        
    
        
        
    
    

    
end


                

