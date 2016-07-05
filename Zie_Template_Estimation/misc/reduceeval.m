function pts=reduceeval(pts,seg,siz,voxelSize)

%if ~exist(siz,'var')
%    siz=6;
%end

pts = [pts(:,1)./voxelSize(1) pts(:,2)./voxelSize(2) pts(:,3)./voxelSize(3)];

seg=seg~=0;
se = createStrel(siz, siz);
dilated = imdilate(seg, se);
[x1,y1,z1] = ind2sub(size(seg), find(dilated==1));
pts = intersect(pts,[x1(:) y1(:) z1(:)],'rows');

pts = [pts(:,1).*voxelSize(1) pts(:,2).*voxelSize(2) pts(:,3).*voxelSize(3)];
% [xTmp,yTmp,zTmp] = ind2sub(size(seg), find(calcBandVol(seg, dilationLvl, erosionLvl)==1));
% ptsBand = [xTmp yTmp zTmp];
% clear dilated xTmp yTmp zTmp;
% pts = union(ptsBand, pts, 'rows');
end




function se = createStrel(radius, halfWidth)
  % function se = createStrel(radius, halfWidth)
  %
  % Generate binary 3d ball structuring element with radius=radius embedded and centred in a halfWidth*2+1 cube
  % 
  % INPUT:
  % - radius    : radius of the ball
  % - halfWidth : the embedding cube is has side size halfWidth*2+1
  %
  % OUTPUT: 
  % - se : the result structuring element

  [x,y,z] = meshgrid(-halfWidth:halfWidth,-halfWidth:halfWidth,-halfWidth:halfWidth); 
  idx = ( (x/radius).^2 + (y/radius).^2 + (z/radius).^2 ) <= 1;
  se = zeros(size(x)); 
  se(idx) = 1;
end

function bandVol = calcBandVol(binary, dilationLvl, erosionLvl)
% function calcBandPoints(binary, erosion, dilation)
%
% Subtract binary eroded volume from ditto dilated to volume with 1s in all band points and 0s otherwise
%
% INPUT:
% - binary      : binary input volume
% - dilationLvl : amount of dilation to define band
% - erosionLvl  : amounts of erosion the
se = createStrel(dilationLvl, dilationLvl);
dilated = imdilate(binary, se);
se = createStrel(erosionLvl, erosionLvl);
eroded = imerode(double(dilated), se);
bandVol = dilated-eroded;
end
