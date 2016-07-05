function points = generateEvaluationPoints(voxelSize, dim, density)
  % function points = generateEvaluationPoints(voxelSize, dim, density)
  %
  % Creates a number of uniformly spaced evaluation points.
  %
  % INPUT:
  % - voxelSize: The size of a voxel in mm (1x3 vector).
  % - dim      : The size of the scan in mm (1x3 vector).
  % - density  : The sampling rate in mm (1x3 vector).
  %
  % OUTPUT:
  % - points    : Matrix of evaluation points.
  
  
  [x1, x2, x3] = ndgrid(...
    voxelSize(1) : density(1) : dim(1),...
    voxelSize(2) : density(2) : dim(2),...
    voxelSize(3) : density(3) : dim(3));

  points = [x1(:) x2(:) x3(:)];
end
