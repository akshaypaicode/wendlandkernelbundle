function [im,voxSize] = loadImage(filename,ext)

switch ext
    
    case '.hdr'
        
        [im,targetInfo] = ReadData3D(filename);im=double(im);
        voxSize = targetInfo.PixelDimensions;
        
    case '.nii'
        im=MRIread(filename);
        voxSize=im.volres; im=im.vol;
        
    case '.mgz'
        
        im=MRIread(filename);
        voxSize=im.volres; im=im.vol;
      
    case '.nii.gz'
        im=MRIread(filename);  
        voxSize=im.volres; im=im.vol;


end

end