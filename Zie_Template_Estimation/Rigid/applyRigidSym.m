function [fwpoints,d,id,R]=applyRigidSym(p,points,center,inverse)

if inverse
    points(:,1) = points(:,1) - center(1) - p(10);
    points(:,2) = points(:,2) - center(2) - p(11);
    points(:,3) = points(:,3) - center(3) - p(12);
else
    points(:,1) = points(:,1) - center(1);
    points(:,2) = points(:,2) - center(2);
    points(:,3) = points(:,3) - center(3);
    
end


R = [p(1)+1 p(2) p(3); p(4) p(5)+1 p(6); p(7) p(8) p(9)+1];


if inverse
    R = inv(R);
else
    R = R;
end

fwpoints = points*R;


d{1} = [1 0 0; 0 0 0; 0 0 0];
d{2} = [0 1 0; 0 0 0; 0 0 0];
d{3} = [0 0 1; 0 0 0; 0 0 0];
d{4} = [0 0 0; 1 0 0; 0 0 0];
d{5} = [0 0 0; 0 1 0; 0 0 0];
d{6} = [0 0 0; 0 0 1; 0 0 0];
d{7} = [0 0 0; 0 0 0; 1 0 0];
d{8} = [0 0 0; 0 0 0; 0 1 0];
d{9} = [0 0 0; 0 0 0; 0 0 1];

if inverse
    id{1} = -R*[1 0 0; 0 0 0; 0 0 0]*R;
    id{2} = -R*[0 1 0; 0 0 0; 0 0 0]*R;
    id{3} = -R*[0 0 1; 0 0 0; 0 0 0]*R;
    id{4} = -R*[0 0 0; 1 0 0; 0 0 0]*R;
    id{5} = -R*[0 0 0; 0 1 0; 0 0 0]*R;
    id{6} = -R*[0 0 0; 0 0 1; 0 0 0]*R;
    id{7} = -R*[0 0 0; 0 0 0; 1 0 0]*R;
    id{8} = -R*[0 0 0; 0 0 0; 0 1 0]*R;
    id{9} = -R*[0 0 0; 0 0 0; 0 0 1]*R;
    
end



if inverse
    fwpoints(:,1) = fwpoints(:,1) + center(1);
    fwpoints(:,2) = fwpoints(:,2) + center(2);
    fwpoints(:,3) = fwpoints(:,3) + center(3);
else
    fwpoints(:,1) = fwpoints(:,1) + center(1) + p(10);
    fwpoints(:,2) = fwpoints(:,2) + center(2) + p(11);
    fwpoints(:,3) = fwpoints(:,3) + center(3) + p(12);
    
end

    
end

