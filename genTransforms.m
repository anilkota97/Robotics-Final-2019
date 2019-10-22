function [T] = genTransforms(a,d,alpha,theta)
%UNTITLED6 Generate Transformation matrices from DH parameters
%   Detailed explanation goes here

T = zeros(4,4,length(a));
A = zeros(size(T));

for i = 1:length(a)
    %Create the series of A matrices for the transforms
    A(:,:,i) = [cosd(theta(i)), -sind(theta(i))*cosd(alpha(i)), sind(theta(i))*sind(alpha(i)), a(i)*cosd(theta(i));
        sind(theta(i)), cosd(theta(i))*cosd(alpha(i)), -cosd(theta(i))*sind(alpha(i)), a(1)*sind(theta(i));
        0 sind(alpha(i)), cosd(alpha(i)), d(i); 0 0 0 1];
end

%Create a starting matrix holding the first transformation matrix
prev = A(:,:,1);

for i = 1:length(a)
    
    if(i == 1)
        T(:,:,i) = prev;
    else 
        %Find the transformation matrix up to the current joint
        T(:,:,i) = prev*A(:,:,i);
        %Update the holder for the previous transformation matrix
        prev = T(:,:,i);

    end
end
                


end

