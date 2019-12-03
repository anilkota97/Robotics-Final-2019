function [T,symT] = genTransforms2(a,d,alpha,theta)
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

%Create a symbollic matrix
symTheta = sym('Theta',[1, length(a)]);
symTheta = transpose(symTheta);
symD = sym('D',[1,length(a)]);
symD = transpose(symD);
%Create a symbollic DH matrix
symA = sym('a',[1 length(a)]);
symA = transpose(symA);
symAlpha = sym('alpha', [1 length(a)]);
symAlpha = transpose(symAlpha);

for i = 1:length(a)
    symAmat(:,:,i) = [cosd(symTheta(i)), -sind(symTheta(i))*cosd(symAlpha(i)), sind(symTheta(i))*sind(symAlpha(i)), symA(i)*cosd(symTheta(i));
        sind(symTheta(i)), cosd(symTheta(i))*cosd(symAlpha(i)), -cosd(symTheta(i))*sind(symAlpha(i)), symA(1)*sind(symTheta(i));
        0 sind(symAlpha(i)), cosd(symAlpha(i)), symD(i); 0 0 0 1];
end




%Create a starting matrix holding the first transformation matrix
prev = A(:,:,1);
prevsym = symAmat(:,:,1);

for i = 1:length(a)
    
    if(i == 1)
        T(:,:,i) = prev;
        symT(:,:,i) = prevsym;
    else 
        %Find the transformation matrix up to the current joint
        T(:,:,i) = prev*A(:,:,i);
        symT(:,:,i) = prevsym*symAmat(:,:,i);
        %Update the holder for the previous transformation matrix
        prev = T(:,:,i);
        prevsym = symT(:,:,i);

    end
end
                


end

