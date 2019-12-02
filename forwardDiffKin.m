function [sing,symv, Ve, Jold] = forwardDiffKin(a, d, alpha, theta, jointTypes, jointVel)
%UNTITLED3 Compute robot forward kinematics
%   Detailed explanation goes here

%Generate the set of transformation matricies for the entire robot
[T, symTmat] = genTransforms(a,d,alpha,theta);

%Generate the Jacobian matrix from the transforms
[J, symJ] = GenerateJacobian(T,jointTypes, symTmat);

%Store the J matrix
Jold = J;

%Make the matrix square
temp = size(J);
numCol = temp(2);

colTracker = zeros(1,numCol);
rowtracker = zeros(1,6);

%Strip the extra rows so the jacobian can be performed
tempJ = rref(J);
for i = 1:6
    if(sum(tempJ(i,:)) == 0)
        rowtracker(i) = 1;
    end
end

counter = 1;
for i = 1:6
    if(rowtracker(counter) ~= 0)
        J(counter,:) = [];
        symJ(counter,:) = [];
        rowtracker(counter) = [];
    else
        counter = counter +1;
    end
end
    


% if(numCol ~= 6)
%     for i = 1:6
%         for j = 1:numCol
%             if(J(i,j) == 0)
%                 colTracker(j) = 1;
%             end
%         end
%         if(sum(colTracker) == numCol)
%             rowtracker(i) = 1;
%         end
%         
%         colTracker = zeros(1,numCol);
%     end
%     
%     temp = 1;
%     for i = 1:length(rowtracker)
%         if(rowtracker(temp) == 1)
%             J(temp,:) = [];
%             symJ(temp,:) = [];
%             rowtracker(temp) = [];
%             
%         else
%             temp = temp+1;
%         end
%         
%     end
% end



%Find the singularities
%first find the determinant of the Jacobian matrix
tempJ = det(J);
tempsym = det(symJ);

tempsym = simplify(tempsym);

eqn = tempsym == 0;

symv = symvar(tempsym);

sol = solve(eqn,symv);

for i = 1:length(symv)
    sing(:,:,i) = getfield(sol,char(symv(i)));
end


%Find the end effector velocity
%calculates the end effector linear and rotational velocity
Ve = Jold*jointVel;


end

