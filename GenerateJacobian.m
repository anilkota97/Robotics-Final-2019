function [J] = GenerateJacobian(Tmat,joints)
%UNTITLED4 Generate the Jacobian matrix for a general robotic arm
%   Tmat is a 4x4xn transformation matrix, with n being the number of
%   joints, and each layer of Tmat being the transformation matrix from 0
%   to n-1.
%   Joints is a char array containing the types of joints, either R or P

%Create the Jacobian Matrix
J = zeros(6,length(joints));

%Find the end effector pose
Pe = Tmat(1:3,4,end);


for i = 1:length(joints)
    %Check if the joint is prismatic
    if(strcmpi(joints(i),'p'))
        if(i == 1)
            Zcomp = [0;0;1];
        else
            Zcomp = Tmat(1:3,3,i-1);
        end
        %Fill in the Jacobian
        J(:,i) = [Zcomp;0;0;0];
    %Revolute Joints
    elseif(strcmpi(joints(i),'r'))
        
        %Find the P and the Z components
        if(i == 1)
            P = [0;0;0];
            Zcomp = [0;0;1];            
        else
            P = Tmat(1:3,4,i-1);
            Zcomp = Tmat(1:3,3,i-1);            
        end
        %Find the Jacobian component
        Jcomp = cross(Zcomp,(Pe - P));
        %Create the Jacobian entry
        J(:,i) = [Jcomp;Zcomp];
    end
end
        




end

