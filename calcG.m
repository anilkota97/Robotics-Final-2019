function [G] = calcG(gdir, linkMass, motorMass, JpLi, JpMi, numJoints)
%UNTITLED6 Calculates the Gravity component for the manipulator dynamics 
%   Input:
%        :  gdir - A unit vector for the direction of gravity. Include the
%                 negative sign if it is needed.
%        :  linkMass - A Vector containing the mass of all the links of the
%                      system
%        :  motorMass - A vector containing the mass of all the motors of
%                       the system
%        :  JpLi - A Three Dimensional Matrix containing the Jp matrix of
%                  the links, and the third dimension contains the further
%                  ith matrices
%        :  JpMi - A Three Dimensional matrix containing the Jp matrix of
%                  the motors, and the third dimension contains the futher
%                  ith matricies
%        :  NumJoints - The number of joints in the system. Used to define
%                       the legth of the G output vector and the number of
%                       iterations
%   Output:
%         : G - A 1 x numJoints matrix containing the symbolic G values 

%Create a symbolic gravity acceleration
syms g
%Create the symbolic gravity vector
g0 = g*gdir;

%Initialize G
G = 0;
%Loop over the number of joints
for i = 1:numJoints
    %Perform the summation over G
    G = G + -linkMass(i)*transpose(g0)*JpLi(:,:,i) - motorMass(i)*transpose(g0)*JpMi(:,:,i);
    
end

end

