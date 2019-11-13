function [symDH] = calcDHfromRobot(numJoints, jointTypes, linkLengths, zAxis)
%UNTITLED12 Calculate DH parameters from Robot parameters
%   Goes from the base joint to the end effector

%Maybe we need to ask for the link direction as well? This will help
%fetermine o_i and o_i prime so that a and d can be determined more easily.
%If we know the link Direction, o_i is the current zAxis coordinate, and
%o_i prime is the curent z coordinate - linkLength*linkDirection. This can
%also help with directly calculating the layout of the robot, which will
%help find o_i and _i prime




%Assign the variables to symbolic parameters
symThetas = sym('Theta',[1, numJoints-1]);
symThetas = transpose(symThetas);
symD = sym('D',[1,numJoints-1]);
symD = transpose(symD);
%Create a symbollic DH matrix
symA = sym('a',[1 numJoints-1]);
symA = transpose(symA);
symAlpha = sym('alpha', [1 numJoints-1]);
symAlpha = transpose(symAlpha);


symDH = [symA symD symAlpha symThetas];

 

        
% d = zeros(numJoints-1,1);
% a = zeros(numJoints-1,1);
% alpha = zeros(numJoints-1, 1);
% theta = zeros(numJoints-1, 1);

%Calculate the parameters for each link
for i = 1:numJoints
    
    %For the first joint, there is no previous values yet, so it needs to
    %be handled separetly.
    if(i == 1)
        %Assign the initial X, Y, Z based on the first value
        Z_i = zAxis(:,1,1);
        
        %Take the previous Z vector to be the unit vector in the Z
        %direction
        prevZ = [0;0;1];
        prevX = [1;0;0];
        prevY = [0;1;0];
    else
        
        %Assign the z axis vector
        z_i = zAxis(:,:,i);
        
        fprintf('The Z axis is:\n')
        disp(z_i)
        
        %Find the common normal
        cn = cross(z_i,prevZ);
        
        fprintf('The common normal is:\n')
        disp(cn)
        
        %Check if the z axis are parallel
        if(sum(cn) == 0)
            %Assign an arbitrary x_i
            x_i = [1;0;0];
        else
            %X is defined along the common normal, with direction from
            %joint i to i+1
            x_i = cn;
        end
        
        fprintf('The X axis was assigned to be:\n')
        disp(x_i)
        
        %Find Y_i
        y_i = cross(z_i,x_i);
        
        
%----------------Finding O_i and O_i prime
        %O_i is located at the intersection of z_i with the common normal
        %of z_i and z_i-1
        %If neither of the common normals are parallel to the Z axis or z-1
        %axis
        if(sum(cross(cn,z_i))~= 0 && sum(cross(cn,prevZ)) ~= 0)
            %If the axes are not parallel, o_i is located at the current
            %joint, and O_i prime is located at the previos joint. This
            %means that a is the distance between the joints.
            
            symDH = subs(symDH,symDH(i-1,1),linkLengths(i));
        else
            
            symDH = subs(symDH, symDH(i-1,1),linkLengths(i));
        end
        
        fprintf('A was assigned to be: \n')
        disp(symDH(i-1,1))
        
        
        %Need to check this logic for finding the values of d, I am not
        %sure if it is entirely correct
        if(jointTypes(i-1) == 'P')
            
            %since it is already a variable, Don't assign a value
            %continue
            
            
        elseif(sum(cn) == 0)
            
            symDH = subs(symDH,symDH(i-1,2),0);
        else
            
            symDH = subs(symDH,symDH(i-1,2),linkLengths(i));
        end
            
        fprintf('D was assigned to be: \n')
        disp(symDH(i-1,2))
                
%----------------Finding Alpha
        %alpha is zero if the axes are parallel
        if(z_i == zAxis(i-1))
            symDH = subs(symDH,symDH(i-1,3),0);
            
        else
            %Calculate the cosine of the angle between two vectors
            ca = (dot(z_i,prevZ))/(sqrt(sum(z_i.^2))*sqrt(sum(prevZ.^2)));
            %Find the angle between two vectors
            
            symDH = subs(symDH,symDH(i-1,3),acosd(ca));
        end
        
        fprintf('Alpha was assigned to be: %d\n', symDH(i-1,3))
        

%----------------Finding Theta        
        %Theta will be zero for all prismatic joints, and a variable for
        %all revolute joints
        if(jointTypes(i-1) == 'P')
            %Set to zero because prismatic
            
            symDH = subs(symDH,symDH(i-1,4),0);
        else
            %Revolute joints, therefore it is a variable
            
            %Since the DH table is already a symbolic variable, leave it
            %continue
        end
        
        fprintf('Theta was chosen to be: \n')
        disp(symDH(i-1,4));
        
        %Update the previous values with the new values
        prevZ = z_i;
        prevX = x_i;
        prevY = y_i;
        
        
    end
%     if(i == 1)
%         %Initialize the origin
%         %prevx = [1;0;0];
%         %prevy = [0;1;0];
%         %prevz = [0;0;1];
%         %prevo = [0;0;0];
%         %Plot the first frame
%         %Plot the origin
% %         origin = [0,0,0];
% %         plot3(origin(1),origin(2),origin(3),'r*')
% %         grid on
% %         hold on
% %         
% %         %This behavior is incorrect, I want to move along the defined
% %         %z-axis
% %         newOrigin = origin;
% %         newOrigin(3) = newOrigin(3)+linkLengths(i);
% %         plot3(newOrigin(1), newOrigin(2), newOrigin(3),'b*')
% %         plot3([origin(1) prevx(1)], [origin(2) prevx(2)], [origin(3) prevx(3)],'b')
% %         plot3([origin(1) prevy(1)], [origin(2) prevy(2)], [origin(3) prevy(3)],'b')
% %         plot3([origin(1) prevz(1)], [origin(2) prevz(2)], [origin(3) prevz(3)],'b')
%         
%     else
%         %Find x_i
%         z_i = zAxis(:,:,i);
%         %have to take care of case when zi and zi-1 are parallel
%         x_i = cross(z_i,prevz);
%         %Take care of if the z Axis vectors are parallel
%         if(sum(x_i) == 0)
%             x_i = [1;0;0];
%         end
%         y_i = cross(z_i,x_i);
%         
%         %Migrate the frames
% %         newOrigin(3) = prevz(3) + linkLengths(i);
% %         plot3(newOrigin(1),newOrigin(2),newOrigin(3),'g*')
% %         plot3([newOrigin(1) x_i(1)], [newOrigin(2) x_i(2)], [newOrigin(3) x_i(3)],'g')
% %         plot3([newOrigin(1) y_i(1)], [newOrigin(2) y_i(2)], [newOrigin(3) y_i(3)],'g')
% %         plot3([newOrigin(1) z_i(1)], [newOrigin(2) z_i(2)], [newOrigin(3) z_i(3)],'g')        
%         
%         %Finding theta
%         if(jointTypes(i) == 'P')
%             theta(i-1) = 0;
%         else
%             ct = (dot(x_i,prevx))/(sqrt(sum(x_i.^2))*sqrt(sum(prevx.^2)));
%             theta(i-1) = acosd(ct);
%         end
%         
%         %Finding alpha
%         ca = (dot(z_i,prevz))/(sqrt(sum(z_i.^2))*sqrt(sum(prevz.^2)));
%         alpha(i-1) = acosd(ca);
%         
%         
%         %Find o_i: Intersection of zi, and the common normals of zi and
%         %zi-1
%         cn = cross(z_i,prevz);
%         
%         %Z axis are parallel
%         if(sum(cn) == 0)
%             o_i_prime = prevo;
%             o_i = prevo + z_i(i);
%             a(i-1) = linkLengths(i);
%         else
%             
%         end
%            
%         
%         prevx = x_i;
%         prevz = z_i(i);
%         prevy = y_i;
%     end

end

