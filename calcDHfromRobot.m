function [a,d,alpha,theta] = calcDHfromRobot(numJoints, jointTypes, linkLengths, zAxis)
%UNTITLED12 Calculate DH parameters from Robot parameters
%   Detailed explanation goes here

%Display a reference frame and take input of the z axis about the reference
%frame



%figure(1)
%Calculate the parameters for each link
for i = 1:numJoints
    
    if( i == 1)
        %Assign the initial X, Y, Z based on the first value
        Z_i = zAxis(:,1,1);
        
        %Take the previous Z vector to be the unit vector in the Z
        %direction
        prevZ = [0;0;1];

        prevX = [1;0;0];
        prevY = [0;1;0];
    else
        
        z_i = zAxis(:,:,i);
        %Find the common normal
        cn = cross(z_i,prevZ);
        
        %Check if the z axis are parallel
        if(sum(cn) == 0)
            %Assign an arbitrary x_i
            x_i = [1;0;0];
        else
            %X is defined along the common normal, with direction from
            %joint i to i+1
            x_i = cn;
        end
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
            a(i-1) = linkLengths(i);
        else
            a(i-1) = linkLengths(i);
        end
        
        %Need to check this logic for finding the values of d, I am not
        %sure if it is entirely correct
        if(sum(cn) == 0)
            d(i-1) = 0;
        else
            d(i-1) = linkLengths(i);
        end
            
                
%----------------Finding Alpha        
        %Calculate the cosine of the angle between two vectors
        ca = (dot(z_i,prevZ))/(sqrt(sum(z_i.^2))*sqrt(sum(prevZ.^2)));
        %Find the angle between two vectors
        alpha(i-1) = acosd(ca);
        

%----------------Finding Theta        
        %Theta will be zero for all prismatic joints, and a variable for
        %all revolute joints
        if(jointTypes(i) == 'P')
            %Set to zero because prismatic
            theta(i-1) = 0;
        else
            %Revolute joints
            ct = (dot(x_i,prevx))/(sqrt(sum(x_i.^2))*sqrt(sum(prevx.^2)));
            theta(i-1) = acosd(ct);            
        end
        
        
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

