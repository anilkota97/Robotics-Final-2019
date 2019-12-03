clc; close all;
inverse_diff_kine_using_jacobian(pd,phid,tmax,numJoints, jointTypes, linkLengths, zAxis, linkDir,linklens)
function inverse_diff_kine_using_jacobian(pd,phid,tmax,numJoints, jointTypes, linkLengths, zAxis, linkDir)
if isempty(pd)
    disp('pd matrix is empty !!');
elseif isempty(phid)
    disp('phid matrix is empty !!');
elseif isempty(tmax)
    disp('Maximum time is empty !!');
elseif isempty(a)
    disp('a matrix of DH parameters table is empty !!');
elseif isempty(d)
    disp('d matrix of DH parameters table is empty !!');
elseif isempty(alpha)
    disp('alpha matrix of DH parameters table is empty !!');
elseif isempty(theta)
    disp('theta matrix of DH parameters table is empty !!');
elseif isempty(linkLengths)
    disp('linklengths matrix is empty !!');
else
    size_of_a = size(a); size_of_d = size(d);
    size_of_alpha = size(alpha); size_of_theta = size(theta);
    if size_of_a(1,1) == 1
        a = a';
    end
    if size_of_d(1,1) == 1
        d = d';
    end
    if size_of_alpha(1,1) == 1
        alpha = alpha';
    end
    if size_of_theta(1,1) == 1
        theta = theta';
    end
    
    [dh] = calcDHfromRobot(numJoints, jointTypes, linkLengths, zAxis, linkDir);
    step = 2;
    iter = 0:step:tmax;
    syms t; tval = 0;
    
    % Extracting the inputs
    pd = pd(2:length(pd)-1); phid = phid(2:length(phid)-1);
    pd = strsplit(pd); phid = strsplit(phid);
    pd = str2sym(pd); phid = str2sym(phid);
    pd = pd'; phid = phid';
    % Extracting the pd_dot and phid_dot from pd and phid
    pd_dot = diff(pd,t); phid_dot = diff(phid,t);
    
    q = zeros((length(pd)+length(phid)),length(iter));
    K = eye((length(pd)+length(phid)));
    
    % Obtaining the initial configuration
    init_joint_config = [pd;phid];
    init_joint_config_ = subs(init_joint_config,t,0);
    q(:,1) = double(init_joint_config_);
    
    % Obtaining the initial Position and Rotation of end effector
    T = genTransforms(numJoints, jointTypes, linkLengths, zAxis, linkDir);
    pos = T(1:3,end,end);
    R = T(1:3,1:3,end);
    angle_phi = atan2(R(3,2),R(3,3));
    angle_theta = atan2(-R(3,1),real(sqrt((R(3,2)^2) + (R(3,3)^2))));
    angle_psi = atan2(R(2,1),R(1,1));
    rot = [angle_phi;angle_theta;angle_psi];
    Pv = [pos;rot];
    
    
    for i = 1:length(iter)
        %Getting the forward kinematics for the robot
        xe(:,i)= get_ee_fkine(Pv, q(:,1), jointTypes);
        
        % Getting the analytical jacobian
        Ja = get_analytical_jacobian(Pv, q(:,1), jointTypes);
        
        % substituting the tvalues in the following
        pd__(:,i) = subs(pd(:,1),tval);
        phid__(:,i) = subs(phid(:,1),tval);
        pd_dot__(:,i) = subs(pd_dot(:,1),tval);
        phid_dot__(:,i) = subs(phid_dot(:,1),tval);
        
        % Converting from sym to double
        pd_(:,i) = double(pd__(:,i));
        phid_(:,i) = double(phid__(:,i));
        pd_dot_(:,i) = double(pd_dot__(:,1));
        phid_dot_(:,i) = double(phid_dot__(:,i));
        
        % Logic for evaluating q
        e(:,i) = [pd_(:,i); phid_(:,i)] - xe(:,i);
        xd_dot = [pd_dot_(:,i); phid_dot_(:,i)];
        qdot = inv(Ja) * (xd_dot + K * e(:,i));
        q(:,i+1) = q(:,i) + qdot * step;
        tval = tval + step;
    end
    %     disp(q(:,end));
end
end

function Pvect = get_ee_fkine(pv, q, jointTypes)
Pvect = [];
for i = 1:1:length(pv)
    exp_to_sub = pv(i);
    variables = symvar(exp_to_sub);
    values = zeros(1,length(variables));
    theta = zeros(1,length(variables));
    d = zeros(1,length(variables));
    sym_theta = sym('theta', [1,length(theta)]);
    sym_d = sym('d', [1,length(d)]);
    for j = 1:1:length(jointTypes)
        if strcmpi(jointTypes(j),'R')
            theta(i) = q(i);
        elseif strcmpi(jointTypes(j),'P')
            d(i) = q(i);
        end
    end
    for i = 1:1:length(variables)
        if variables(i) == sym_theta(i)
            values(i) = theta(i);
        elseif variables(i) == sym_d(i)
            values(i) = d(i);
    end
    sub_eqn = subs(exp_to_sub,variables,values);
    sub_eqn_ = double(sub_eqn);
    Pvect = [Pvect;sub_eqn_];
end
end

function J = get_analytical_jacobian(Pv, q, jointTypes)
if isempty(dh)
    disp('The DH Matrix is empty !!');
else
    sym_theta = sym('theta',[length(dh(:,4)),1],'rational');
    for i = 1:1:length(dh(:,4))
        if strcmpi(type_of_joint(i),'R')
            differential = [zeros(length(Pv)-1,1);1];
        elseif strcmpi(type_of_joint(i),'P')
            differential = zeros(length(Pv),1);
        end
        var_list = symvar(Pv);
        theta_in_var_list = false;
        for j = 1:1:length(var_list)
            if sym_theta(i,1) == var_list(j)
                theta_in_var_list = true;
                break;
            end
        end
        if theta_in_var_list
            differential = diff(Pv,sym_theta(i,1));
        end
        if i == 1
            J = differential;
        else
            J = [J differential];
        end
    end
    J = substitute_values_in_jacobian(J, q, jointTypes);
end
end

function J = substitute_values_in_jacobian(symb_J, q, jointTypes)
if isempty(symb_J)
    disp('The symbolic Jacobian is empty !!');
else
    J = symb_J;
    for col = 1:1:symb_J(1,:)
        substituted_col = get_ee_fkine(symb_J(:,col), q, jointTypes);
        J(:,col) = substituted_col;
    end
%     a = dh(:,1); d = dh(:,2); alpha = dh(:,3); theta = dh(:,4);
%     sym_a = sym('a',[length(a),1],'rational');
%     sym_d = sym('d',[length(d),1],'rational');
%     sym_alpha = sym('alpha',[length(alpha),1],'rational');
%     sym_theta = sym('theta',[length(theta),1],'rational');
%     [m,n] = size(symb_J);
%     row_after_which_there_are_no_zeros = 2;
%     J = zeros(m,n); rows_to_remove = [];
%     for i = 1:1:m
%         logic_count = 0;
%         for j = 1:1:n
%             eqn_to_sub = symb_J(i,j);
%             variables_in_eqn = symvar(eqn_to_sub);
%             values_to_sub = zeros(1,length(variables_in_eqn));
%             for k1 = 1:1:length(variables_in_eqn)
%                 for k2 = 1:1:length(theta)
%                     if variables_in_eqn(k1) == sym_a(k2)
%                         values_to_sub(k1) = a(k2);
%                         break;
%                     elseif variables_in_eqn(k1) == sym_d(k2)
%                         values_to_sub(k1) = d(k2);
%                         break;
%                     elseif variables_in_eqn(k1) == sym_alpha(k2)
%                         values_to_sub(k1) = alpha(k2);
%                         break;
%                     elseif variables_in_eqn(k1) == sym_theta(k2)
%                         values_to_sub(k1) = theta(k2);
%                         break;
%                     end
%                 end
%             end
%             substituted_eqn = subs(eqn_to_sub,variables_in_eqn,values_to_sub);
%             substituted_eqn_ = double(substituted_eqn);
%             J(i,j) = substituted_eqn_;
%             if i > row_after_which_there_are_no_zeros
%                 if J(i,j) == 0
%                     logic_count = logic_count + 1;
%                 end
%             end
%         end
%         % logic for recording which row to delete
%         if i > row_after_which_there_are_no_zeros
%             if logic_count == length(J(i,:))
%                 rows_to_remove = [rows_to_remove i];
%             end
%         end
%     end
%     % logic for deleting a row in angles if all are zeros
%     for i = 1:1:length(rows_to_remove)
%         offset = 0;
%         J(rows_to_remove(i) - offset,:) = [];
%         for j = i+1:1:length(rows_to_remove)
%             rows_to_remove(j) = rows_to_remove(j) - 1;
%         end
%     end
end
end

function [symDH] = calcDHfromRobot(numJoints, jointTypes, linkLengths, zAxis, linkDir)
%UNTITLED12 Calculate DH parameters from Robot parameters
%   Takes input from the base frame to the end effector from. The base
%   frame is taken as co-incidental with the frame of the first joint.
%   Inputs:
%         : Zaxis - The Z Axis of the joints. For revolute joints it is the
%                   axis of rotation, For prismatic joints it is the axis
%                   of translations. Runs from first joint to the end
%                   effector (inclusive)
%         : jointTypes - The types of joints, From the first joint to the
%                        end effector (inclusive)
%         : linkLengths - The length of the links between the joints, Only
%                         number of links
%         : numJoints - The number of joints in the system, runs from the
%                       first joint to the end effector (inclusive)
%         : linkDir - The direction of the links between the joints. From
%                     first joint to the end effector (inclusive). Only the
%                     number of links
%
%   Outputs:
%          : Outputs a DH table with symbols in the place of the variables

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
        
        %-------Assign the z axis vector
        z_i = zAxis(:,:,i);
        
        %         fprintf('The Z axis is:\n')
        %         disp(z_i)
        
        %-------Find the common normal
        cn = cross(z_i,prevZ);
        
        %         fprintf('The common normal is:\n')
        %         disp(cn)
        
        %-------Find the X axis
        if(sum(cn) == 0)
            %Assign an arbitrary x_i for parallel z axis
            %x_i = [1;0;0];
            if(z_i(3) ~= 0)
                x_i = [1;0;0];
            elseif(z_i(2) ~= 0)
                x_i = [1;0;0];
            elseif(z_i(1) ~= 0)
                x_i = [0;0;1];
            end
        else
            %X is defined along the common normal, with direction from
            %joint i to i+1
            x_i = cn;
            %multiply by the sign of the of the link to go from joint i to
            %i+1
            x_i = x_i.*sign(linkDir(i));
            
        end
        
        %         fprintf('The X axis was assigned to be:\n')
        %         disp(x_i)
        
        %-------Find Y_i
        y_i = cross(z_i,x_i);
        
        
        % %----------------Finding O_i and O_i prime
        %         %O_i is located at the intersection of z_i with the common normal
        %         %of z_i and z_i-1
        %         %If neither of the common normals are parallel to the Z axis or z-1
        %         %axis
        %         if(jointTypes(i-1) == 'P')
        %             symDH = subs(symDH,symDH(i-1,1),0);
        %
        %         elseif(sum(cross(cn,z_i))~= 0 && sum(cross(cn,prevZ)) ~= 0)
        %             %If the axes are not parallel, o_i is located at the current
        %             %joint, and O_i prime is located at the previos joint. This
        %             %means that a is the distance between the joints.
        %
        %             symDH = subs(symDH,symDH(i-1,1),linkLengths(i));
        %         else
        %
        %             symDH = subs(symDH, symDH(i-1,1),linkLengths(i));
        %         end
        %
        %         fprintf('A was assigned to be: \n')
        %         disp(symDH(i-1,1))
        
        %-------Calculating a
        tempLinkDir = linkDir(:,:,i-1);
        if(jointTypes(i-1) == 'P')
            %If the joint is prismatic, the value of a should be zero
            %because the dh table should have the value in d instead
            symDH = subs(symDH,symDH(i-1,1),0);
        elseif(tempLinkDir(1) ~= 0)
            %Substitute in the value of a if the distance along the xi axis
            %is not equal to zero.
            symDH = subs(symDH,symDH(i-1,1),linkLengths(i-1));
        else
            symDH = subs(symDH,symDH(i-1,1),0);
        end
        
        %         fprintf('A was assigned to be: \n')
        %         disp(symDH(i-1,1))
        
        
        
        
        %-------Finding the values of d
        tempLinkDir2 = linkDir(:,:,i-1);
        if(jointTypes(i-1) == 'P')
            
            %since it is already a variable, Don't assign a value
            %continue
            %Check the sign
            temp = sign(tempLinkDir2(3));
            temp2 = temp*symDH(i-1,2);
            %subs(symDH,symDH(i-1,2),temp*symDH(i-1,2));
            symDH(i-1,2) = temp2;
            
        elseif(tempLinkDir2(3) ~= 0)
            %symDH = subs(symDH,symDH(i-1,2),linkLengths(i-1));
            %Since it is already a variable, don't assign and value and
            %simply continue
            %Check the direction first
            if(jointTypes(i-1) == 'R')
                temp = -1;% sign(tempLinkDir2(3));
                fprintf('%d\n',temp);
                temp2 = temp*linkLengths(i-1);
                fprintf('%d\n',temp2);
                disp(symDH)
                %subs(symDH,symDH(i-1,2),temp2);
                symDH(i-1,2) = temp2;
                disp(symDH)
            else
                temp = sign(tempLinkDir2(3));
                temp2 = temp*symDH(i-1,2);
                subs(symDH,symDH(i-1,2),temp2);
            end
            
        elseif(sum(cn) == 0)
            %If the common normal is zero, d should be zero
            symDH = subs(symDH,symDH(i-1,2),0);
        else
            
            symDH = subs(symDH,symDH(i-1,2),linkLengths(i));
        end
        
        %         fprintf('D was assigned to be: \n')
        %         disp(symDH(i-1,2))
        
        %----------------Finding Alpha
        %alpha is zero if the axes are parallel
        if(z_i == zAxis(i-1))
            symDH = subs(symDH,symDH(i-1,3),0);
            
        else
            %Calculate the cosine of the angle between two vectors
            ca = (dot(z_i,prevZ))/(sqrt(sum(z_i.^2))*sqrt(sum(prevZ.^2)));
            %Find the angle between two vectors
            tempangle = acosd(ca);
            symDH = subs(symDH,symDH(i-1,3),acosd(ca));
        end
        
        %         fprintf('Alpha was assigned to be: %d\n', symDH(i-1,3))
        
        
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
        
        %         fprintf('Theta was chosen to be: \n')
        %         disp(symDH(i-1,4));
        
        %Update the previous values with the new values
        prevZ = z_i;
        prevX = x_i;
        prevY = y_i;
        
        
    end
    
    
end
end

function T = genTransforms(numJoints, jointTypes, linkLengths, zAxis, linkDir)
%UNTITLED6 Generate Transformation matrices from DH parameters
%   Detailed explanation goes here

    [symDH] = calcDHfromRobot(numJoints, jointTypes, linkLengths, zAxis, linkDir);
    
    % A contains the transformation matrices between successive frames(confirm if alpha,theta are in radians in symDH)
    A = sym(zeros(4,4,numJoints-1));
    for i = 1:numJoints-1
        A(:,:,i) = [cosd(symDH(i,4)), -sind(symDH(i,4))*cosd(symDH(i,3)), sind(symDH(i,4))*sind(symDH(i,3)), symDH(i,1)*cosd(symDH(i,4));
                    sind(symDH(i,4)), cosd(symDH(i,4))*cosd(symDH(i,3)), -cosd(symDH(i,4))*sind(symDH(i,3)), symDH(1,1)*sind(symDH(i,4));
                    0, sind(symDH(i,3)), cosd(symDH(i,3)), symDH(i,2); 
                    0, 0, 0, 1];
    end
    
    % T contains the transformation matrices of frames w.r.t the base frame 
    T = sym(zeros(4,4,numJoints-1));
    prev = A(:,:,1);
    for i = 1:numJoints-1
        if i == 1
            T(:,:,i) = prev;
        else
            T(:,:,i) = prev * A(:,:,i);
            prev = T(:,:,i);
        end
    end
                


end

function [pos_e, T] = f_kinematics(numJoints, jointTypes, linkLengths, zAxis, linkDir, q)
[symDH] = calcDHfromRobot(numJoints, jointTypes, linkLengths, zAxis, linkDir);

% Substituting user input joint variables in DH parameters (confirm if numJoints includes end effector or not)
k = symvar(symDH);
[rows,cols] = size(symDH);
for i = 1:rows
    [row,col] = find(symDH == k(i));
    symDH(row,col) = q(i);
end

% A contains the transformation matrices between successive frames(confirm if alpha,theta are in radians in symDH)
A = sym(zeros(4,4,rows));
for i = 1:rows
    A(:,:,i) = [cosd(symDH(i,4)), -sind(symDH(i,4))*cosd(symDH(i,3)), sind(symDH(i,4))*sind(symDH(i,3)), symDH(i,1)*cosd(symDH(i,4));
        sind(symDH(i,4)), cosd(symDH(i,4))*cosd(symDH(i,3)), -cosd(symDH(i,4))*sind(symDH(i,3)), symDH(1,1)*sind(symDH(i,4));
        0, sind(symDH(i,3)), cosd(symDH(i,3)), symDH(i,2);
        0, 0, 0, 1];
end

% T contains the transformation matrices of frames w.r.t the base frame
T = sym(zeros(4,4,rows));
prev = A(:,:,1);
for i = 1:rows
    if i == 1
        T(:,:,i) = prev;
    else
        T(:,:,i) = prev * A(:,:,i);
        prev = T(:,:,i);
    end
end

% Finding the end-effector position w.r.t base frame
pos_e = T(1:3,end,end);

end