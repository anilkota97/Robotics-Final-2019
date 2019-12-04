clc; close all;
disp('Enter the following in matrix form:');
fprintf('\n\t(Note: Enter the initial configuration of the arm\n\t DO NOT ENTER CONSTANTS OR VARIABLES AS SYMBOLS!)\n');
a = input('Enter the "a" values of the DH parameter table: ');
d = input('Enter the "d" values of the DH parameter table: ');
alpha = input("Enter the '\alpha' values of the DH parameter table: ");
theta = input("Enter the '\theta' values of the DH parameter table: ");
ll = input('Enter the link lengths of the arm: ');
q_mins = input("Enter the minimum joint angle / distance matrix: ");
q_maxes = input("Enter the maximum joint angle / distance matrix: ");
workspace(a,d,alpha,theta,q_mins,q_maxes,ll,'draw');

function workspace(a,d,alpha,theta,qmins,qmaxes,linklens,ret_type)
if isempty(a)
    disp('The matrix containing the "a" values of DH parameters is empty !!');
elseif isempty(d)
    disp('The matrix containing the "d" values of DH parameters is empty !!');
elseif isempty(alpha)
    disp('The matrix containing the "alpha" values of DH parameters is empty !!');
elseif isempty(theta)
    disp('The matrix containing the "theta" values of DH parameters is empty !!');
elseif isempty(qmins)
    disp('The matrix containing minimum joint angles/distances is empty !!');
elseif isempty(qmaxes)
    disp('The matrix containing maximum joint angles/distances is empty !!');
elseif isempty(linklens)
    disp('The matrix containing link lengths is empty !!');
else
    sia = size(a); sid = size(d); sial = size(alpha); sit = size(theta);
    siaqmi = size(qmins); siqma = size(qmaxes); sill = size(linklens);
    if sia(1,1) == 1
        a = a';
    end
    if sid(1,1) == 1
        d = d';
    end
    if sial(1,1) == 1
        alpha = alpha';
    end
    if sit(1,1) == 1
        theta = theta';
    end
    if siaqmi(1,1) == 1
        qmins = qmins';
    end
    if siqma(1,1) == 1
        qmaxes = qmaxes';
    end
    if sill(1,1) == 1
        linklens = linklens';
    end
    dh_table = [a d alpha theta];
    ws_mat = get_workspace_coordinates(dh_table,qmins,qmaxes,linklens);
    X_Coordinates = ws_mat(:,1); Y_Coordinates = ws_mat(:,2); Z_Coordinates = ws_mat(:,3);
    if strcmpi(ret_type,'draw')
        grid on; hold on;
        [xc, yc, zc] = cylinder([0.2 0.5]);
        surface(xc, yc, zc,'EdgeColor','none','FaceColor','red');
        plot3(X_Coordinates,Y_Coordinates,Z_Coordinates,'k*','LineWidth',3);
        xlabel('X-AXIS');ylabel('Y-AXIS');zlabel('Z-AXIS');
        title('WORKSPACE OF THE MANIPULATOR / ARM');
        legend('Base of the arm','Workspace of the manipulator / arm');
    elseif strcmpi(ret_type,'display')
        fprintf('\nThe position co-ordinates of the end_effector in cartesian space is:\n\n');
        T = table(X_Coordinates,Y_Coordinates,Z_Coordinates);
        disp(T);
    else
        fprintf('\n\tThe entered return type is unknown!\nThe available types are:\n\t"draw"\n\t"display"\n');
    end
end
end

function ws = get_workspace_coordinates(dh_tab,qmin,qmax,linklen)
if isempty(dh_tab)
    disp('The matrix containing DH Parameters is empty !!');
elseif isempty(qmin)
    disp('The matrix containing minimum joint angles/distances is empty !!');
elseif isempty(qmax)
    disp('The matrix containing maximum joint angles/distances is empty !!');
elseif isempty(linklen)
    disp('The matrix containing link lengths is empty !!');
else
    resolution = 10; Jl = get_type_of_joints(dh_tab);
    qminsi = size(qmin); qmaxs = size(qmax);
    if qminsi(1,1) == 1
        qmin = qmin';
    elseif qmaxs(1,1) == 1
        qmax = qmax';
    end
    if length(qmin(:,1)) == length(qmax(:,1))
        [X,Y,Z] = get_coordinates_equations(dh_tab);
        ws = []; count = 1; q = qmin;
        for s = 1:1:length(qmin)
            while count <= length(min(qmin):1:max(qmax))
                for i = s:1:length(q)
                    if q(i,1)+resolution <= qmax(i,1)
                        q(i,1) = q(i,1)+resolution;
                    elseif q(i,1)+resolution > qmax(i,1)
                        q(i,1) = qmax(i,1);
                    else
                        disp('Joint angles/distances are not updating');
                    end
                end
                dh_tab = update_dh(dh_tab,q,Jl,linklen);
                [x,y,z] = get_coordinates(dh_tab,X,Y,Z);
                ws = [ws;x y z];
                count = count + 1;
            end
        end
    end
end
end

function [X,Y,Z] = get_coordinates(dh,x_eqn,y_eqn,z_eqn)
if isempty(dh)
    disp('The matrix containing DH Parameters is empty !!');
else
    a = dh(:,1); d = dh(:,2); alpha = dh(:,3); theta = dh(:,4);
    sym_a = sym('a',[length(a),1],'rational');
    sym_d = sym('d',[length(d),1],'rational');
    sym_alpha = sym('alpha',[length(alpha),1],'rational');
    sym_theta = sym('theta',[length(theta),1],'rational');
    variables_in_x = symvar(x_eqn);
    variables_in_y = symvar(y_eqn);
    variables_in_z = symvar(z_eqn);
    values_to_sub_in_x = zeros(length(theta),1);
    values_to_sub_in_y = zeros(length(theta),1);
    values_to_sub_in_z = zeros(length(theta),1);
    for i = 1:1:length(variables_in_x)
        for j = 1:1:length(theta)
            if variables_in_x(i) == sym_a(j)
                values_to_sub_in_x(i) = a(j);
            elseif variables_in_x(i) == sym_d(j)
                values_to_sub_in_x(i) = d(j);
            elseif variables_in_x(i) == sym_alpha(j)
                values_to_sub_in_x(i) = alpha(j);
            elseif variables_in_x(i) == sym_theta(j)
                values_to_sub_in_x(i) = theta(j);
            end
        end
    end
    for i = 1:1:length(variables_in_y)
        for j = 1:1:length(theta)
            if variables_in_y(i) == sym_a(j)
                values_to_sub_in_y(i) = a(j);
            elseif variables_in_y(i) == sym_d(j)
                values_to_sub_in_y(i) = d(j);
            elseif variables_in_y(i) == sym_alpha(j)
                values_to_sub_in_y(i) = alpha(j);
            elseif variables_in_y(i) == sym_theta(j)
                values_to_sub_in_y(i) = theta(j);
            end
        end
    end
    for i = 1:1:length(variables_in_z)
        for j = 1:1:length(theta)
            if variables_in_z(i) == sym_a(j)
                values_to_sub_in_z(i) = a(j);
            elseif variables_in_z(i) == sym_d(j)
                values_to_sub_in_z(i) = d(j);
            elseif variables_in_z(i) == sym_alpha(j)
                values_to_sub_in_z(i) = alpha(j);
            elseif variables_in_z(i) == sym_theta(j)
                values_to_sub_in_z(i) = theta(j);
            end
        end
    end
    X = x_eqn; Y = y_eqn; Z = z_eqn;
    for i = 1:1:length(values_to_sub_in_x)
        X = subs(X,variables_in_x(i),values_to_sub_in_x(i));
    end
    for i = 1:1:length(values_to_sub_in_y)
        Y = subs(Y,variables_in_y(i),values_to_sub_in_y(i));
    end
    for i = 1:1:length(values_to_sub_in_z)
        Z = subs(Z,variables_in_z(i),values_to_sub_in_z(i));
    end
    X = double(X); Y = double(Y); Z = double(Z);
end
end

function [X,Y,Z] = get_coordinates_equations(dh)
if isempty(dh)
    disp('The matrix containing DH Parameters is empty !!');
else
    a = dh(:,1); d = dh(:,1); alpha = dh(:,3); theta = dh(:,4);
    sym_a = sym('a',[length(a),1],'rational');
    sym_d = sym('d',[length(d),1],'rational');
    sym_alpha = sym('alpha',[length(alpha),1],'rational');
    sym_theta = sym('theta',[length(theta),1],'rational');
    % Calculating the Transformation matrix
    for i = 1:1:length(sym_theta)
        A{i} = [[cosd(sym_theta(i)) -(sind(sym_theta(i))*cosd(sym_alpha(i)))  (sind(sym_theta(i))*sind(sym_alpha(i))) (sym_a(i)*cosd(sym_theta(i)))];
                [sind(sym_theta(i))  (cosd(sym_theta(i))*cosd(sym_alpha(i))) -(cosd(sym_theta(i))*sind(sym_alpha(i))) (sym_a(i)*sind(sym_theta(i)))];
                [0 sind(sym_alpha(i)) cosd(sym_alpha(i)) sym_d(i)];
                [0 0 0 1]];
    end
    for i = 1:1:length(A)-1
        T = A{i}*A{i+1};
    end
    %Extracting the coordinates
    X = T(1,end);
    Y = T(2,end);
    Z = T(3,end);
end
end

function new_dh = update_dh(old_dh, q, Jl, linklens)
if isempty(old_dh)
    disp('The matrix containing DH Parameters is empty !!');
elseif isempty(q)
    fprintf('\nThe matrix containing updating q is empty !!\nError in minimum and maximum joint angles/distances !!\n');
elseif isempty(linklens)
    disp('The matrix containing link lengths is empty !!');
else
    a = old_dh(:,1); d = old_dh(:,2);
    alpha = old_dh(:,3); theta = old_dh(:,4);
    sa = size(a);sd = size(d);
    sal = size(alpha);st = size(theta);
    sll = size(linklens);
    if sa(1,1) == 1
        a = a';
    end
    if sd(1,1) == 1
        d = d';
    end
    if sal(1,1) == 1
        alpha = alpha';
    end
    if st(1,1) == 1
        theta = theta';
    end
    if sll(1,1) == 1
        linklens = linklens';
    end
    for i = 1:1:length(Jl)
        if strcmpi(Jl(i,1),'R')
            theta(i,1) = q(i,1);
        elseif strcmpi(Jl(i,1),'P')
            d(i,1) = (cos(alpha(i,1)))*(q(i,1) + linklens(i,1));
        end
    end
    new_dh = [a d alpha theta];
end
end

function Jllis = get_type_of_joints(dh)
if isempty(dh)
    disp('The matrix containing DH Parameters is empty !!');
else
    Jllis = []; a = dh(:,1); d = dh(:,2);
    thetas = dh(:,4); st = size(dh(:,4));
    if st(1,1) == 1
        thetas = thetas';
    end
    for i = 1:1:length(thetas)
        if thetas(i,1) ~= 0 && d(i,1) == 0
            Jllis = [Jllis;'R'];
        elseif thetas(i,1) == 0 && (a(i,1) ~= 0 || d(i,1) ~= 0)
            Jllis = [Jllis;'P'];
        end
    end
end
end
