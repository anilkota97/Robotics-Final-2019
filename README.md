# Robotics-Final-2019
A repository for the Final Project for MAE 547 Robotics

MATLAB Guide tutorial: https://www.youtube.com/watch?v=cl0AcnN3Bmk&t=1640s

#### a) Homogeneous:
Methods<br />

1) About vector:<br />
  i) Pure translation<br /> 
  ii) Pure Rotation<br />
  iii) Both<br />
  
2) About frame<br /> 
    i) Pure translation<br /> 
    ii) Pure Rotation<br />
    iii) Both<br />
  
#### b) Euler:<br />
User input:<br />

1) Forward kinematics:<br />
 ZYZ<br />
 ZYX<br />
 angles<br />
        OR<br />

2) Inverse kinematics:<br />
rotation matrix<br />

#### c) Forward Kinematics:<br />
User input:<br />

1) DH parameters<br />

  OR<br />

2) Robot Parameters<br /> 

i) Number joints <br />
ii) Types of joints<br />
iii) Link lengths<br />
iv) z axis: Steps:  User input ----> vector form input with respect to the reference frame that we provide to the user<br />
                    we calculate x axis <br />
                    then we can find alpha <br />
 
#### d) Workspace: <br />
User input:
Limits of joint configurations (q values)
Types of joint
Link Lengths

#### e) Inverse Kinematics: <br />
User input:
End effector pose
End effector orientation w.r.t base frame
(Hint: use ikine funtion)

#### f) Differential Kinematics: <br />
User input:<br />

1) DH parameters<br />

  OR<br />

2) Robot Parameters<br /> 

i) Number joints <br />
ii) Types of joints<br />
iii) Link lengths<br />
iv) z axis: Steps:  User input ----> vector form input with respect to the reference frame that we provide to the user<br />
                    we calculate x axis <br />
                    then we can find alpha <br />
3) Translational velocities or Angular velocities or the user should specify q interms of t<br /> 

#### g) Inverse differential kinematics and inverse kinematics using Jacobians: <br />
Homework4 problem: user should give (pd and phid)

#### i) Manipulator control: <br />
1) Operational space motion control (Inverse Dynamics Model)
2) Indirect force control via compliance control 
