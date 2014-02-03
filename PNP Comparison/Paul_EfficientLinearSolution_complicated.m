clear
close all
%依据文献Efficient Linear Solution of Exterior Orientation，进行的摄像机位姿估计算法
%文献作者：Paul D. Fiore, Member, IEEE
%上一个版本的程序是直接调用SVD进行的零空间求解，但是对于含噪的数据而言，需要较多的点进行一个最小二乘的求解
%在具有较多点的情况下，最小二乘解的SVD求解比较复杂，因而需要采用其他方式进行分析，本程序按照作者的步骤完整进行仿真
addpath 'D:\计算机视觉测量\张志龙老师合作\EPnP_matlab\data';
alpha=5/180*pi;
beta=5/180*pi;
gamma=1/180*pi;
tx=-20;
ty=-20;
tz=500;
num_point=6;
R=return_Rt_matrix(alpha,beta,gamma,tx,ty,tz);
xworld=400*rand(1,num_point)-200;
yworld=400*rand(1,num_point)-200;
zworld=400*rand(1,num_point)-200;
point_world_coordinate=[xworld;yworld;zworld];
point_world_coordinate_h=[point_world_coordinate;ones(1,num_point)];
point_camera_coordinate=zeros(4,num_point);
point_image=zeros(2,num_point);
point_image_h=zeros(3,num_point);
for point_index=1:num_point%将世界坐标系转化为摄像机坐标系下的坐标以及图像坐标
    point_camera_coordinate(:,point_index)=R*point_world_coordinate_h(:,point_index);
%     distance(point_index)=norm(point_camera_coordinate(1:3,point_index));
    point_image(:,point_index)=point_camera_coordinate(1:2,point_index)/point_camera_coordinate(3,point_index);
    point_image_h(1:2,point_index)=point_image(:,point_index);
    point_image_h(3,point_index)=1;
end
[UU,DD,VV]=svd(point_world_coordinate_h);
[aa,dd]=qr(point_world_coordinate_h')
matrix2=eye(6)-aa(:,1:4)*aa(:,1:4)'
matrix_w=VV(:,5:end);%在此得到加权矩阵，以便将文中的深度值进行更好的求解
for point_index=1:num_point
    matrix_for_l(1:2,point_index)=matrix_w(point_index,1)*point_image(:,point_index);
    matrix_for_l(3:4,point_index)=matrix_w(point_index,2)*point_image(:,point_index);
end
matrix_cc=point_world_coordinate_h*((matrix_w*matrix_w').*(point_image'*point_image))*point_world_coordinate_h'
[Ucc,Dcc,Vcc]=svd(matrix_cc)
l_complicate=point_world_coordinate_h'*Vcc(:,4);
l_complicate=l_complicate/l_complicate(1)
pause(0.1)
[Ul,Dl,Vl]=svd(matrix_for_l);
matrix_C=matrix_for_l*point_world_coordinate_h';
[UC,DC,VC]=svd(matrix_C);
l=point_world_coordinate_h'*VC(:,4)
l=l/l(1);
% distance=distance/distance(1);
point_camera_coordinate_ratio=point_camera_coordinate(3,:)/point_camera_coordinate(3,1)
pause(0.1)