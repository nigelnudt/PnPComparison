clear
close all
%依据文献Efficient Linear Solution of Exterior Orientation，进行的摄像机位姿估计算法
%文献作者：Paul D. Fiore, Member, IEEE
%上一个版本的程序是直接调用SVD进行的零空间求解，但是对于含噪的数据而言，需要较多的点进行一个最小二乘的求解
%在具有较多点的情况下，最小二乘解的SVD求解比较复杂，因而需要采用其他方式进行分析，本程序按照作者的步骤完整进行仿真
addpath 'D:\计算机视觉测量\张志龙老师合作\EPnP_matlab\data';
close all
clear
%线性求解摄像机的位姿参数，参考文献
%首先介绍四点法，通过任选三个点构造一个一元四次方程，将这些方程的
%Linear N-Point Camera Pose Determination，Long Quan and Zhongdan Lan
addpath 'D:\计算机视觉测量\张志龙老师合作\EPnP_matlab\data';
% clc;
% clear all;
%此处将张志龙老师的仿真数据利用Long Quan的线性方法进行求解，看两者之间的误差
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cFileName = 'D:\计算机视觉测量\张志龙老师合作\N点位姿数据 2013-09-13\Pose_PointPattern2.dat';
%Pose_PointPattern1.dat
%Pose_PointPattern2.dat
%Pose_PointPattern3.dat
%Pose_PointPattern4.dat
%Pose_PointPattern5.dat
fid = fopen(cFileName,'r');
[Alfa, nCount] = fread(fid,1,'double');
[Beta, nCount] = fread(fid,1,'double');
[Gama, nCount] = fread(fid,1,'double');
[RR,   nCount] = fread(fid,[3,3],'double');
[tt,   nCount] = fread(fid, 3,'double');
[num,  nCount] = fread(fid, 1,'int');
[DA,   nCount] = fread(fid,[5,num],'double');
fclose(fid);
RR = RR';
DA = DA';
num_point=length(DA);
point_world_coordinate_h=[DA(:,3:5)';ones(1,num_point)];
for point_index=1:num_point%将世界坐标系转化为摄像机坐标系下的坐标以及图像坐标
    point_camera_coordinate(:,point_index)=RR*(point_world_coordinate_h(1:3,point_index)-[0,3,0]');
    point_camera_distance(point_index)=norm( point_camera_coordinate(:,point_index));
    point_image(:,point_index)=point_camera_coordinate([1,2],point_index)/point_camera_coordinate(3,point_index);
    point_image_h(1:2,point_index)=point_image(:,point_index);
    point_image_h(3,point_index)=1;
end
% for point_index=1:num_point%将世界坐标系转化为摄像机坐标系下的坐标以及图像坐标
%     point_camera_coordinate(:,point_index)=R*point_world_coordinate_h(:,point_index);
% %     distance(point_index)=norm(point_camera_coordinate(1:3,point_index));
%     point_image(:,point_index)=point_camera_coordinate(1:2,point_index)/point_camera_coordinate(3,point_index);
%     point_image_h(1:2,point_index)=point_image(:,point_index);
%     point_image_h(3,point_index)=1;
% end
[UU,DD,VV]=svd(point_world_coordinate_h);
[aa,dd]=qr(point_world_coordinate_h');
matrix2=eye(10)-aa(:,1:4)*aa(:,1:4)';
matrix_w=VV(:,4:end);%在此得到加权矩阵，以便将文中的深度值进行更好的求解
point_image=DA(:,1:2)';

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
point_camera_coordinate_ratio=point_camera_coordinate(3,:)/point_camera_coordinate(3,1);
b_matrix1=[l';l';l'].*point_image_h;
sum_ab=0;
sum_aa=0;
b_mean=mean(b_matrix1,2);
b_matrix=b_matrix1-repmat(b_mean,1,10);
a_mean=mean(DA(:,3:5)',2);
a_matrix=DA(:,3:5)'-repmat(a_mean,1,10);
for point_index=1:num_point
    sum_ab=sum_ab+norm(b_matrix(:,point_index))*norm(a_matrix(:,point_index));
    sum_aa=sum_aa+norm(norm(a_matrix(:,point_index)))*norm(a_matrix(:,point_index));
end
scale=sum_ab/sum_aa;
[UR,DR,VR]=svd(a_matrix*b_matrix');
% [UR1,DR1,VR1]=svd(DA(:,3:5)'*b_matrix1');
rotation_matrix1=VR*UR'*[1,0,0;0,-1,0;0,0,1];
rotation_matrix2=VR*[1,0,0;0,1,0;0,0,-1]*UR';
diff1=rotation_matrix1-RR;
diff2=rotation_matrix2-RR;
% rotation_matrix1=VR1*UR1';
% rotation_matrix=rotation_matrix*[1,0,0;0,-1,0;0,0,1];
T1=1/scale*rotation_matrix1'*b_mean-a_mean;
T2=1/scale*rotation_matrix2'*b_mean-a_mean;
pause(0.1)