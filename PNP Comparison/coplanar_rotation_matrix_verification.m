%测试同一平面上的点旋转后的结果是否与预设的旋转矩阵一致
clear
close all
addpath 'D:\计算机视觉测量\张志龙老师合作\EPnP_matlab\data';
alpha=50/180*pi;
beta=20/180*pi;
gamma=10/180*pi;
scale=1/450;
tx=-20;
ty=-20;
tz=500;
num_point=5;
R=return_Rt_matrix(alpha,beta,gamma,tx,ty,tz);
xworld=20*rand(1,num_point)-10;
yworld=20*rand(1,num_point)-10;
zworld=-(2*xworld+7*xworld)/4

% zworld=20*rand(1,num_point)-10;
point_world_coordinate=[xworld;yworld;zworld];
point_world_coordinate_h=[point_world_coordinate;ones(1,num_point)];
point_image=zeros(2,num_point);
point_image_h=zeros(3,num_point);
point_camera_coordinate=zeros(4,num_point);
for point_index=1:num_point%将世界坐标系转化为摄像机坐标系下的坐标以及图像坐标
    point_camera_coordinate(:,point_index)=R*point_world_coordinate_h(:,point_index);
end
b_matrix1=point_camera_coordinate(1:3,:);
sum_ab=0;
sum_aa=0;
b_mean=mean(b_matrix1,2);
b_matrix=b_matrix1-repmat(b_mean,1,num_point)
a_mean=mean(point_world_coordinate,2);
a_matrix=point_world_coordinate-repmat(a_mean,1,num_point)
for point_index=1:num_point
    sum_ab=sum_ab+norm(b_matrix(:,point_index))*norm(a_matrix(:,point_index));
    sum_aa=sum_aa+norm(norm(a_matrix(:,point_index)))*norm(a_matrix(:,point_index));
end
scale=sum_ab/sum_aa;
[UR,DR,VR]=svd(a_matrix*b_matrix');
% [UR1,DR1,VR1]=svd(DA(:,3:5)'*b_matrix1');
 AA=det(UR)*det(VR)
if AA<0
    rotation_matrix=VR*[1,0,0;0,1,0;0,0,-1]*UR';
else
    rotation_matrix=VR*UR';
end
diff=rotation_matrix-R(1:3,1:3);
% rotation_matrix1=VR1*UR1';
% rotation_matrix=rotation_matrix*[1,0,0;0,-1,0;0,0,1];
T=1/scale*rotation_matrix'*b_mean-a_mean;

%%
%%新的旋转由小杨来写，已经证明可行
num_point = length( a_matrix ); 
B_yang =  zeros( 3*num_point, 1 );
A_yang = zeros( 3*num_point, 9 );
a_matrix=a_matrix'
b_matrix=b_matrix'
for i=1:num_point,
    A_yang(3*i-2,:) =[a_matrix(i,1),a_matrix(i,2),a_matrix(i,3),0,0,0,0,0,0];
    A_yang(3*i-1,:) =[0,0,0,a_matrix(i,1),a_matrix(i,2),a_matrix(i,3),0,0,0];
    A_yang(3*i,:) =[0,0,0,0,0,0,a_matrix(i,1),a_matrix(i,2),a_matrix(i,3)];
    B_yang(3*i-2:3*i,1)=b_matrix(i,:)';
end
x_yang = inv(A_yang'*A_yang)*A_yang'*B_yang;
y_yang = sqrt(3) * x_yang / norm( x_yang );
Rotate_yang = [ y_yang(1) y_yang(2) y_yang(3) 
           y_yang(4) y_yang(5) y_yang(6);
           y_yang(7) y_yang(8) y_yang(9)];
diff_yang=Rotate_yang-R(1:3,1:3);
pause(0.1)