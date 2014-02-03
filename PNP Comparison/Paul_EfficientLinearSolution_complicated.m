clear
close all
%��������Efficient Linear Solution of Exterior Orientation�����е������λ�˹����㷨
%�������ߣ�Paul D. Fiore, Member, IEEE
%��һ���汾�ĳ�����ֱ�ӵ���SVD���е���ռ���⣬���Ƕ��ں�������ݶ��ԣ���Ҫ�϶�ĵ����һ����С���˵����
%�ھ��н϶�������£���С���˽��SVD���Ƚϸ��ӣ������Ҫ����������ʽ���з����������������ߵĲ����������з���
addpath 'D:\������Ӿ�����\��־����ʦ����\EPnP_matlab\data';
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
for point_index=1:num_point%����������ϵת��Ϊ���������ϵ�µ������Լ�ͼ������
    point_camera_coordinate(:,point_index)=R*point_world_coordinate_h(:,point_index);
%     distance(point_index)=norm(point_camera_coordinate(1:3,point_index));
    point_image(:,point_index)=point_camera_coordinate(1:2,point_index)/point_camera_coordinate(3,point_index);
    point_image_h(1:2,point_index)=point_image(:,point_index);
    point_image_h(3,point_index)=1;
end
[UU,DD,VV]=svd(point_world_coordinate_h);
[aa,dd]=qr(point_world_coordinate_h')
matrix2=eye(6)-aa(:,1:4)*aa(:,1:4)'
matrix_w=VV(:,5:end);%�ڴ˵õ���Ȩ�����Ա㽫���е����ֵ���и��õ����
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