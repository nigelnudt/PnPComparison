clear
close all
%��������Efficient Linear Solution of Exterior Orientation�����е������λ�˹����㷨
%�������ߣ�Paul D. Fiore, Member, IEEE
%��һ���汾�ĳ�����ֱ�ӵ���SVD���е���ռ���⣬���Ƕ��ں�������ݶ��ԣ���Ҫ�϶�ĵ����һ����С���˵����
%�ھ��н϶�������£���С���˽��SVD���Ƚϸ��ӣ������Ҫ����������ʽ���з����������������ߵĲ����������з���
addpath 'D:\������Ӿ�����\��־����ʦ����\EPnP_matlab\data';
close all
clear
%��������������λ�˲������ο�����
%���Ƚ����ĵ㷨��ͨ����ѡ�����㹹��һ��һԪ�Ĵη��̣�����Щ���̵�
%Linear N-Point Camera Pose Determination��Long Quan and Zhongdan Lan
addpath 'D:\������Ӿ�����\��־����ʦ����\EPnP_matlab\data';
% clc;
% clear all;
%�˴�����־����ʦ�ķ�����������Long Quan�����Է���������⣬������֮������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cFileName = 'D:\������Ӿ�����\��־����ʦ����\N��λ������ 2013-09-13\Pose_PointPattern2.dat';
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
for point_index=1:num_point%����������ϵת��Ϊ���������ϵ�µ������Լ�ͼ������
    point_camera_coordinate(:,point_index)=RR*(point_world_coordinate_h(1:3,point_index)-[0,3,0]');
    point_camera_distance(point_index)=norm( point_camera_coordinate(:,point_index));
    point_image(:,point_index)=point_camera_coordinate([1,2],point_index)/point_camera_coordinate(3,point_index);
    point_image_h(1:2,point_index)=point_image(:,point_index);
    point_image_h(3,point_index)=1;
end
% for point_index=1:num_point%����������ϵת��Ϊ���������ϵ�µ������Լ�ͼ������
%     point_camera_coordinate(:,point_index)=R*point_world_coordinate_h(:,point_index);
% %     distance(point_index)=norm(point_camera_coordinate(1:3,point_index));
%     point_image(:,point_index)=point_camera_coordinate(1:2,point_index)/point_camera_coordinate(3,point_index);
%     point_image_h(1:2,point_index)=point_image(:,point_index);
%     point_image_h(3,point_index)=1;
% end
[UU,DD,VV]=svd(point_world_coordinate_h);
[aa,dd]=qr(point_world_coordinate_h');
matrix2=eye(10)-aa(:,1:4)*aa(:,1:4)';
matrix_w=VV(:,4:end);%�ڴ˵õ���Ȩ�����Ա㽫���е����ֵ���и��õ����
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