close all
clear
%线性求解摄像机的位姿参数，参考文献
%首先介绍四点法，通过任选三个点构造一个一元四次方程，将这些方程的
%Linear N-Point Camera Pose Determination，Long Quan and Zhongdan Lan
addpath 'D:\计算机视觉测量\张志龙老师合作\EPnP_matlab\data';
addpath 'C:\Users\nigel\Documents\Multiprecision Computing Toolbox'
% clc;
% clear all;
%此处将张志龙老师的仿真数据利用Long Quan的线性方法进行求解，看两者之间的误差
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cFileName = 'D:\计算机视觉测量\张志龙老师合作\N点位姿数据 2013-09-13\Pose_PointPattern1.dat';
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
% num_point=5;
% R=return_Rt_matrix(alpha,beta,gamma,tx,ty,tz);
% xworld=400*rand(1,num_point)-200;
% yworld=400*rand(1,num_point)-200;
% zworld=400*rand(1,num_point)-200;
num_point=10;
point_world_coordinate_h=[DA(:,3:5)';ones(1,num_point)];
R=[RR,[0,3,0]'];
for point_index=1:num_point%将世界坐标系转化为摄像机坐标系下的坐标以及图像坐标
    point_camera_coordinate(:,point_index)=RR*(point_world_coordinate_h(1:3,point_index)-[0,3,0]');
    point_camera_distance(point_index)=norm( point_camera_coordinate(:,point_index));
    point_image(:,point_index)=point_camera_coordinate([1,2],point_index)/point_camera_coordinate(3,point_index);
    point_image_h(1:2,point_index)=point_image(:,point_index);
    point_image_h(3,point_index)=1;
end
figure(1)
plot3( DA(:,3),DA(:,4),DA(:,5), 'o'); %  xlabel('xw');   ylabel('zw');    zlabel('yw');  title('物点分布');
hold on
% axis([-1,1,-1,1,0,5])
M1 = RR* ( DA(:,3:5)' - repmat( tt, 1, num ) );
% figure(2)
plot3(M1(1,:),M1(2,:),M1(3,:),'*');grid on;   xlabel('xc');   ylabel('yc');    zlabel('zc');  title('物点分布');
axis([-5,5,-5,5,-5,5])
for i=1:num,
    M1(1,i) = M1(1,i) / M1(3,i);
    M1(2,i) = M1(2,i) / M1(3,i);
    M1(3,i) = 0.0;
end
set_all_point=[1,2,3,4,5,6,7,8,9,10];
num_point=length(DA);
scale=5;
scale_matrix=[scale^8,0,0,0,0;0,scale^6,0,0,0;0,0,scale^4,0,0;0,0,0,scale^2,0;0,0,0,0,1];
for point_index1=1:9
    for point_index2=point_index1+1:10
        residual_original(point_index1,point_index2)=point_camera_distance(point_index1)^2+point_camera_distance(point_index2)^2-...
            2*point_camera_distance(point_index1)*point_camera_distance(point_index2)*sum(point_image_h(:,point_index1).*point_image_h(:,point_index2))/...
            norm(point_image_h(:,point_index1))/norm(point_image_h(:,point_index2))-(norm(point_camera_coordinate(:,point_index1)-point_camera_coordinate(:,point_index2)))^2;
    end
end

for prime_point_index=1:num_point
    set_prime_point=[prime_point_index];
    set_of_other=setdiff(set_all_point,set_prime_point);
    group_of_coef=(num_point-1)*(num_point-2)/2;
    coef_matrix=zeros(group_of_coef,5);
    coef_matrix_index=1;
    for group_index_point1=1:3%num_point-2
        for group_index_point2=group_index_point1+1:num_point-1
            triple_point_set=[prime_point_index,set_of_other(group_index_point1),set_of_other(group_index_point2)];
            %             coef_matrix(coef_matrix_index,:)=get_forth_time_parameter(DA(triple_point_set,3),DA(triple_point_set,4),DA(triple_point_set,5),RR)
            [aa8,aa6,aa4,aa2,aa0]=get_forth_time_parameter_zhangdata(DA(triple_point_set,3),DA(triple_point_set,4),DA(triple_point_set,5),DA(triple_point_set,1),DA(triple_point_set,2));
%             point_camera_distance(prime_point_index)=point_camera_distance(prime_point_index)*10;
            residual(coef_matrix_index)=(aa8*point_camera_distance(prime_point_index)^8+aa6*point_camera_distance(prime_point_index)^6+...
                aa4*point_camera_distance(prime_point_index)^4+aa2*point_camera_distance(prime_point_index)^2+aa0)/aa0;
            
            coef_matrix(coef_matrix_index,:)=[aa0,aa2,aa4,aa6,aa8]/aa0;
            coef_matrix_index=coef_matrix_index+1;
        end
    end
    %     coef_matrix=coef_matrix*scale_matrix;
    [U_5p,d_5p,V_5p]=svd(coef_matrix(1:36,:));
    depth(prime_point_index)=sqrt(V_5p(2,5)/V_5p(1,5));
    
    v4=V_5p(:,4);
    v5=V_5p(:,5);
    ijkl_index_matrix=[4,2,3,3;2,0,1,1;4,1,3,2;4,0,2,2;3,1,2,2;4,0,3,1;3,0,2,1;]+1;
    ijkl_index_matrix_extended=[ijkl_index_matrix];
    matrix_B_for_null_vector=zeros(7,3);
    for iter_index=1:7
        matrix_B_for_null_vector(iter_index,1)=v4(ijkl_index_matrix_extended(iter_index,1))*v4(ijkl_index_matrix_extended(iter_index,2))-v4(ijkl_index_matrix_extended(iter_index,3))*v4(ijkl_index_matrix_extended(iter_index,4));
        matrix_B_for_null_vector(iter_index,2)=v4(ijkl_index_matrix_extended(iter_index,1))*v5(ijkl_index_matrix_extended(iter_index,2))+v5(ijkl_index_matrix_extended(iter_index,1))*v4(ijkl_index_matrix_extended(iter_index,2))-...
            (v4(ijkl_index_matrix_extended(iter_index,3))*v5(ijkl_index_matrix_extended(iter_index,4))+v5(ijkl_index_matrix_extended(iter_index,3))*v4(ijkl_index_matrix_extended(iter_index,4)));
        matrix_B_for_null_vector(iter_index,3)=v5(ijkl_index_matrix_extended(iter_index,1))*v5(ijkl_index_matrix_extended(iter_index,2))-v5(ijkl_index_matrix_extended(iter_index,3))*v5(ijkl_index_matrix_extended(iter_index,4));
    end
    [UB,dB,Vb]=svd(matrix_B_for_null_vector);
    % for iter_index=1:21
    %     matrix_B_for_null_vector(iter_index,1)=v4(ijkl_index_matrix_extended(iter_index,1))*v4(ijkl_index_matrix_extended(iter_index,2))-v4(ijkl_index_matrix_extended(iter_index,3))*v4(ijkl_index_matrix_extended(iter_index,4));
    %     matrix_B_for_null_vector(iter_index,2)=v4(ijkl_index_matrix_extended(iter_index,1))*v5(ijkl_index_matrix_extended(iter_index,2))+v5(ijkl_index_matrix_extended(iter_index,1))*v4(ijkl_index_matrix_extended(iter_index,2))-...
    %         (v4(ijkl_index_matrix_extended(iter_index,3))*v5(ijkl_index_matrix_extended(iter_index,4))+v5(ijkl_index_matrix_extended(iter_index,4))*v4(ijkl_index_matrix_extended(iter_index,3)));
    %     matrix_B_for_null_vector(iter_index,3)=v5(ijkl_index_matrix_extended(iter_index,1))*v5(ijkl_index_matrix_extended(iter_index,2))-v5(ijkl_index_matrix_extended(iter_index,3))*v5(ijkl_index_matrix_extended(iter_index,4));
    % end
    % [UB1,dB1,Vb1]=svd(matrix_B_for_null_vector);
    ratio_nameda_rio=Vb(1,3)/Vb(2,3);
    scale2=1/(ratio_nameda_rio*v4(1)+v5(1))
    t5=scale2*(ratio_nameda_rio*v4+v5)
    x1=t5(2);
    
end
tic

% [aa81,aa61,aa41,aa21,aa01]=get_forth_time_parameter(xworld(1:3),yworld(1:3),zworld(1:3),R);
% [aa82,aa62,aa42,aa22,aa02]=get_forth_time_parameter(xworld([1,2,4]),yworld([1,2,4]),zworld([1,2,4]),R);
% [aa83,aa63,aa43,aa23,aa03]=get_forth_time_parameter(xworld([1,3,4]),yworld([1,3,4]),zworld([1,3,4]),R);
% [aa84,aa64,aa44,aa24,aa04]=get_forth_time_parameter(xworld([1,2,5]),yworld([1,2,5]),zworld([1,2,5]),R);
% [aa85,aa65,aa45,aa25,aa05]=get_forth_time_parameter(xworld([1,4,5]),yworld([1,4,5]),zworld([1,4,5]),R);
% [aa86,aa66,aa46,aa26,aa06]=get_forth_time_parameter(xworld([1,3,5]),yworld([1,3,5]),zworld([1,3,5]),R);
% solution1=roots([aa81,aa61,aa41,aa21,aa01]);
% solution2=roots([aa82,aa62,aa42,aa22,aa02]);
% solution3=roots([aa83,aa63,aa43,aa23,aa03]);
% solution10=roots([aa81,aa61*scale^2,aa41*scale^4,aa21*scale^6,aa01*scale^8]);
% solution20=roots([aa82,aa62*scale^2,aa42*scale^4,aa22*scale^6,aa02*scale^8]);
% solution30=roots([aa83,aa63*scale^2,aa43*scale^4,aa23*scale^6,aa03*scale^8]);
% matrix_for_distance_from_equations=[[aa01,aa21,aa41,aa61,aa81];[aa02,aa22,aa42,aa62,aa82];[aa03,aa23,aa43,aa63,aa83]];
% matrix_for_distance_from_equations_5p=[[aa01,aa21,aa41,aa61,aa81];[aa02,aa22,aa42,aa62,aa82];[aa03,aa23,aa43,aa63,aa83];...
%     [aa04,aa24,aa44,aa64,aa84];[aa05,aa25,aa45,aa65,aa85];[aa06,aa26,aa46,aa66,aa86]];
%为了改善数据的病态性，将x=scalex带入原方程，以便得到比较优良的解，相应的需要将估计的尺度信息放入里面，在此处我们是否可以写出一篇文章，也就是说尺度规范化在摄像机位姿线性求解中的影响

% scale=200;
% scale_matrix=[scale^8,0,0,0,0;0,scale^6,0,0,0;0,0,scale^4,0,0;0,0,0,scale^2,0;0,0,0,0,1];
% matrix_for_distance_from_equations_data_normlize_5p=coef_matrix*scale_matrix;
% matrix_for_distance_from_equations_data_normlize_5p=matrix_for_distance_from_equations_5p*scale_matrix;
% matrix_for_distance_from_equations_data_normlize(1,:)=matrix_for_distance_from_equations_data_normlize_1(1,:)/norm(matrix_for_distance_from_equations_data_normlize_1(1,:));
% matrix_for_distance_from_equations_data_normlize(2,:)=matrix_for_distance_from_equations_data_normlize_1(2,:)/norm(matrix_for_distance_from_equations_data_normlize_1(2,:));
% matrix_for_distance_from_equations_data_normlize(3,:)=matrix_for_distance_from_equations_data_normlize_1(3,:)/norm(matrix_for_distance_from_equations_data_normlize_1(3,:));
% % matrix_for_distance_from_equations_data_normlize=matrix_for_distance_from_equations_data_normlize_1
% [U,d,V]=svd(matrix_for_distance_from_equations_data_normlize);

% [U1,d1,V1]=svd(matrix_for_distance_from_equations_data_normlize,0);



