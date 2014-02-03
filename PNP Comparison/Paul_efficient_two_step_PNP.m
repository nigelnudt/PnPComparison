function [RPaul,TPaul]=Paul_efficient_two_step_PNP(point_world_coor_h,point_image_h)
% author;Shuqiang Yang;
% %acording to "Efficient Linear Solution of Exterior Orientation"
% %reference author£ºPaul D. Fiore, Member, IEEE
%
num_point=length(point_image_h);
point_world_coordinate=point_world_coor_h(:,1:3);
[UU,DD,VV]=svd(point_world_coor_h');
matrix_w=VV(:,5:end);
for point_index=1:num_point
    matrix_for_l(1:2,point_index)=matrix_w(point_index,1)*point_image_h(point_index,1:2)';
    matrix_for_l(3:4,point_index)=matrix_w(point_index,2)*point_image_h(point_index,1:2)';
end
[Ul,Dl,Vl]=svd(matrix_for_l);
matrix_C=matrix_for_l*point_world_coor_h;
[UC,DC,VC]=svd(matrix_C);
l=point_world_coor_h*VC(:,4);
l=l/l(1);
for point_index=1:num_point
    point_image_h(point_index,:)=point_image_h(point_index,:)./norm(point_image_h(point_index,:));
    x3d_camera(point_index,:)=l(point_index)*point_image_h(point_index,:);
end
% distance_3p=[solution_confirmed_x1,solution_confirmed_x2,solution_confirmed_x3,solution_confirmed_x4];
% x3d_camera_p1=solution_confirmed_x1*x2d_h_normlized(1,:)./norm(x2d_h_normlized(1,:));
% x3d_camera_p2=solution_confirmed_x2*x2d_h_normlized(2,:)./norm(x2d_h_normlized(2,:));
% x3d_camera_p3=solution_confirmed_x3*x2d_h_normlized(3,:)./norm(x2d_h_normlized(3,:));
% x3d_camera_p4=solution_confirmed_x4*x2d_h_normlized(4,:)./norm(x2d_h_normlized(4,:));
% x3d_camera=[x3d_camera_p1;x3d_camera_p2;x3d_camera_p3;x3d_camera_p4];
% Xc3=x3d_camera;
sum_world=0;
sum_camera_world=0;
world_coor_mean=mean(point_world_coor_h(:,1:3),1);
world_coor_matrix=point_world_coor_h(:,1:3)-repmat(world_coor_mean,num_point,1);
camera_coor_mean=mean(x3d_camera,1);
camera_coor_matrix=x3d_camera-repmat(camera_coor_mean,num_point,1);
for point_index=1:num_point
    sum_camera_world=sum_camera_world+norm(camera_coor_matrix(point_index,:))*norm(world_coor_matrix(point_index,:));
    sum_world=sum_world+norm(norm(world_coor_matrix(point_index,:)))*norm(world_coor_matrix(point_index,:));
end
scale_camera_world=sum_camera_world/sum_world;
[UR,DR,VR]=svd(world_coor_matrix'*camera_coor_matrix);
% [UR1,DR1,VR1]=svd(DA(:,3:5)'*b_matrix1');
AA=det(UR)*det(VR);
if AA<0
    RPaul=VR*[1,0,0;0,1,0;0,0,-1]*UR';
else
    RPaul=VR*UR';
end
% diff=rotation_matrix-R(1:3,1:3);
% rotation_matrix1=VR1*UR1';
% rotation_matrix=rotation_matrix*[1,0,0;0,-1,0;0,0,1];
TPaul=1/scale_camera_world*camera_coor_mean'-RPaul*world_coor_mean';
% distance=distance/distance(1);
% point_camera_coordinate_ratio=point_camera_coordinate(3,:)/point_camera_coordinate(3,1)
pause(0.1)