function [Rp3,tp3,Xc3,Rp4,tp4,Xc4]=LongQuan3Point4PointPnPv2(x3d_h,x2d_h_normlized)
%吃方程主要比较了利用3组分别解方程求出X后平均的解与原文的解，以及5点平均的解在不同性能下的比较，并估算计算时间。
%该方程返回的时候三点法和4点法返回的首先建立一个
% x1^2+x2^2-2*x1*x2*cos(theta12)-d12^2=0
% x1^2+x3^2-2*x1*x3*cos(theta13)-d13^2=0
% x2^2+x3^2-2*x2*x3*cos(theta23)-d23^2=0
%三个未知数，基于Sylvester约减

num_point=4;
xworld=x3d_h(1:4,1);
yworld=x3d_h(1:4,2);
zworld=x3d_h(1:4,3);
point_image_h=x2d_h_normlized;
d12=norm([xworld(1)-xworld(2),yworld(1)-yworld(2),zworld(1)-zworld(2)]);
d13=norm([xworld(1)-xworld(3),yworld(1)-yworld(3),zworld(1)-zworld(3)]);
d14=norm([xworld(1)-xworld(4),yworld(1)-yworld(4),zworld(1)-zworld(4)]);
d23=norm([xworld(2)-xworld(3),yworld(2)-yworld(3),zworld(2)-zworld(3)]);
d24=norm([xworld(2)-xworld(4),yworld(2)-yworld(4),zworld(2)-zworld(4)]);
d34=norm([xworld(3)-xworld(4),yworld(3)-yworld(4),zworld(3)-zworld(4)]);
% d121=norm(point_camera_coordinate(1:3,1)-point_camera_coordinate(1:3,2));
% d231=norm(point_camera_coordinate(1:3,2)-point_camera_coordinate(1:3,3));
% d131=norm(point_camera_coordinate(1:3,1)-point_camera_coordinate(1:3,3));
theta12=acos(sum(point_image_h(1,:).*point_image_h(2,:))/norm(point_image_h(1,:))/norm(point_image_h(2,:)));
theta13=acos(sum(point_image_h(1,:).*point_image_h(3,:))/norm(point_image_h(1,:))/norm(point_image_h(3,:)));
theta14=acos(sum(point_image_h(1,:).*point_image_h(4,:))/norm(point_image_h(1,:))/norm(point_image_h(4,:)));
theta23=acos(sum(point_image_h(2,:).*point_image_h(3,:))/norm(point_image_h(2,:))/norm(point_image_h(3,:)));
theta24=acos(sum(point_image_h(2,:).*point_image_h(4,:))/norm(point_image_h(2,:))/norm(point_image_h(4,:)));
theta34=acos(sum(point_image_h(3,:).*point_image_h(4,:))/norm(point_image_h(3,:))/norm(point_image_h(4,:)));
% [aa8,aa6,aa4,aa2,aa0]=get_forth_time_parameter(x3d_h(1:3,1),x3d_h(1:3,2),x3d_h(1:3,3),x2d_h_normlized(1:3,1),x2d_h_normlized(1:3,2));

% solution_for_x1=zeros(4,4);%4个组合，123,124,134,234，每个组合可以求解出四个解，所以是4*4的一个组合，每行代表每个组合的四个解%总共四个组合
% solution_for_x2=zeros(4,4);
% solution_for_x3=zeros(4,4);
% solution_for_x4=zeros(4,4);
% [aa8,aa6,aa4,aa2,aa0]=get_forth_time_parameter_zhangdata(x3d_h(1:3,1),x3d_h(1:3,2),x3d_h(1:3,3),x2d_h_normlized(1:3,1),x2d_h_normlized(1:3,2));
% coefficient_class8_group1=[aa8,aa6,aa4,aa2,aa0];
% solution_group11=roots(coefficient_class8_group1);
% solution_group1=real(solution_group11);
% 
% [aa8,aa6,aa4,aa2,aa0]=get_forth_time_parameter_zhangdata(x3d_h([1,2,4],1),x3d_h([1,2,4],2),x3d_h([1,2,4],3),x2d_h_normlized([1,2,4],1),x2d_h_normlized([1,2,4],2));
% coefficient_class8_group2=[aa8,aa6,aa4,aa2,aa0];
% solution_group22=roots(coefficient_class8_group2);
% solution_group2=real(solution_group22);
% 
% [aa8,aa6,aa4,aa2,aa0]=get_forth_time_parameter_zhangdata(x3d_h([1,3,4],1),x3d_h([1,3,4],2),x3d_h([1,3,4],3),x2d_h_normlized([1,3,4],1),x2d_h_normlized([1,3,4],2));
% coefficient_class8_group3=[aa8,aa6,aa4,aa2,aa0];
% solution_group33=roots(coefficient_class8_group3);
% solution_group3=real(solution_group33);
% 
% [aa8,aa6,aa4,aa2,aa0]=get_forth_time_parameter_zhangdata(x3d_h(2:4,1),x3d_h(2:4,2),x3d_h(2:4,3),x2d_h_normlized(2:4,1),x2d_h_normlized(2:4,2));
% coefficient_class8_group4=[aa8,aa6,aa4,aa2,aa0];
% solution_group44=roots(coefficient_class8_group4);
% solution_group4=real(solution_group44);
% 
% residual_total=zeros(4,4);
% for solution_index=1:length(solution_group1)%group1,2,3
%     if isreal(solution_group1(solution_index))
%         if solution_group1(solution_index)>0
%             solution_for_x1(1,solution_index)=sqrt(solution_group1(solution_index));
%             coefficient_for_x2=[1,-2*solution_for_x1(1,solution_index)*cos(theta12),solution_for_x1(1,solution_index)^2-d12^2];
%             possible_solution_for_x2=roots(coefficient_for_x2);
%             residual_2=zeros(1,2);
%             residual_3=zeros(1,2);
%             for index_distance2=1:2
%                 x1=solution_for_x1(1,solution_index);
%                 a4= 1;
%                 a3=- 4*x1*cos(theta13)*cos(theta23);
%                 a2=- 4*d13^2*cos(theta23)^2 + 2*d13^2 - 2*d23^2+ 4*x1^2*cos(theta13)^2 + 4*x1^2*cos(theta23)^2 - 2*x1^2 ;
%                 a1= - 4*x1^3*cos(theta13)*cos(theta23)+ 4*d23^2*x1*cos(theta13)*cos(theta23)+ 4*d13^2*x1*cos(theta13)*cos(theta23);
%                 a0=d13^4 - 2*d13^2*d23^2 - 2*d13^2*x1^2  + d23^4- 4*d23^2*x1^2*cos(theta13)^2 + 2*d23^2*x1^2  + x1^4;
%                 residual_2(index_distance2)=a4*possible_solution_for_x2(index_distance2)^4+a3*possible_solution_for_x2(index_distance2)^3+...
%                     a2*possible_solution_for_x2(index_distance2)^2+a1*possible_solution_for_x2(index_distance2)+a0;
%             end
%             if abs(residual_2(1))<abs(residual_2(2))
%                 solution_for_x2(1,solution_index)=possible_solution_for_x2(1);
%                 residual2=residual_2(1);
%             else
%                 solution_for_x2(1,solution_index)=possible_solution_for_x2(2);
%                 residual2=residual_2(2);
%             end
%             coefficient_for_x3=[1,-2*solution_for_x1(1,solution_index)*cos(theta13),solution_for_x1(1,solution_index)^2-d13^2];
%             possible_solution_for_x3=roots(coefficient_for_x3);
%             for index_distance3=1:2
%                 residual_3(index_distance3)=solution_for_x2(1,solution_index)^2+possible_solution_for_x3(index_distance3)^2-...
%                     2*solution_for_x2(1,solution_index)*possible_solution_for_x3(index_distance3)*cos(theta23)-d23^2;
%             end
%             if abs(residual_3(1))<abs(residual_3(2))
%                 solution_for_x3(1,solution_index)=possible_solution_for_x3(1);
%                 residual3=abs(residual_3(1));
%             else
%                 solution_for_x3(1,solution_index)=possible_solution_for_x3(2);
%                 residual3=abs(residual_3(2));
%             end
% %             solution_for3=roots(coefficient_for_3);
%             residual_total(1,solution_index)=residual2+residual3;
%         end
%     end
% end
% 
% for solution_index=1:length(solution_group2)%group1,2,4
%     if isreal(solution_group2(solution_index))
%         if solution_group2(solution_index)>0
%             solution_for_x1(2,solution_index)=sqrt(solution_group2(solution_index));
%             coefficient_for_x2=[1,-2*solution_for_x1(2,solution_index)*cos(theta12),solution_for_x1(2,solution_index)^2-d12^2];
%             possible_solution_for_x2=roots(coefficient_for_x2);
%             residual_2=zeros(1,2);
%             residual_4=zeros(1,2);
%             for index_distance2=1:2
%                 x1=solution_for_x1(2,solution_index);
%                 a4= 1;
%                 a3=- 4*x1*cos(theta14)*cos(theta24);
%                 a2=- 4*d14^2*cos(theta24)^2 + 2*d14^2 - 2*d24^2+ 4*x1^2*cos(theta14)^2 + 4*x1^2*cos(theta24)^2 - 2*x1^2 ;
%                 a1= - 4*x1^3*cos(theta14)*cos(theta24)+ 4*d24^2*x1*cos(theta14)*cos(theta24)+ 4*d14^2*x1*cos(theta14)*cos(theta24);
%                 a0=d14^4 - 2*d14^2*d24^2 - 2*d14^2*x1^2  + d24^4- 4*d24^2*x1^2*cos(theta14)^2 + 2*d24^2*x1^2  + x1^4;
%                 residual_2(index_distance2)=a4*possible_solution_for_x2(index_distance2)^4+a3*possible_solution_for_x2(index_distance2)^3+...
%                     a2*possible_solution_for_x2(index_distance2)^2+a1*possible_solution_for_x2(index_distance2)+a0;
%             end
%             if abs(residual_2(1))<abs(residual_2(2))
%                 solution_for_x2(2,solution_index)=possible_solution_for_x2(1);
%                 residual2=residual_2(1);
%             else
%                 solution_for_x2(2,solution_index)=possible_solution_for_x2(2);
%                 residual2=residual_2(2);
%             end
%             coefficient_for_x4=[1,-2*solution_for_x1(2,solution_index)*cos(theta14),solution_for_x1(2,solution_index)^2-d14^2];
%             possible_solution_for_x4=roots(coefficient_for_x4);
%             for index_distance4=1:2
%                 residual_4(index_distance4)=solution_for_x2(2,solution_index)^2+possible_solution_for_x4(index_distance4)^2-...
%                     2*solution_for_x2(2,solution_index)*possible_solution_for_x4(index_distance4)*cos(theta24)-d24^2;
%             end
%             if abs(residual_4(1))<abs(residual_4(2))
%                 solution_for_x4(2,solution_index)=possible_solution_for_x4(1);
%                 residual4=abs(residual_4(1));
%             else
%                 solution_for_x4(2,solution_index)=possible_solution_for_x4(2);
%                 residual4=abs(residual_4(2));
%             end
% %             solution_for4=roots(coefficient_for_4);
%             residual_total(2,solution_index)=residual2+residual4;
%         end
%     end
% end
% 
% for solution_index=1:length(solution_group3)%group1,3,4
%     if isreal(solution_group3(solution_index))
%         if solution_group3(solution_index)>0
%             solution_for_x1(3,solution_index)=sqrt(solution_group3(solution_index));
%             coefficient_for_x3=[1,-2*solution_for_x1(3,solution_index)*cos(theta13),solution_for_x1(3,solution_index)^2-d13^2];
%             possible_solution_for_x3=roots(coefficient_for_x3);
%             residual_3=zeros(1,2);
%             residual_4=zeros(1,2);
%             for index_distance3=1:2
%                 x1=solution_for_x1(3,solution_index);
%                 a4= 1;
%                 a3=- 4*x1*cos(theta14)*cos(theta34);
%                 a2=- 4*d14^2*cos(theta34)^2 + 2*d14^2 - 2*d34^2+ 4*x1^2*cos(theta14)^2 + 4*x1^2*cos(theta34)^2 - 2*x1^2 ;
%                 a1= - 4*x1^3*cos(theta14)*cos(theta34)+ 4*d34^2*x1*cos(theta14)*cos(theta34)+ 4*d14^2*x1*cos(theta14)*cos(theta34);
%                 a0=d14^4 - 2*d14^2*d34^2 - 2*d14^2*x1^2  + d34^4- 4*d34^2*x1^2*cos(theta14)^2 + 2*d34^2*x1^2  + x1^4;
%                 residual_3(index_distance3)=a4*possible_solution_for_x3(index_distance3)^4+a3*possible_solution_for_x3(index_distance3)^3+...
%                     a2*possible_solution_for_x3(index_distance3)^2+a1*possible_solution_for_x3(index_distance3)+a0;
%             end
%             if abs(residual_3(1))<abs(residual_3(2))
%                 solution_for_x3(3,solution_index)=possible_solution_for_x3(1);
%                 residual3=residual_3(1);
%             else
%                 solution_for_x3(3,solution_index)=possible_solution_for_x3(2);
%                 residual3=residual_3(2);
%             end
%             coefficient_for_x4=[1,-2*solution_for_x1(3,solution_index)*cos(theta14),solution_for_x1(3,solution_index)^2-d14^2];
%             possible_solution_for_x4=roots(coefficient_for_x4);
%             for index_distance4=1:2
%                 residual_4(index_distance4)=solution_for_x3(3,solution_index)^2+possible_solution_for_x4(index_distance4)^2-...
%                     2*solution_for_x3(3,solution_index)*possible_solution_for_x4(index_distance4)*cos(theta34)-d34^2;
%             end
%             if abs(residual_4(1))<abs(residual_4(2))
%                 solution_for_x4(3,solution_index)=possible_solution_for_x4(1);
%                 residual4=abs(residual_4(1));
%             else
%                 solution_for_x4(3,solution_index)=possible_solution_for_x4(2);
%                 residual4=abs(residual_4(2));
%             end
% %             solution_for4=roots(coefficient_for_4);
%             residual_total(3,solution_index)=residual3+residual4;
%         end
%     end
% end
% for solution_index=1:length(solution_group4)%group2,3,4
%     if isreal(solution_group4(solution_index))
%         if solution_group4(solution_index)>0
%             solution_for_x2(4,solution_index)=sqrt(solution_group4(solution_index));
%             coefficient_for_x3=[1,-2*solution_for_x2(4,solution_index)*cos(theta23),solution_for_x2(4,solution_index)^2-d23^2];
%             possible_solution_for_x3=roots(coefficient_for_x3);
%             residual_3=zeros(1,2);
%             residual_4=zeros(1,2);
%             for index_distance3=1:2
%                 x1=solution_for_x2(4,solution_index);
%                 a4= 1;
%                 a3=- 4*x1*cos(theta24)*cos(theta34);
%                 a2=- 4*d24^2*cos(theta34)^2 + 2*d24^2 - 2*d34^2+ 4*x1^2*cos(theta24)^2 + 4*x1^2*cos(theta34)^2 - 2*x1^2 ;
%                 a1= - 4*x1^3*cos(theta24)*cos(theta34)+ 4*d34^2*x1*cos(theta24)*cos(theta34)+ 4*d24^2*x1*cos(theta24)*cos(theta34);
%                 a0=d24^4 - 2*d24^2*d34^2 - 2*d24^2*x1^2  + d34^4- 4*d34^2*x1^2*cos(theta24)^2 + 2*d34^2*x1^2  + x1^4;
%                 residual_3(index_distance3)=a4*possible_solution_for_x3(index_distance3)^4+a3*possible_solution_for_x3(index_distance3)^3+...
%                     a2*possible_solution_for_x3(index_distance3)^2+a1*possible_solution_for_x3(index_distance3)+a0;
%             end
%             if abs(residual_3(1))<abs(residual_3(2))
%                 solution_for_x3(4,solution_index)=possible_solution_for_x3(1);
%                 residual3=residual_3(1);
%             else
%                 solution_for_x3(4,solution_index)=possible_solution_for_x3(2);
%                 residual3=residual_3(2);
%             end
%             coefficient_for_x4=[1,-2*solution_for_x2(4,solution_index)*cos(theta24),solution_for_x2(4,solution_index)^2-d24^2];
%             possible_solution_for_x4=roots(coefficient_for_x4);
%             for index_distance4=1:2
%                 residual_4(index_distance4)=solution_for_x3(4,solution_index)^2+possible_solution_for_x4(index_distance4)^2-...
%                     2*solution_for_x3(4,solution_index)*possible_solution_for_x4(index_distance4)*cos(theta34)-d34^2;
%             end
%             if abs(residual_4(1))<abs(residual_4(2))
%                 solution_for_x4(4,solution_index)=possible_solution_for_x4(1);
%                 residual4=abs(residual_4(1));
%             else
%                 solution_for_x4(4,solution_index)=possible_solution_for_x4(2);
%                 residual4=abs(residual_4(2));
%             end
% %             solution_for4=roots(coefficient_for_4);
%             residual_total(4,solution_index)=residual3+residual4;
%         end
%     end
% end
% min_distance=100000;
% for solution_index_group1=1:4
%     for solution_index_group2=1:4
%         for solution_index_group3=1:4
%             for solution_index_group4=1:4
%                 solution_set_1=...
%                     [solution_for_x1(1,solution_index_group1),solution_for_x1(2,solution_index_group2),solution_for_x1(3,solution_index_group3)];
%                 solution_set_2=...
%                     [solution_for_x2(1,solution_index_group1),solution_for_x2(2,solution_index_group2),solution_for_x2(4,solution_index_group4)];
%                 solution_set_3=...
%                     [solution_for_x3(1,solution_index_group1),solution_for_x3(3,solution_index_group3),solution_for_x3(4,solution_index_group4)];
%                 solution_set_4=...
%                     [solution_for_x4(2,solution_index_group2),solution_for_x4(3,solution_index_group3),solution_for_x4(4,solution_index_group4)];
%                 isall_valide=length(find([solution_set_1,solution_set_2,solution_set_3,solution_set_4]==0));
%                 var_this_set=sum([var(solution_set_1),var(solution_set_2),var(solution_set_3),var(solution_set_4)]);
%                 if isall_valide==0&&var_this_set<min_distance
%                     min_distance=var_this_set;
%                     best_set_index=[solution_index_group1,solution_index_group2,solution_index_group3,solution_index_group4];
%                 end
%             end
%         end
%     end
% end
% solution_confirmed_x1=mean([solution_for_x1(1,best_set_index(1)),solution_for_x1(2,best_set_index(2)),solution_for_x1(3,best_set_index(3))]);
% solution_confirmed_x2=mean([solution_for_x2(1,best_set_index(1)),solution_for_x2(2,best_set_index(2)),solution_for_x2(4,best_set_index(4))]);
% solution_confirmed_x3=mean([solution_for_x3(1,best_set_index(1)),solution_for_x3(3,best_set_index(3)),solution_for_x3(4,best_set_index(4))]);
% solution_confirmed_x4=mean([solution_for_x4(2,best_set_index(2)),solution_for_x4(3,best_set_index(3)),solution_for_x4(4,best_set_index(4))]);
% % distance_3p=[solution_confirmed_x1,solution_confirmed_x2,solution_confirmed_x3,solution_confirmed_x4];
% x3d_camera_p1=solution_confirmed_x1*x2d_h_normlized(1,:)./norm(x2d_h_normlized(1,:));
% x3d_camera_p2=solution_confirmed_x2*x2d_h_normlized(2,:)./norm(x2d_h_normlized(2,:));
% x3d_camera_p3=solution_confirmed_x3*x2d_h_normlized(3,:)./norm(x2d_h_normlized(3,:));
% x3d_camera_p4=solution_confirmed_x4*x2d_h_normlized(4,:)./norm(x2d_h_normlized(4,:));
% x3d_camera=[x3d_camera_p1;x3d_camera_p2;x3d_camera_p3;x3d_camera_p4];
% Xc3=x3d_camera;
% sum_world=0;
% sum_camera_world=0;
% world_coor_mean=mean(x3d_h(1:4,1:3),1);
% world_coor_matrix=x3d_h(1:4,1:3)-repmat(world_coor_mean,num_point,1);
% camera_coor_mean=mean(x3d_camera,1);
% camera_coor_matrix=x3d_camera-repmat(camera_coor_mean,num_point,1);
% for point_index=1:num_point
%     sum_camera_world=sum_camera_world+norm(camera_coor_matrix(point_index,:))*norm(camera_coor_matrix(point_index,:));
%     sum_world=sum_world+norm(world_coor_matrix(point_index,:))*norm(world_coor_matrix(point_index,:));
% end
% scale_camera_world=sum_camera_world/sum_world;
% [UR,DR,VR]=svd(world_coor_matrix'*camera_coor_matrix);
% % [UR1,DR1,VR1]=svd(DA(:,3:5)'*b_matrix1');
% AA=det(UR)*det(VR);
% if AA<0
%     Rp3=VR*[1,0,0;0,1,0;0,0,-1]*UR';
% else
%     Rp3=VR*UR';
% end
% % diff=rotation_matrix-R(1:3,1:3);
% % rotation_matrix1=VR1*UR1';
% % rotation_matrix=rotation_matrix*[1,0,0;0,-1,0;0,0,1];
% tp3=1/scale_camera_world*Rp3'*camera_coor_mean'-world_coor_mean';
% % Rp3=Rp3'
% tp3=Rp3*tp3;
% 

%%%%%%%%%%%
%下面是采用Long方法，运用四个点进行最小二乘得到深度值的求解方法
set_all_point=[1:num_point];
for prime_point_index=1:num_point
    set_prime_point=[prime_point_index];
    set_of_other=setdiff(set_all_point,set_prime_point);%从主点以外
    group_of_coef=(num_point-1)*(num_point-2)/2;
    coef_matrix=zeros(group_of_coef,5);
    coef_matrix_index=1;
    for group_index_point1=1:3%num_point-2
        for group_index_point2=group_index_point1+1:num_point-1
            triple_point_set=[prime_point_index,set_of_other(group_index_point1),set_of_other(group_index_point2)];
            %             coef_matrix(coef_matrix_index,:)=get_forth_time_parameter(DA(triple_point_set,3),DA(triple_point_set,4),DA(triple_point_set,5),RR)
            [aa8,aa6,aa4,aa2,aa0]=get_forth_time_parameter_zhangdata(x3d_h(triple_point_set,1),x3d_h(triple_point_set,2),x3d_h(triple_point_set,3),x2d_h_normlized(triple_point_set,1),x2d_h_normlized(triple_point_set,2));
%             [aa8,aa6,aa4,aa2,aa0]=get_forth_time_parameter_zhangdata(xworld(triple_point_set),yworld(triple_point_set),zworld(triple_point_set),point_image_h(triple_point_set),point_image_h(triple_point_set));
%             point_camera_distance(prime_point_index)=point_camera_distance(prime_point_index)*10;
%             residual(coef_matrix_index)=(aa8*point_camera_distance(prime_point_index)^8+aa6*point_camera_distance(prime_point_index)^6+...
%                 aa4*point_camera_distance(prime_point_index)^4+aa2*point_camera_distance(prime_point_index)^2+aa0)/aa0;
            
            coef_matrix(coef_matrix_index,:)=[aa0,aa2,aa4,aa6,aa8]/aa0;
            coef_matrix_index=coef_matrix_index+1;
        end
    end
    [U_5p1,d_5p1,V_5p1]=svd(coef_matrix);
    depth1(prime_point_index)=mean(sqrt(V_5p1(2:5,5)./V_5p1(1:4,5)));
%     scale=1/5;
    scale=1;
    scale_matrix=[scale^8,0,0,0,0;0,scale^6,0,0,0;0,0,scale^4,0,0;0,0,0,scale^2,0;0,0,0,0,1];
    coef_matrix=coef_matrix*scale_matrix;
    [U_5p,d_5p,V_5p]=svd(coef_matrix);
    depth(prime_point_index)=mean(sqrt(V_5p(2:5,5)./V_5p(1:4,5)))/scale;
    
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
    scale2=1/(ratio_nameda_rio*v4(1)+v5(1));
    t5=scale2*(ratio_nameda_rio*v4+v5);
%     x1=mean(sqrt(t5(2:5)./t5(1:4)))/scale;
    distance_MS1(prime_point_index)=mean(sqrt(t5(2:5)./t5(1:4)))/scale;
    distance_MS2(prime_point_index)=sqrt(t5(2)./t5(1))/scale; 

end

x3d_camera4_p1=distance_MS2(1)*x2d_h_normlized(1,:)./norm(x2d_h_normlized(1,:));
x3d_camera4_p2=distance_MS2(2)*x2d_h_normlized(2,:)./norm(x2d_h_normlized(2,:));
x3d_camera4_p3=distance_MS2(3)*x2d_h_normlized(3,:)./norm(x2d_h_normlized(3,:));
x3d_camera4_p4=distance_MS2(4)*x2d_h_normlized(4,:)./norm(x2d_h_normlized(4,:));
x3d_camera4=[x3d_camera4_p1;x3d_camera4_p2;x3d_camera4_p3;x3d_camera4_p4];
Xc4=x3d_camera4;
sum_world4=0;    
sum_camera_world4=0;
world_coor_mean4=mean(x3d_h(1:4,1:3),1);
world_coor_matrix4=x3d_h(1:4,1:3)-repmat(world_coor_mean4,num_point,1);
camera_coor_mean4=mean(x3d_camera4,1);
camera_coor_matrix4=x3d_camera4-repmat(camera_coor_mean4,num_point,1);
for point_index=1:num_point
    sum_camera_world4=sum_camera_world4+norm(camera_coor_matrix4(point_index,:))*norm(world_coor_matrix4(point_index,:));
    sum_world4=sum_world4+norm(norm(world_coor_matrix4(point_index,:)))*norm(world_coor_matrix4(point_index,:));
end
scale_camera_world4=sum_camera_world4/sum_world4;
[UR4,DR4,VR4]=svd(world_coor_matrix4'*camera_coor_matrix4);
% [UR1,DR1,VR1]=svd(DA(:,3:5)'*b_matrix1');
AA4=det(UR4)*det(VR4);
if AA4<0
    Rp4=VR4*[1,0,0;0,1,0;0,0,-1]*UR4';
else
    Rp4=VR4*UR4';
end
% diff=rotation_matrix-R(1:3,1:3);
% rotation_matrix1=VR1*UR1';
% rotation_matrix=rotation_matrix*[1,0,0;0,-1,0;0,0,1];
tp4=1/scale_camera_world4*Rp4'*camera_coor_mean4'-world_coor_mean4';
% Rp4=Rp4';
tp4=Rp4*tp4;
Rp3=Rp4;
tp3=tp4;
Xc3=Xc4;
% disp length(solution_group1) length(solution_group2) length(solution_group3) length(solution_group4)
% tic
% Rp3=1;
% tp3=1;
% Xc3=1;
% Rp4=1;
% tp4=1;
% Xc4=1;

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



