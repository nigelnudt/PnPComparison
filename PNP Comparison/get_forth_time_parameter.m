function [aa8,aa6,aa4,aa2,aa0]=get_forth_time_parameter(xworld,yworld,zworld,R)
%通过本函数，将输入的三个点进行转换为一个有关医院四次方程，其中以第一个点的距离为所有方程的变量，通过Sylvester消元法将第二个以及第三个点的距离消去，
%xworld，yworld和zworld为三个点的1*3维坐标向量
num_point=3;
point_world_coordinate=[xworld;yworld;zworld];
point_world_coordinate_h=[point_world_coordinate;ones(1,num_point)];
point_camera_coordinate=zeros(4,num_point);
point_image=zeros(2,num_point);
point_image_h=zeros(3,num_point);
for point_index=1:num_point%将世界坐标系转化为摄像机坐标系下的坐标以及图像坐标
    point_camera_coordinate(:,point_index)=R*point_world_coordinate_h(:,point_index);
    point_image(:,point_index)=point_camera_coordinate(1:2,point_index)/point_camera_coordinate(3,point_index);
    point_image_h(1:2,point_index)=point_image(:,point_index);
    point_image_h(3,point_index)=1;
end
d12=norm([xworld(1)-xworld(2),yworld(1)-yworld(2),zworld(1)-zworld(2)]);
d23=norm([xworld(2)-xworld(3),yworld(2)-yworld(3),zworld(2)-zworld(3)]);
d13=norm([xworld(1)-xworld(3),yworld(1)-yworld(3),zworld(1)-zworld(3)]);
d121=norm(point_camera_coordinate(1:3,1)-point_camera_coordinate(1:3,2));
d231=norm(point_camera_coordinate(1:3,2)-point_camera_coordinate(1:3,3));
d131=norm(point_camera_coordinate(1:3,1)-point_camera_coordinate(1:3,3));
%接着计算三个点之间的夹角
theta12=acos(sum(point_image_h(:,1).*point_image_h(:,2))/norm(point_image_h(:,1))/norm(point_image_h(:,2)));
theta13=acos(sum(point_image_h(:,1).*point_image_h(:,3))/norm(point_image_h(:,1))/norm(point_image_h(:,3)));
theta23=acos(sum(point_image_h(:,2).*point_image_h(:,3))/norm(point_image_h(:,2))/norm(point_image_h(:,3)));
%然后得到各个系数
aa8= 16*cos(theta12)^4 - 64*cos(theta12)^3*cos(theta13)*cos(theta23) + 64*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 +...
    32*cos(theta12)^2*cos(theta13)^2+ 32*cos(theta12)^2*cos(theta23)^2 - 32*cos(theta12)^2 - 64*cos(theta12)*cos(theta13)^3*cos(theta23)- ...
    64*cos(theta12)*cos(theta13)*cos(theta23)^3 + 64*cos(theta12)*cos(theta13)*cos(theta23) + 16*cos(theta13)^4 + 32*cos(theta13)^2*cos(theta23)^2 -...
    32*cos(theta13)^2 + 16*cos(theta23)^4 - 32*cos(theta23)^2 + 16

aa6= - 64*d23^2*cos(theta12)^4*cos(theta13)^2 + 32*d23^2*cos(theta12)^4 + 128*d23^2*cos(theta12)^3*cos(theta13)^3*cos(theta23) - ...
    32*d23^2*cos(theta12)^3*cos(theta13)*cos(theta23) - 64*d23^2*cos(theta12)^2*cos(theta13)^4 - 128*d23^2*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 + ...
    128*d23^2*cos(theta12)^2*cos(theta13)^2 + 32*d23^2*cos(theta12)^2*cos(theta23)^2 - 64*d23^2*cos(theta12)^2 - 32*d23^2*cos(theta12)*cos(theta13)^3*cos(theta23) +...
    32*d23^2*cos(theta12)*cos(theta13)*cos(theta23)^3 + 32*d23^2*cos(theta12)*cos(theta13)*cos(theta23) + 32*d23^2*cos(theta13)^4+ ...
    32*d23^2*cos(theta13)^2*cos(theta23)^2 - 64*d23^2*cos(theta13)^2 - 32*d23^2*cos(theta23)^2 + 32*d23^2 - 32*d13^2*cos(theta12)^4 +...
    96*d13^2*cos(theta12)^3*cos(theta13)*cos(theta23) - 64*d13^2*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 - 32*d13^2*cos(theta12)^2*cos(theta13)^2 - ...
    64*d13^2*cos(theta12)^2*cos(theta23)^2 + 64*d13^2*cos(theta12)^2 + 32*d13^2*cos(theta12)*cos(theta13)^3*cos(theta23) +...
    96*d13^2*cos(theta12)*cos(theta13)*cos(theta23)^3 - 96*d13^2*cos(theta12)*cos(theta13)*cos(theta23) - 32*d13^2*cos(theta13)^2*cos(theta23)^2 + ...
    32*d13^2*cos(theta13)^2 - 32*d13^2*cos(theta23)^4 + 64*d13^2*cos(theta23)^2 - 32*d13^2  + 32*d12^2*cos(theta12)^3*cos(theta13)*cos(theta23) - ...
    64*d12^2*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 - 32*d12^2*cos(theta12)^2*cos(theta13)^2 - 32*d12^2*cos(theta12)^2*cos(theta23)^2+...
    32*d12^2*cos(theta12)^2 + 96*d12^2*cos(theta12)*cos(theta13)^3*cos(theta23) + 96*d12^2*cos(theta12)*cos(theta13)*cos(theta23)^3 - ...
    96*d12^2*cos(theta12)*cos(theta13)*cos(theta23) - 32*d12^2*cos(theta13)^4 - 64*d12^2*cos(theta13)^2*cos(theta23)^2+ 64*d12^2*cos(theta13)^2 -...
    32*d12^2*cos(theta23)^4 + 64*d12^2*cos(theta23)^2 - 32*d12^2 

aa4= 16*d23^4*cos(theta12)^4 - 32*d23^4*cos(theta12)^3*cos(theta13)*cos(theta23) + 48*d23^4*cos(theta12)^2*cos(theta13)^2 + ...
    16*d23^4*cos(theta12)^2*cos(theta23)^2 - 40*d23^4*cos(theta12)^2 - 32*d23^4*cos(theta12)*cos(theta13)^3*cos(theta23) + ...
    16*d23^4*cos(theta12)*cos(theta13)*cos(theta23) + 16*d23^4*cos(theta13)^4 + 16*d23^4*cos(theta13)^2*cos(theta23)^2 - ...
    40*d23^4*cos(theta13)^2 - 8*d23^4*cos(theta23)^2 + 24*d23^4 - 32*d13^2*d23^2*cos(theta12)^4 + ...
    64*d13^2*d23^2*cos(theta12)^3*cos(theta13)*cos(theta23)+ 64*d13^2*d23^2*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 - ...
    64*d13^2*d23^2*cos(theta12)^2*cos(theta13)^2 - 64*d13^2*d23^2*cos(theta12)^2*cos(theta23)^2 + 80*d13^2*d23^2*cos(theta12)^2 -...
    32*d13^2*d23^2*cos(theta12)*cos(theta13)^3*cos(theta23) - 32*d13^2*d23^2*cos(theta12)*cos(theta13)*cos(theta23)^3 - ...
    32*d13^2*d23^2*cos(theta12)*cos(theta13)*cos(theta23) + 48*d13^2*d23^2*cos(theta13)^2 + 48*d13^2*d23^2*cos(theta23)^2 - ...
    48*d13^2*d23^2 + 16*d13^4*cos(theta12)^4- 32*d13^4*cos(theta12)^3*cos(theta13)*cos(theta23) + 16*d13^4*cos(theta12)^2*cos(theta13)^2 + ...
    48*d13^4*cos(theta12)^2*cos(theta23)^2 - 40*d13^4*cos(theta12)^2 - 32*d13^4*cos(theta12)*cos(theta13)*cos(theta23)^3 + ...
    16*d13^4*cos(theta12)*cos(theta13)*cos(theta23) + 16*d13^4*cos(theta13)^2*cos(theta23)^2 - 8*d13^4*cos(theta13)^2+ 16*d13^4*cos(theta23)^4 - ...
    40*d13^4*cos(theta23)^2 + 24*d13^4  - 32*d12^2*d23^2*cos(theta12)^3*cos(theta13)*cos(theta23) + 64*d12^2*d23^2*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 -...
    64*d12^2*d23^2*cos(theta12)^2*cos(theta13)^2 + 48*d12^2*d23^2*cos(theta12)^2+ 64*d12^2*d23^2*cos(theta12)*cos(theta13)^3*cos(theta23) - ...
    32*d12^2*d23^2*cos(theta12)*cos(theta13)*cos(theta23)^3 - 32*d12^2*d23^2*cos(theta12)*cos(theta13)*cos(theta23) - 32*d12^2*d23^2*cos(theta13)^4 -...
    64*d12^2*d23^2*cos(theta13)^2*cos(theta23)^2+ 80*d12^2*d23^2*cos(theta13)^2 + 48*d12^2*d23^2*cos(theta23)^2 -...
    48*d12^2*d23^2 - 32*d12^2*d13^2*cos(theta12)^3*cos(theta13)*cos(theta23) + 64*d12^2*d13^2*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 + ...
    32*d12^2*d13^2*cos(theta12)^2*cos(theta23)^2 - 48*d12^2*d13^2*cos(theta12)^2 - 32*d12^2*d13^2*cos(theta12)*cos(theta13)^3*cos(theta23) - ...
    128*d12^2*d13^2*cos(theta12)*cos(theta13)*cos(theta23)^3+ 160*d12^2*d13^2*cos(theta12)*cos(theta13)*cos(theta23) + ...
    32*d12^2*d13^2*cos(theta13)^2*cos(theta23)^2 - 48*d12^2*d13^2*cos(theta13)^2 + 64*d12^2*d13^2*cos(theta23)^4 - 112*d12^2*d13^2*cos(theta23)^2 +...
    48*d12^2*d13^2 + 16*d12^4*cos(theta12)^2*cos(theta13)^2 + 16*d12^4*cos(theta12)^2*cos(theta23)^2 - 8*d12^4*cos(theta12)^2 - ...
    32*d12^4*cos(theta12)*cos(theta13)^3*cos(theta23) - 32*d12^4*cos(theta12)*cos(theta13)*cos(theta23)^3 + 16*d12^4*cos(theta12)*cos(theta13)*cos(theta23) +...
    16*d12^4*cos(theta13)^4 + 48*d12^4*cos(theta13)^2*cos(theta23)^2 - 40*d12^4*cos(theta13)^2 + 16*d12^4*cos(theta23)^4 - 40*d12^4*cos(theta23)^2 + 24*d12^4;

aa2=- 8*d23^6*cos(theta12)^2 + 8*d23^6*cos(theta12)*cos(theta13)*cos(theta23) - 8*d23^6*cos(theta13)^2 + 8*d23^6  - ...
    16*d13^2*d23^4*cos(theta12)^2*cos(theta23)^2 + 24*d13^2*d23^4*cos(theta12)^2 - 8*d13^2*d23^4*cos(theta12)*cos(theta13)*cos(theta23) +...
    16*d13^2*d23^4*cos(theta13)^2+ 8*d13^2*d23^4*cos(theta23)^2 - 24*d13^2*d23^4 - 16*d13^6*cos(theta12)^2*cos(theta23)^2 + ...
    8*d13^6*cos(theta12)^2 + 8*d13^6*cos(theta12)*cos(theta13)*cos(theta23) + 8*d13^6*cos(theta23)^2 - 8*d13^6  + ...
    32*d13^4*d23^2*cos(theta12)^2*cos(theta23)^2 - 24*d13^4*d23^2*cos(theta12)^2 - 8*d13^4*d23^2*cos(theta12)*cos(theta13)*cos(theta23) - ...
    8*d13^4*d23^2*cos(theta13)^2 - 16*d13^4*d23^2*cos(theta23)^2 + 24*d13^4*d23^2+ 16*d12^2*d23^4*cos(theta12)^2 - ...
    8*d12^2*d23^4*cos(theta12)*cos(theta13)*cos(theta23) - 16*d12^2*d23^4*cos(theta13)^2*cos(theta23)^2 + 24*d12^2*d23^4*cos(theta13)^2 +...
    8*d12^2*d23^4*cos(theta23)^2 - 24*d12^2*d23^4+ 16*d12^2*d13^4*cos(theta12)^2+ 32*d12^2*d13^4*cos(theta12)*cos(theta13)*cos(theta23)^3 - ...
    40*d12^2*d13^4*cos(theta12)*cos(theta13)*cos(theta23) - 16*d12^2*d13^4*cos(theta13)^2*cos(theta23)^2 + 8*d12^2*d13^4*cos(theta13)^2 -...
    32*d12^2*d13^4*cos(theta23)^4 + 56*d12^2*d13^4*cos(theta23)^2 - 24*d12^2*d13^4  - 32*d12^2*d13^2*d23^2*cos(theta12)^2 + ...
    32*d12^2*d13^2*d23^2*cos(theta12)*cos(theta13)*cos(theta23)^3 + 48*d12^2*d13^2*d23^2*cos(theta12)*cos(theta13)*cos(theta23) - ...
    32*d12^2*d13^2*d23^2*cos(theta13)^2- 64*d12^2*d13^2*d23^2*cos(theta23)^2 + 48*d12^2*d13^2*d23^2 - 16*d12^4*d13^2*cos(theta12)^2*cos(theta23)^2 +...
    8*d12^4*d13^2*cos(theta12)^2 + 32*d12^4*d13^2*cos(theta12)*cos(theta13)*cos(theta23)^3 - 40*d12^4*d13^2*cos(theta12)*cos(theta13)*cos(theta23) + ...
    16*d12^4*d13^2*cos(theta13)^2 - 32*d12^4*d13^2*cos(theta23)^4 + 56*d12^4*d13^2*cos(theta23)^2 - 24*d12^4*d13^2 - 8*d12^4*d23^2*cos(theta12)^2 - ...
    8*d12^4*d23^2*cos(theta12)*cos(theta13)*cos(theta23) + 32*d12^4*d23^2*cos(theta13)^2*cos(theta23)^2 - 24*d12^4*d23^2*cos(theta13)^2 - ...
    16*d12^4*d23^2*cos(theta23)^2+ 24*d12^4*d23^2+ 8*d12^6*cos(theta12)*cos(theta13)*cos(theta23) - 16*d12^6*cos(theta13)^2*cos(theta23)^2 +...
    8*d12^6*cos(theta13)^2 + 8*d12^6*cos(theta23)^2 - 8*d12^6 ;


aa0=d12^8 - 8*d12^6*d13^2*cos(theta23)^2 + 4*d12^6*d13^2 - 4*d12^6*d23^2 + 16*d12^4*d13^4*cos(theta23)^4 - 16*d12^4*d13^4*cos(theta23)^2 +...
    6*d12^4*d13^4 + 16*d12^4*d13^2*d23^2*cos(theta23)^2 - 12*d12^4*d13^2*d23^2  - 8*d12^2*d13^6*cos(theta23)^2+ 4*d12^2*d13^6 + ...
    16*d12^2*d13^4*d23^2*cos(theta23)^2 - 12*d12^2*d13^4*d23^2 - 4*d12^2*d23^6 + d13^8 - 4*d13^6*d23^2 - 4*d13^2*d23^6+ d23^8 + ...
    6*d13^4*d23^4- 8*d12^2*d13^2*d23^4*cos(theta23)^2 + 12*d12^2*d13^2*d23^4+ 6*d12^4*d23^4 ;