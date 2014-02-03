function [aa8,aa6,aa4,aa2,aa0]=get_forth_time_parameter_zhangdata(xworld,yworld,zworld,u_image,v_image)


%ͨ������������������������ͼ�������Լ������������ת��Ϊһ���й�һԪ�Ĵη��̣������Ե�һ����ľ���Ϊ���з��̵ı�����ͨ��Sylvester��Ԫ�����ڶ����Լ���������ľ�����ȥ��
%xworld��yworld��zworldΪ�������1*3ά��������
%u_image��v_image��ʾ��������ϵ�µ�����
num_point=3;
% point_world_coordinate=[xworld;yworld;zworld];
% point_world_coordinate_h=[point_world_coordinate;ones(1,num_point)];
% point_camera_coordinate=zeros(4,num_point);
% point_image=[u_image;v_image];'
% point_image_h=[u_image';v_image';ones(1,num_point)];
point_image_h=[u_image';v_image';ones(1,num_point)];
% for point_index=1:num_point%����������ϵת��Ϊ���������ϵ�µ������Լ�ͼ������
%     point_camera_coordinate(:,point_index)=RR*point_world_coordinate_h(:,point_index);
%     point_image(:,point_index)=point_camera_coordinate(1:2,point_index)/point_camera_coordinate(3,point_index);
%     point_image_h(1:2,point_index)=point_image(:,point_index);
%     point_image_h(3,point_index)=1;
% end
% d12=mp(norm([xworld(1)-xworld(2),yworld(1)-yworld(2),zworld(1)-zworld(2)]));
% d23=mp(norm([xworld(2)-xworld(3),yworld(2)-yworld(3),zworld(2)-zworld(3)]));
% d13=mp(norm([xworld(1)-xworld(3),yworld(1)-yworld(3),zworld(1)-zworld(3)]));
d12=norm([xworld(1)-xworld(2),yworld(1)-yworld(2),zworld(1)-zworld(2)]);
d23=norm([xworld(2)-xworld(3),yworld(2)-yworld(3),zworld(2)-zworld(3)]);
d13=norm([xworld(1)-xworld(3),yworld(1)-yworld(3),zworld(1)-zworld(3)]);

% d12=norm(point_camera_coordinate(1:3,1)-point_camera_coordinate(1:3,2));
% d23=norm(point_camera_coordinate(1:3,2)-point_camera_coordinate(1:3,3));
% d13=norm(point_camera_coordinate(1:3,1)-point_camera_coordinate(1:3,3));
%���ż���������֮��ļн�
% theta12=mp(sum(point_image_h(:,1).*point_image_h(:,2))/norm(point_image_h(:,1))/norm(point_image_h(:,2)));
% theta13=mp(sum(point_image_h(:,1).*point_image_h(:,3))/norm(point_image_h(:,1))/norm(point_image_h(:,3)));
% theta23=mp(sum(point_image_h(:,2).*point_image_h(:,3))/norm(point_image_h(:,2))/norm(point_image_h(:,3)));
theta12=sum(point_image_h(:,1).*point_image_h(:,2))/norm(point_image_h(:,1))/norm(point_image_h(:,2));
theta13=sum(point_image_h(:,1).*point_image_h(:,3))/norm(point_image_h(:,1))/norm(point_image_h(:,3));
theta23=sum(point_image_h(:,2).*point_image_h(:,3))/norm(point_image_h(:,2))/norm(point_image_h(:,3));
%Ȼ��õ�����ϵ��
% aa8= 16*cos(theta12)^4 - 64*cos(theta12)^3*cos(theta13)*cos(theta23) + 64*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 +...
%     32*cos(theta12)^2*cos(theta13)^2+ 32*cos(theta12)^2*cos(theta23)^2 - 32*cos(theta12)^2 - 64*cos(theta12)*cos(theta13)^3*cos(theta23)- ...
%     64*cos(theta12)*cos(theta13)*cos(theta23)^3 + 64*cos(theta12)*cos(theta13)*cos(theta23) + 16*cos(theta13)^4 + 32*cos(theta13)^2*cos(theta23)^2 -...
%     32*cos(theta13)^2 + 16*cos(theta23)^4 - 32*cos(theta23)^2 + 16
% 
% aa6= - 64*d23^2*cos(theta12)^4*cos(theta13)^2 + 32*d23^2*cos(theta12)^4 + 128*d23^2*cos(theta12)^3*cos(theta13)^3*cos(theta23) - ...
%     32*d23^2*cos(theta12)^3*cos(theta13)*cos(theta23) - 64*d23^2*cos(theta12)^2*cos(theta13)^4 - 128*d23^2*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 + ...
%     128*d23^2*cos(theta12)^2*cos(theta13)^2 + 32*d23^2*cos(theta12)^2*cos(theta23)^2 - 64*d23^2*cos(theta12)^2 - 32*d23^2*cos(theta12)*cos(theta13)^3*cos(theta23) +...
%     32*d23^2*cos(theta12)*cos(theta13)*cos(theta23)^3 + 32*d23^2*cos(theta12)*cos(theta13)*cos(theta23) + 32*d23^2*cos(theta13)^4+ ...
%     32*d23^2*cos(theta13)^2*cos(theta23)^2 - 64*d23^2*cos(theta13)^2 - 32*d23^2*cos(theta23)^2 + 32*d23^2 - 32*d13^2*cos(theta12)^4 +...
%     96*d13^2*cos(theta12)^3*cos(theta13)*cos(theta23) - 64*d13^2*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 - 32*d13^2*cos(theta12)^2*cos(theta13)^2 - ...
%     64*d13^2*cos(theta12)^2*cos(theta23)^2 + 64*d13^2*cos(theta12)^2 + 32*d13^2*cos(theta12)*cos(theta13)^3*cos(theta23) +...
%     96*d13^2*cos(theta12)*cos(theta13)*cos(theta23)^3 - 96*d13^2*cos(theta12)*cos(theta13)*cos(theta23) - 32*d13^2*cos(theta13)^2*cos(theta23)^2 + ...
%     32*d13^2*cos(theta13)^2 - 32*d13^2*cos(theta23)^4 + 64*d13^2*cos(theta23)^2 - 32*d13^2  + 32*d12^2*cos(theta12)^3*cos(theta13)*cos(theta23) - ...
%     64*d12^2*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 - 32*d12^2*cos(theta12)^2*cos(theta13)^2 - 32*d12^2*cos(theta12)^2*cos(theta23)^2+...
%     32*d12^2*cos(theta12)^2 + 96*d12^2*cos(theta12)*cos(theta13)^3*cos(theta23) + 96*d12^2*cos(theta12)*cos(theta13)*cos(theta23)^3 - ...
%     96*d12^2*cos(theta12)*cos(theta13)*cos(theta23) - 32*d12^2*cos(theta13)^4 - 64*d12^2*cos(theta13)^2*cos(theta23)^2+ 64*d12^2*cos(theta13)^2 -...
%     32*d12^2*cos(theta23)^4 + 64*d12^2*cos(theta23)^2 - 32*d12^2 
% 
% aa4= 16*d23^4*cos(theta12)^4 - 32*d23^4*cos(theta12)^3*cos(theta13)*cos(theta23) + 48*d23^4*cos(theta12)^2*cos(theta13)^2 + ...
%     16*d23^4*cos(theta12)^2*cos(theta23)^2 - 40*d23^4*cos(theta12)^2 - 32*d23^4*cos(theta12)*cos(theta13)^3*cos(theta23) + ...
%     16*d23^4*cos(theta12)*cos(theta13)*cos(theta23) + 16*d23^4*cos(theta13)^4 + 16*d23^4*cos(theta13)^2*cos(theta23)^2 - ...
%     40*d23^4*cos(theta13)^2 - 8*d23^4*cos(theta23)^2 + 24*d23^4 - 32*d13^2*d23^2*cos(theta12)^4 + ...
%     64*d13^2*d23^2*cos(theta12)^3*cos(theta13)*cos(theta23)+ 64*d13^2*d23^2*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 - ...
%     64*d13^2*d23^2*cos(theta12)^2*cos(theta13)^2 - 64*d13^2*d23^2*cos(theta12)^2*cos(theta23)^2 + 80*d13^2*d23^2*cos(theta12)^2 -...
%     32*d13^2*d23^2*cos(theta12)*cos(theta13)^3*cos(theta23) - 32*d13^2*d23^2*cos(theta12)*cos(theta13)*cos(theta23)^3 - ...
%     32*d13^2*d23^2*cos(theta12)*cos(theta13)*cos(theta23) + 48*d13^2*d23^2*cos(theta13)^2 + 48*d13^2*d23^2*cos(theta23)^2 - ...
%     48*d13^2*d23^2 + 16*d13^4*cos(theta12)^4- 32*d13^4*cos(theta12)^3*cos(theta13)*cos(theta23) + 16*d13^4*cos(theta12)^2*cos(theta13)^2 + ...
%     48*d13^4*cos(theta12)^2*cos(theta23)^2 - 40*d13^4*cos(theta12)^2 - 32*d13^4*cos(theta12)*cos(theta13)*cos(theta23)^3 + ...
%     16*d13^4*cos(theta12)*cos(theta13)*cos(theta23) + 16*d13^4*cos(theta13)^2*cos(theta23)^2 - 8*d13^4*cos(theta13)^2+ 16*d13^4*cos(theta23)^4 - ...
%     40*d13^4*cos(theta23)^2 + 24*d13^4  - 32*d12^2*d23^2*cos(theta12)^3*cos(theta13)*cos(theta23) + 64*d12^2*d23^2*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 -...
%     64*d12^2*d23^2*cos(theta12)^2*cos(theta13)^2 + 48*d12^2*d23^2*cos(theta12)^2+ 64*d12^2*d23^2*cos(theta12)*cos(theta13)^3*cos(theta23) - ...
%     32*d12^2*d23^2*cos(theta12)*cos(theta13)*cos(theta23)^3 - 32*d12^2*d23^2*cos(theta12)*cos(theta13)*cos(theta23) - 32*d12^2*d23^2*cos(theta13)^4 -...
%     64*d12^2*d23^2*cos(theta13)^2*cos(theta23)^2+ 80*d12^2*d23^2*cos(theta13)^2 + 48*d12^2*d23^2*cos(theta23)^2 -...
%     48*d12^2*d23^2 - 32*d12^2*d13^2*cos(theta12)^3*cos(theta13)*cos(theta23) + 64*d12^2*d13^2*cos(theta12)^2*cos(theta13)^2*cos(theta23)^2 + ...
%     32*d12^2*d13^2*cos(theta12)^2*cos(theta23)^2 - 48*d12^2*d13^2*cos(theta12)^2 - 32*d12^2*d13^2*cos(theta12)*cos(theta13)^3*cos(theta23) - ...
%     128*d12^2*d13^2*cos(theta12)*cos(theta13)*cos(theta23)^3+ 160*d12^2*d13^2*cos(theta12)*cos(theta13)*cos(theta23) + ...
%     32*d12^2*d13^2*cos(theta13)^2*cos(theta23)^2 - 48*d12^2*d13^2*cos(theta13)^2 + 64*d12^2*d13^2*cos(theta23)^4 - 112*d12^2*d13^2*cos(theta23)^2 +...
%     48*d12^2*d13^2 + 16*d12^4*cos(theta12)^2*cos(theta13)^2 + 16*d12^4*cos(theta12)^2*cos(theta23)^2 - 8*d12^4*cos(theta12)^2 - ...
%     32*d12^4*cos(theta12)*cos(theta13)^3*cos(theta23) - 32*d12^4*cos(theta12)*cos(theta13)*cos(theta23)^3 + 16*d12^4*cos(theta12)*cos(theta13)*cos(theta23) +...
%     16*d12^4*cos(theta13)^4 + 48*d12^4*cos(theta13)^2*cos(theta23)^2 - 40*d12^4*cos(theta13)^2 + 16*d12^4*cos(theta23)^4 - 40*d12^4*cos(theta23)^2 + 24*d12^4;
% 
% aa2=- 8*d23^6*cos(theta12)^2 + 8*d23^6*cos(theta12)*cos(theta13)*cos(theta23) - 8*d23^6*cos(theta13)^2 + 8*d23^6  - ...
%     16*d13^2*d23^4*cos(theta12)^2*cos(theta23)^2 + 24*d13^2*d23^4*cos(theta12)^2 - 8*d13^2*d23^4*cos(theta12)*cos(theta13)*cos(theta23) +...
%     16*d13^2*d23^4*cos(theta13)^2+ 8*d13^2*d23^4*cos(theta23)^2 - 24*d13^2*d23^4 - 16*d13^6*cos(theta12)^2*cos(theta23)^2 + ...
%     8*d13^6*cos(theta12)^2 + 8*d13^6*cos(theta12)*cos(theta13)*cos(theta23) + 8*d13^6*cos(theta23)^2 - 8*d13^6  + ...
%     32*d13^4*d23^2*cos(theta12)^2*cos(theta23)^2 - 24*d13^4*d23^2*cos(theta12)^2 - 8*d13^4*d23^2*cos(theta12)*cos(theta13)*cos(theta23) - ...
%     8*d13^4*d23^2*cos(theta13)^2 - 16*d13^4*d23^2*cos(theta23)^2 + 24*d13^4*d23^2+ 16*d12^2*d23^4*cos(theta12)^2 - ...
%     8*d12^2*d23^4*cos(theta12)*cos(theta13)*cos(theta23) - 16*d12^2*d23^4*cos(theta13)^2*cos(theta23)^2 + 24*d12^2*d23^4*cos(theta13)^2 +...
%     8*d12^2*d23^4*cos(theta23)^2 - 24*d12^2*d23^4+ 16*d12^2*d13^4*cos(theta12)^2+ 32*d12^2*d13^4*cos(theta12)*cos(theta13)*cos(theta23)^3 - ...
%     40*d12^2*d13^4*cos(theta12)*cos(theta13)*cos(theta23) - 16*d12^2*d13^4*cos(theta13)^2*cos(theta23)^2 + 8*d12^2*d13^4*cos(theta13)^2 -...
%     32*d12^2*d13^4*cos(theta23)^4 + 56*d12^2*d13^4*cos(theta23)^2 - 24*d12^2*d13^4  - 32*d12^2*d13^2*d23^2*cos(theta12)^2 + ...
%     32*d12^2*d13^2*d23^2*cos(theta12)*cos(theta13)*cos(theta23)^3 + 48*d12^2*d13^2*d23^2*cos(theta12)*cos(theta13)*cos(theta23) - ...
%     32*d12^2*d13^2*d23^2*cos(theta13)^2- 64*d12^2*d13^2*d23^2*cos(theta23)^2 + 48*d12^2*d13^2*d23^2 - 16*d12^4*d13^2*cos(theta12)^2*cos(theta23)^2 +...
%     8*d12^4*d13^2*cos(theta12)^2 + 32*d12^4*d13^2*cos(theta12)*cos(theta13)*cos(theta23)^3 - 40*d12^4*d13^2*cos(theta12)*cos(theta13)*cos(theta23) + ...
%     16*d12^4*d13^2*cos(theta13)^2 - 32*d12^4*d13^2*cos(theta23)^4 + 56*d12^4*d13^2*cos(theta23)^2 - 24*d12^4*d13^2 - 8*d12^4*d23^2*cos(theta12)^2 - ...
%     8*d12^4*d23^2*cos(theta12)*cos(theta13)*cos(theta23) + 32*d12^4*d23^2*cos(theta13)^2*cos(theta23)^2 - 24*d12^4*d23^2*cos(theta13)^2 - ...
%     16*d12^4*d23^2*cos(theta23)^2+ 24*d12^4*d23^2+ 8*d12^6*cos(theta12)*cos(theta13)*cos(theta23) - 16*d12^6*cos(theta13)^2*cos(theta23)^2 +...
%     8*d12^6*cos(theta13)^2 + 8*d12^6*cos(theta23)^2 - 8*d12^6 ;
% 
% 
% aa0=d12^8 - 8*d12^6*d13^2*cos(theta23)^2 + 4*d12^6*d13^2 - 4*d12^6*d23^2 + 16*d12^4*d13^4*cos(theta23)^4 - 16*d12^4*d13^4*cos(theta23)^2 +...
%     6*d12^4*d13^4 + 16*d12^4*d13^2*d23^2*cos(theta23)^2 - 12*d12^4*d13^2*d23^2  - 8*d12^2*d13^6*cos(theta23)^2+ 4*d12^2*d13^6 + ...
%     16*d12^2*d13^4*d23^2*cos(theta23)^2 - 12*d12^2*d13^4*d23^2 - 4*d12^2*d23^6 + d13^8 - 4*d13^6*d23^2 - 4*d13^2*d23^6+ d23^8 + ...
%     6*d13^4*d23^4- 8*d12^2*d13^2*d23^4*cos(theta23)^2 + 12*d12^2*d13^2*d23^4+ 6*d12^4*d23^4 ;
% aa81= 16*theta12^4 - 64*theta12^3*theta13*theta23 + 64*theta12^2*theta13^2*theta23^2 +...
%     32*theta12^2*theta13^2+ 32*theta12^2*theta23^2 - 32*theta12^2 - 64*theta12*theta13^3*theta23- ...
%     64*theta12*theta13*theta23^3 + 64*theta12*theta13*theta23 + 16*theta13^4 + 32*theta13^2*theta23^2 -...
%     32*theta13^2 + 16*theta23^4 - 32*theta23^2 + 16
% 
% aa61= - 64*d23^2*theta12^4*theta13^2 + 32*d23^2*theta12^4 + 128*d23^2*theta12^3*theta13^3*theta23 - ...
%     32*d23^2*theta12^3*theta13*theta23 - 64*d23^2*theta12^2*theta13^4 - 128*d23^2*theta12^2*theta13^2*theta23^2 + ...
%     128*d23^2*theta12^2*theta13^2 + 32*d23^2*theta12^2*theta23^2 - 64*d23^2*theta12^2 - 32*d23^2*theta12*theta13^3*theta23 +...
%     32*d23^2*theta12*theta13*theta23^3 + 32*d23^2*theta12*theta13*theta23 + 32*d23^2*theta13^4+ ...
%     32*d23^2*theta13^2*theta23^2 - 64*d23^2*theta13^2 - 32*d23^2*theta23^2 + 32*d23^2 - 32*d13^2*theta12^4 +...
%     96*d13^2*theta12^3*theta13*theta23 - 64*d13^2*theta12^2*theta13^2*theta23^2 - 32*d13^2*theta12^2*theta13^2 - ...
%     64*d13^2*theta12^2*theta23^2 + 64*d13^2*theta12^2 + 32*d13^2*theta12*theta13^3*theta23 +...
%     96*d13^2*theta12*theta13*theta23^3 - 96*d13^2*theta12*theta13*theta23 - 32*d13^2*theta13^2*theta23^2 + ...
%     32*d13^2*theta13^2 - 32*d13^2*theta23^4 + 64*d13^2*theta23^2 - 32*d13^2  + 32*d12^2*theta12^3*theta13*theta23 - ...
%     64*d12^2*theta12^2*theta13^2*theta23^2 - 32*d12^2*theta12^2*theta13^2 - 32*d12^2*theta12^2*theta23^2+...
%     32*d12^2*theta12^2 + 96*d12^2*theta12*theta13^3*theta23 + 96*d12^2*theta12*theta13*theta23^3 - ...
%     96*d12^2*theta12*theta13*theta23 - 32*d12^2*theta13^4 - 64*d12^2*theta13^2*theta23^2+ 64*d12^2*theta13^2 -...
%     32*d12^2*theta23^4 + 64*d12^2*theta23^2 - 32*d12^2 
% 
% aa41= 16*d23^4*theta12^4 - 32*d23^4*theta12^3*theta13*theta23 + 48*d23^4*theta12^2*theta13^2 + ...
%     16*d23^4*theta12^2*theta23^2 - 40*d23^4*theta12^2 - 32*d23^4*theta12*theta13^3*theta23 + ...
%     16*d23^4*theta12*theta13*theta23 + 16*d23^4*theta13^4 + 16*d23^4*theta13^2*theta23^2 - ...
%     40*d23^4*theta13^2 - 8*d23^4*theta23^2 + 24*d23^4 - 32*d13^2*d23^2*theta12^4 + ...
%     64*d13^2*d23^2*theta12^3*theta13*theta23+ 64*d13^2*d23^2*theta12^2*theta13^2*theta23^2 - ...
%     64*d13^2*d23^2*theta12^2*theta13^2 - 64*d13^2*d23^2*theta12^2*theta23^2 + 80*d13^2*d23^2*theta12^2 -...
%     32*d13^2*d23^2*theta12*theta13^3*theta23 - 32*d13^2*d23^2*theta12*theta13*theta23^3 - ...
%     32*d13^2*d23^2*theta12*theta13*theta23 + 48*d13^2*d23^2*theta13^2 + 48*d13^2*d23^2*theta23^2 - ...
%     48*d13^2*d23^2 + 16*d13^4*theta12^4- 32*d13^4*theta12^3*theta13*theta23 + 16*d13^4*theta12^2*theta13^2 + ...
%     48*d13^4*theta12^2*theta23^2 - 40*d13^4*theta12^2 - 32*d13^4*theta12*theta13*theta23^3 + ...
%     16*d13^4*theta12*theta13*theta23 + 16*d13^4*theta13^2*theta23^2 - 8*d13^4*theta13^2+ 16*d13^4*theta23^4 - ...
%     40*d13^4*theta23^2 + 24*d13^4  - 32*d12^2*d23^2*theta12^3*theta13*theta23 + 64*d12^2*d23^2*theta12^2*theta13^2*theta23^2 -...
%     64*d12^2*d23^2*theta12^2*theta13^2 + 48*d12^2*d23^2*theta12^2+ 64*d12^2*d23^2*theta12*theta13^3*theta23 - ...
%     32*d12^2*d23^2*theta12*theta13*theta23^3 - 32*d12^2*d23^2*theta12*theta13*theta23 - 32*d12^2*d23^2*theta13^4 -...
%     64*d12^2*d23^2*theta13^2*theta23^2+ 80*d12^2*d23^2*theta13^2 + 48*d12^2*d23^2*theta23^2 -...
%     48*d12^2*d23^2 - 32*d12^2*d13^2*theta12^3*theta13*theta23 + 64*d12^2*d13^2*theta12^2*theta13^2*theta23^2 + ...
%     32*d12^2*d13^2*theta12^2*theta23^2 - 48*d12^2*d13^2*theta12^2 - 32*d12^2*d13^2*theta12*theta13^3*theta23 - ...
%     128*d12^2*d13^2*theta12*theta13*theta23^3+ 160*d12^2*d13^2*theta12*theta13*theta23 + ...
%     32*d12^2*d13^2*theta13^2*theta23^2 - 48*d12^2*d13^2*theta13^2 + 64*d12^2*d13^2*theta23^4 - 112*d12^2*d13^2*theta23^2 +...
%     48*d12^2*d13^2 + 16*d12^4*theta12^2*theta13^2 + 16*d12^4*theta12^2*theta23^2 - 8*d12^4*theta12^2 - ...
%     32*d12^4*theta12*theta13^3*theta23 - 32*d12^4*theta12*theta13*theta23^3 + 16*d12^4*theta12*theta13*theta23 +...
%     16*d12^4*theta13^4 + 48*d12^4*theta13^2*theta23^2 - 40*d12^4*theta13^2 + 16*d12^4*theta23^4 - 40*d12^4*theta23^2 + 24*d12^4;
% 
% aa21=- 8*d23^6*theta12^2 + 8*d23^6*theta12*theta13*theta23 - 8*d23^6*theta13^2 + 8*d23^6  - ...
%     16*d13^2*d23^4*theta12^2*theta23^2 + 24*d13^2*d23^4*theta12^2 - 8*d13^2*d23^4*theta12*theta13*theta23 +...
%     16*d13^2*d23^4*theta13^2+ 8*d13^2*d23^4*theta23^2 - 24*d13^2*d23^4 - 16*d13^6*theta12^2*theta23^2 + ...
%     8*d13^6*theta12^2 + 8*d13^6*theta12*theta13*theta23 + 8*d13^6*theta23^2 - 8*d13^6  + ...
%     32*d13^4*d23^2*theta12^2*theta23^2 - 24*d13^4*d23^2*theta12^2 - 8*d13^4*d23^2*theta12*theta13*theta23 - ...
%     8*d13^4*d23^2*theta13^2 - 16*d13^4*d23^2*theta23^2 + 24*d13^4*d23^2+ 16*d12^2*d23^4*theta12^2 - ...
%     8*d12^2*d23^4*theta12*theta13*theta23 - 16*d12^2*d23^4*theta13^2*theta23^2 + 24*d12^2*d23^4*theta13^2 +...
%     8*d12^2*d23^4*theta23^2 - 24*d12^2*d23^4+ 16*d12^2*d13^4*theta12^2+ 32*d12^2*d13^4*theta12*theta13*theta23^3 - ...
%     40*d12^2*d13^4*theta12*theta13*theta23 - 16*d12^2*d13^4*theta13^2*theta23^2 + 8*d12^2*d13^4*theta13^2 -...
%     32*d12^2*d13^4*theta23^4 + 56*d12^2*d13^4*theta23^2 - 24*d12^2*d13^4  - 32*d12^2*d13^2*d23^2*theta12^2 + ...
%     32*d12^2*d13^2*d23^2*theta12*theta13*theta23^3 + 48*d12^2*d13^2*d23^2*theta12*theta13*theta23 - ...
%     32*d12^2*d13^2*d23^2*theta13^2- 64*d12^2*d13^2*d23^2*theta23^2 + 48*d12^2*d13^2*d23^2 - 16*d12^4*d13^2*theta12^2*theta23^2 +...
%     8*d12^4*d13^2*theta12^2 + 32*d12^4*d13^2*theta12*theta13*theta23^3 - 40*d12^4*d13^2*theta12*theta13*theta23 + ...
%     16*d12^4*d13^2*theta13^2 - 32*d12^4*d13^2*theta23^4 + 56*d12^4*d13^2*theta23^2 - 24*d12^4*d13^2 - 8*d12^4*d23^2*theta12^2 - ...
%     8*d12^4*d23^2*theta12*theta13*theta23 + 32*d12^4*d23^2*theta13^2*theta23^2 - 24*d12^4*d23^2*theta13^2 - ...
%     16*d12^4*d23^2*theta23^2+ 24*d12^4*d23^2+ 8*d12^6*theta12*theta13*theta23 - 16*d12^6*theta13^2*theta23^2 +...
%     8*d12^6*theta13^2 + 8*d12^6*theta23^2 - 8*d12^6 ;
% 
% 
% aa01=d12^8 - 8*d12^6*d13^2*theta23^2 + 4*d12^6*d13^2 - 4*d12^6*d23^2 + 16*d12^4*d13^4*theta23^4 - 16*d12^4*d13^4*theta23^2 +...
%     6*d12^4*d13^4 + 16*d12^4*d13^2*d23^2*theta23^2 - 12*d12^4*d13^2*d23^2  - 8*d12^2*d13^6*theta23^2+ 4*d12^2*d13^6 + ...
%     16*d12^2*d13^4*d23^2*theta23^2 - 12*d12^2*d13^4*d23^2 - 4*d12^2*d23^6 + d13^8 - 4*d13^6*d23^2 - 4*d13^2*d23^6+ d23^8 + ...
%     6*d13^4*d23^4- 8*d12^2*d13^2*d23^4*theta23^2 + 12*d12^2*d13^2*d23^4+ 6*d12^4*d23^4 ;


aa8= (16*theta12^4 - 64*theta12^3*theta13*theta23 + 64*theta12^2*theta13^2*theta23^2 + 32*theta12^2*theta13^2 + 32*theta12^2*theta23^2- 32*theta12^2 -...
    64*theta12*theta13^3*theta23- 64*theta12*theta13*theta23^3 + 64*theta12*theta13*theta23 + 16*theta13^4 + 32*theta13^2*theta23^2 - 32*theta13^2 +...
    16*theta23^4 - 32*theta23^2 + 16);

aa6=(- 64*d23^2*theta12^4*theta13^2+ 32*d23^2*theta12^4+ 128*d23^2*theta12^3*theta13^3*theta23- 32*d23^2*theta12^3*theta13*theta23-...
    64*d23^2*theta12^2*theta13^4- 128*d23^2*theta12^2*theta13^2*theta23^2+ 128*d23^2*theta12^2*theta13^2+ 32*d23^2*theta12^2*theta23^2-...
    64*d23^2*theta12^2- 32*d23^2*theta12*theta13^3*theta23+ 32*d23^2*theta12*theta13*theta23^3+ 32*d23^2*theta12*theta13*theta23+ ...
    32*d23^2*theta13^4+ 32*d23^2*theta13^2*theta23^2- 64*d23^2*theta13^2- 32*d23^2*theta23^2+ 32*d23^2 - 32*d13^2*theta12^4+...
    96*d13^2*theta12^3*theta13*theta23- 64*d13^2*theta12^2*theta13^2*theta23^2- 32*d13^2*theta12^2*theta13^2- 64*d13^2*theta12^2*theta23^2+...
    64*d13^2*theta12^2+ 32*d13^2*theta12*theta13^3*theta23+ 96*d13^2*theta12*theta13*theta23^3- 96*d13^2*theta12*theta13*theta23-...
    32*d13^2*theta13^2*theta23^2+ 32*d13^2*theta13^2- 32*d13^2*theta23^4+ 64*d13^2*theta23^2- 32*d13^2+ 32*d12^2*theta12^3*theta13*theta23- ...
    64*d12^2*theta12^2*theta13^2*theta23^2- 32*d12^2*theta12^2*theta13^2- 32*d12^2*theta12^2*theta23^2+ ...
    32*d12^2*theta12^2+ 96*d12^2*theta12*theta13^3*theta23+ 96*d12^2*theta12*theta13*theta23^3- 96*d12^2*theta12*theta13*theta23- 32*d12^2*theta13^4- ...
    64*d12^2*theta13^2*theta23^2+ 64*d12^2*theta13^2- 32*d12^2*theta23^4+ 64*d12^2*theta23^2- 32*d12^2);

%

aa4=(16*d23^4*theta12^4 - 32*d23^4*theta12^3*theta13*theta23 + 48*d23^4*theta12^2*theta13^2+...
    16*d23^4*theta12^2*theta23^2 - 40*d23^4*theta12^2 - 32*d23^4*theta12*theta13^3*theta23 + ...
    16*d23^4*theta12*theta13*theta23 +16*d23^4*theta13^4 + 16*d23^4*theta13^2*theta23^2 - ...
    40*d23^4*theta13^2 - 8*d23^4*theta23^2 + 24*d23^4 - 32*d13^2*d23^2*theta12^4 + ...
    64*d13^2*d23^2*theta12^3*theta13*theta23+ 64*d13^2*d23^2*theta12^2*theta13^2*theta23^2 -...
    64*d13^2*d23^2*theta12^2*theta13^2 - 64*d13^2*d23^2*theta12^2*theta23^2 + 80*d13^2*d23^2*theta12^2 - ...
    32*d13^2*d23^2*theta12*theta13^3*theta23 - 32*d13^2*d23^2*theta12*theta13*theta23^3 - ...
    32*d13^2*d23^2*theta12*theta13*theta23+ 48*d13^2*d23^2*theta13^2 + 48*d13^2*d23^2*theta23^2 - ...
    48*d13^2*d23^2 + 16*d13^4*theta12^4 - 32*d13^4*theta12^3*theta13*theta23 + 16*d13^4*theta12^2*theta13^2 + ...
    48*d13^4*theta12^2*theta23^2 - 40*d13^4*theta12^2 - 32*d13^4*theta12*theta13*theta23^3 +...
    16*d13^4*theta12*theta13*theta23 + 16*d13^4*theta13^2*theta23^2 - 8*d13^4*theta13^2 + 16*d13^4*theta23^4- ...
    40*d13^4*theta23^2 + 24*d13^4 - 32*d12^2*d23^2*theta12^3*theta13*theta23 + 64*d12^2*d23^2*theta12^2*theta13^2*theta23^2- ...
    64*d12^2*d23^2*theta12^2*theta13^2 + 48*d12^2*d23^2*theta12^2 + 64*d12^2*d23^2*theta12*theta13^3*theta23 -...
    32*d12^2*d23^2*theta12*theta13*theta23^3 - 32*d12^2*d23^2*theta12*theta13*theta23 - 32*d12^2*d23^2*theta13^4 -...
    64*d12^2*d23^2*theta13^2*theta23^2 + 80*d12^2*d23^2*theta13^2 + 48*d12^2*d23^2*theta23^2 - ...
    48*d12^2*d23^2- 32*d12^2*d13^2*theta12^3*theta13*theta23 + 64*d12^2*d13^2*theta12^2*theta13^2*theta23^2 +...
    32*d12^2*d13^2*theta12^2*theta23^2- 48*d12^2*d13^2*theta12^2 - 32*d12^2*d13^2*theta12*theta13^3*theta23 -...
    128*d12^2*d13^2*theta12*theta13*theta23^3 + 160*d12^2*d13^2*theta12*theta13*theta23+ ...
    32*d12^2*d13^2*theta13^2*theta23^2 - 48*d12^2*d13^2*theta13^2 + 64*d12^2*d13^2*theta23^4 - 112*d12^2*d13^2*theta23^2 +...
    48*d12^2*d13^2 +  16*d12^4*theta12^2*theta13^2 + 16*d12^4*theta12^2*theta23^2 - 8*d12^4*theta12^2 - ...
    32*d12^4*theta12*theta13^3*theta23 - 32*d12^4*theta12*theta13*theta23^3 + 16*d12^4*theta12*theta13*theta23 +...
    16*d12^4*theta13^4 + 48*d12^4*theta13^2*theta23^2 - 40*d12^4*theta13^2+ 16*d12^4*theta23^4 - 40*d12^4*theta23^2 + 24*d12^4);


aa2=(- 8*d23^6*theta12^2 + 8*d23^6*theta12*theta13*theta23- 8*d23^6*theta13^2 + 8*d23^6 -...
    16*d13^2*d23^4*theta12^2*theta23^2 + 24*d13^2*d23^4*theta12^2 - 8*d13^2*d23^4*theta12*theta13*theta23+ ...
    16*d13^2*d23^4*theta13^2 + 8*d13^2*d23^4*theta23^2 - 24*d13^2*d23^4 - 16*d13^6*theta12^2*theta23^2+ ...
    8*d13^6*theta12^2 + 8*d13^6*theta12*theta13*theta23 + 8*d13^6*theta23^2 - 8*d13^6  +...
    32*d13^4*d23^2*theta12^2*theta23^2- 24*d13^4*d23^2*theta12^2 - 8*d13^4*d23^2*theta12*theta13*theta23 - ...
    8*d13^4*d23^2*theta13^2 - 16*d13^4*d23^2*theta23^2 + 24*d13^4*d23^2+ 16*d12^2*d23^4*theta12^2 -...
    8*d12^2*d23^4*theta12*theta13*theta23 - 16*d12^2*d23^4*theta13^2*theta23^2 + 24*d12^2*d23^4*theta13^2+...
    8*d12^2*d23^4*theta23^2 - 24*d12^2*d23^4 + 16*d12^2*d13^4*theta12^2 + 32*d12^2*d13^4*theta12*theta13*theta23^3 -...
    40*d12^2*d13^4*theta12*theta13*theta23- 16*d12^2*d13^4*theta13^2*theta23^2 + 8*d12^2*d13^4*theta13^2 - ...
    32*d12^2*d13^4*theta23^4 + 56*d12^2*d13^4*theta23^2 - 24*d12^2*d13^4- 32*d12^2*d13^2*d23^2*theta12^2 + ...
    32*d12^2*d13^2*d23^2*theta12*theta13*theta23^3+ 48*d12^2*d13^2*d23^2*theta12*theta13*theta23 - ...
    32*d12^2*d13^2*d23^2*theta13^2 - 64*d12^2*d13^2*d23^2*theta23^2 + 48*d12^2*d13^2*d23^2- 16*d12^4*d13^2*theta12^2*theta23^2 +...
    8*d12^4*d13^2*theta12^2 + 32*d12^4*d13^2*theta12*theta13*theta23^3 - 40*d12^4*d13^2*theta12*theta13*theta23+...
    16*d12^4*d13^2*theta13^2 - 32*d12^4*d13^2*theta23^4 + 56*d12^4*d13^2*theta23^2 - 24*d12^4*d13^2 - 8*d12^4*d23^2*theta12^2-...
    8*d12^4*d23^2*theta12*theta13*theta23 + 32*d12^4*d23^2*theta13^2*theta23^2 - 24*d12^4*d23^2*theta13^2 -...
    16*d12^4*d23^2*theta23^2 + 24*d12^4*d23^2+ 8*d12^6*theta12*theta13*theta23 - 16*d12^6*theta13^2*theta23^2 +...
    8*d12^6*theta13^2+ 8*d12^6*theta23^2 - 8*d12^6);

aa0=(d12^8 - 8*d12^6*d13^2*theta23^2 + 4*d12^6*d13^2 - 4*d12^6*d23^2  + 16*d12^4*d13^4*theta23^4 - 16*d12^4*d13^4*theta23^2 + 6*d12^4*d13^4 + ...
    16*d12^4*d13^2*d23^2*theta23^2 - 12*d12^4*d13^2*d23^2- 8*d12^2*d13^6*theta23^2 + 4*d12^2*d13^6 + 16*d12^2*d13^4*d23^2*theta23^2 - 12*d12^2*d13^4*d23^2 -...
    4*d12^2*d23^6 + d13^8 - 4*d13^6*d23^2 - 4*d13^2*d23^6  + d23^8  +6*d13^4*d23^4 - 8*d12^2*d13^2*d23^4*theta23^2 + 12*d12^2*d13^2*d23^4 + 6*d12^4*d23^4);
% diff=[aa81,aa61,aa41,aa21,aa01]-[aa8,aa6,aa4,aa2,aa0];