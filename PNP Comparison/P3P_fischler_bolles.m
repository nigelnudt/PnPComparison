function [R_fischler,T_fishler]=P3P_fischler_bolles(x3d_h,x2d_h_normlized)
R_fischler=eye(3);
T_fishler=zeros(3,1);
% %x3d_h(1:50,:),x2d_h_normlized(1:50,:)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %˵����
% %    ��δ���ʵ����Fischler�㷨�����㷨�Ļ���˼���ǣ�
% %    (1)��������̺��ľ���Լ��������̺��ĽǶ�Լ�������㡰���-���ġ�����;
% %    (2)���ݡ����-���ġ��������������������ϵ�е�����;
% %    (3)����������������ϵ����������ϵ�е������Ӧ��ȷ����ת���� R �� ƽ��ʸ�� T;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;
% clear all;
% %--------------------------------------------------------------------------
% %%%��һ������ȡ������� ���������ɣ���ֵ������������
% %--------------------------------------------------------------------------
% cFileName = 'D:\������Ӿ�����\��־����ʦ����\N��λ������ 2013-09-13\Pose_PointPattern1.dat';
% fid = fopen(cFileName,'r');
% [Alfa, nCount] = fread(fid,1,'double');        %��ת�Ƕ�
% [Beta, nCount] = fread(fid,1,'double');        %��ת�Ƕ�
% [Gama, nCount] = fread(fid,1,'double');        %��ת�Ƕ�
% [RR,   nCount] = fread(fid,[3,3],'double');    %��ת�任����
% [tt,   nCount] = fread(fid, 3,'double');       %ƽ��ʸ��
% [num,  nCount] = fread(fid, 1,'int');
% [DA,   nCount] = fread(fid,[5,num],'double'); 
% fclose(fid);
% RR = RR';
% DA = DA';
%--------------------------------------------------------------------------
%%%�ڶ���, ѡȡ3�Ե㣬ʵ�� Fischler ���������㷨������3����㵽���ĵľ���
%  ע��: �����е�[a,b,c]�ֱ��Ӧ�㷨�е�[1,2,3]
%--------------------------------------------------------------------------
DataIn(:,1:2) = x2d_h_normlized(1:3,1:2);
DataIn(:,3)   = 1;
DataIn(:,4:6) = x3d_h(1:3,1:3);
%%(1)��������̺��ľ���Լ��R(i,j) �� ����̺��ĽǶ�Լ��A(i,j)
for i=1:3,
    for j=1:3,
        R(i,j) = norm( DataIn(i,4:6) - DataIn(j,4:6) );
        A(i,j) = DataIn(i,1:3) * DataIn(j,1:3)' / norm( DataIn(i,1:3) ) / norm( DataIn(j,1:3) );        
    end
end 
K1 = (R(2,3)^2)/(R(1,3)^2);
K2 = (R(2,3)^2)/(R(1,2)^2);
%%(2)ȷ��5������ʽϵ��
G4 = (K1*K2-K1-K2)^2 - 4*K1*K2*(A(2,3)^2);
G3 = 4*(K1*K2-K1-K2)*K2*(1-K1)*A(1,2) + 4*K1*A(2,3)*( (K1*K2+K2-K1)*A(1,3) + 2*K2*A(1,2)*A(2,3) );
G2 = (2*K2*(1-K1)*A(1,2))^2 + 2*(K1*K2+K1-K2)*(K1*K2-K1-K2)+4*K1*( (K1-K2)*(A(2,3)^2) + (1-K2)*K1*(A(1,3)^2) - 2*K2*(1+K1)*A(1,2)*A(1,3)*A(2,3) );
G1 = 4*(K1*K2+K1-K2)*K2*(1-K1)*A(1,2) + 4*K1*( (K1*K2-K1+K2)*A(1,3)*A(2,3) + 2*K1*K2*A(1,2)* (A(1,3)^2) );
G0 = (K1*K2+K1-K2)^2 - 4*K1*K1*K2*(A(1,3)^2);
%%(3)������ʽ�ĸ�
cof = [G4,G3,G2,G1,G0];
res =real(roots( cof ));
%%(4)���ݴ���0��ʵ��������Ч�⣬�������Ӧ�Ķ��������
Index = 1;
for i=1:4,
    if imag( res(i) ) == 0 & real(res(i)) > 0,   %%����Ǵ���0��ʵ��
        x = res(i);
        %%%%���� x ���� a
        a = R(1,2) / sqrt( x^2 -2*x*A(1,2) + 1 );
        %%%%���� x ���� b
        b = x * a;
        %%%%���� a ���� y
        tmp = sqrt( A(1,3)^2 + (R(1,3)/a)^2 - 1 );
        y1 = A(1,3) + tmp;
        c1 = y1 * a;
        t1 = abs( b^2 + c1^2 - 2*b*c1*A(2,3) - R(2,3)^2 );
        y2 = A(1,3) - tmp;            
        c2 = y2 * a;
        t2 = abs( b^2 + c2^2 - 2*b*c2*A(2,3) - R(2,3)^2 );
        %%%1Ԫ2�η��̻����2���� y1 �� y2����������С��һ����
        if t2 < t1,
            c = c2;
        else
            c = c1;
        end
        %%%���������(a,b,c)��Ӧ�Ķ��������
        Error = abs( a^2 + b^2 - 2*a*b*A(1,2) - R(1,2)^2 ) + abs( a^2 + c^2 - 2*a*c*A(1,3) - R(1,3)^2 ) + abs( b^2 + c^2 - 2*b*c*A(2,3) - R(2,3)^2 );  
        %%%������
        Result(Index,1:4) = [a,b,c,Error];     
        Index = Index + 1;
    end %end if
end
% [min_error,right_index]=min(Result(:,4))
verify_error=100;%set an threshold for the test of which solution is correct by the fourth point 
for solution_index=1:Index-1
    x3d_camera_p1=Result(solution_index,1)*x2d_h_normlized(1,:)./norm(x2d_h_normlized(1,:));
    x3d_camera_p2=Result(solution_index,2)*x2d_h_normlized(2,:)./norm(x2d_h_normlized(2,:));
    x3d_camera_p3=Result(solution_index,3)*x2d_h_normlized(3,:)./norm(x2d_h_normlized(3,:));
    % x3d_camera_p4=solution_confirmed_x4*x2d_h_normlized(4,:)./norm(x2d_h_normlized(4,:));
    x3d_camera=[x3d_camera_p1;x3d_camera_p2;x3d_camera_p3];
    Xc3=x3d_camera;
    sum_world=0;
    sum_camera_world=0;
    world_coor_mean=mean(x3d_h(1:3,1:3),1);
    world_coor_matrix=x3d_h(1:3,1:3)-repmat(world_coor_mean,3,1);
    camera_coor_mean=mean(x3d_camera,1);
    camera_coor_matrix=x3d_camera-repmat(camera_coor_mean,3,1);
    for point_index=1:3
        sum_camera_world=sum_camera_world+norm(camera_coor_matrix(point_index,:))*norm(camera_coor_matrix(point_index,:));
        sum_world=sum_world+norm(norm(world_coor_matrix(point_index,:)))*norm(world_coor_matrix(point_index,:));
    end
    scale_camera_world=sum_camera_world/sum_world;
    [UR,DR,VR]=svd(world_coor_matrix'*camera_coor_matrix);
    % [UR1,DR1,VR1]=svd(DA(:,3:5)'*b_matrix1');
    AA=det(UR)*det(VR);
    if AA<0
        Rp3=VR*[1,0,0;0,1,0;0,0,-1]*UR';
    else
        Rp3=VR*UR';
    end
    tp3=1/scale_camera_world*Rp3'*camera_coor_mean'-world_coor_mean';
    p4_camera=Rp3*((x3d_h(4,1:3))'+tp3);
    p4_error=norm([p4_camera(1)/p4_camera(3),p4_camera(2)/p4_camera(3)]-x2d_h_normlized(4,1:2));
    if p4_error<verify_error
        R_fischler=Rp3;
        T_fishler=Rp3*tp3;
        verify_error=p4_error;
    end
        
    % diff=rotation_matrix-R(1:3,1:3);
    % rotation_matrix1=VR1*UR1';
    % rotation_matrix=rotation_matrix*[1,0,0;0,-1,0;0,0,1];
    
end

% %-----------------------------------------------------------
% %%%������, ���ݡ����-���ġ�������� ��ת���� R �� ƽ��ʸ�� T
% %  ע��: �����е�[a,b,c]�ֱ��Ӧ�㷨�е�[1,2,3]
% %-----------------------------------------------------------
% %%%�������Ʊ任��ȷ���任����
% a = Result(2,1);
% b = Result(2,2);
% c = Result(2,3);
% DataOut(:,4:6) = DataIn(:,4:6);
% % DataOut(1,1:3) = (a * DataIn(1,1:3)')';
% % DataOut(2,1:3) = (b * DataIn(2,1:3)')';
% % DataOut(3,1:3) = (c * DataIn(3,1:3)')';
% DataOut(1,1:3) = (a * DataIn(1,1:3)')'/norm(DataIn(1,1:3));
% DataOut(2,1:3) = (b * DataIn(2,1:3)')'/norm(DataIn(2,1:3));
% DataOut(3,1:3) = (c * DataIn(3,1:3)')'/norm(DataIn(3,1:3));
% MA = zeros(9,12);
% MA([1,4,7],1:3) = DataOut(:,1:3); MA([1,4,7],10) = 1;
% MA([2,5,8],4:6) = DataOut(:,1:3); MA([2,5,8],11) = 1;
% MA([3,6,9],7:9) = DataOut(:,1:3); MA([3,6,9],12) = 1;
% VB = [ DataOut(1,4:6), DataOut(2,4:6), DataOut(3,4:6) ];
% VX = inv( MA'* MA ) * (MA' * VB');