clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cFileName = 'D:\������Ӿ�����\��־����ʦ����\N��λ������ 2013-09-13\Pose_PointPattern1.dat';
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1), plot( DA(:,1), DA(:,2), '*');   grid;   xlabel('u');   ylabel('v');   title('���ֲ�');
figure(2), plot3( DA(:,3),DA(:,5),DA(:,4), 'o');   grid;   xlabel('xw');   ylabel('zw');    zlabel('yw');  title('���ֲ�');
%%����������
M1 = RR* ( DA(:,3:5)' - repmat( tt, 1, num ) );
for i=1:num,
    M1(1,i) = M1(1,i) / M1(3,i);
    M1(2,i) = M1(2,i) / M1(3,i);
    M1(3,i) = 0.0;
end
M1 = M1';
%%�����ͶӰ
figure(3), plot( M1(:,1), M1(:,2),'o'); grid;   xlabel('u');   ylabel('v'); title('�����ͶӰ');
figure(4), plot(DA(:,1), DA(:,2), '*', M1(:,1), M1(:,2),'o');  legend('���', '��ͶӰ��');grid;   xlabel('u');   ylabel('v'); title('��ͶӰ��� 2D');
B  = DA(:,1:2) - M1(:,1:2);
figure(5), plot( B(:,1), '-o'); grid;   xlabel('Point Index');   ylabel('u'); title('��ͶӰ��� Err_u');
figure(6), plot( B(:,2), '-o'); grid;   xlabel('Point Index');   ylabel('v'); title('��ͶӰ��� Err_v');
%%������㵽���ĵľ���
Dt = DA(:,3:5)' - repmat( tt, 1, num );
for i=1:num,
    XX(i) = norm(Dt(:,i));
end
%%%2013-09-03����֤���Ʊ任�ָ�����С���˹�ʽ


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%��ȡ��Ⱥ�㷨����������������
cFileName = 'd:\POSE_Ants.dat';
f2                    = fopen( cFileName,'r' );
[Optimal_Pose,nCount] = fread( f2,  6,     'double' );
[Optimal_RR, nCount]  = fread( f2,  [3,3], 'double');
[Optimal_tt, nCount]  = fread( f2,  3,     'double');
[nMaxTimes,  nCount]  = fread( f2,  1,     'int');
[Fitness,    nCount]  = fread( f2,  nMaxTimes,'double'); 
fclose(f2);
Optimal_RR = Optimal_RR';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%�����Ӧ��
Optimal_Fitness = min( Fitness )
figure, plot(Fitness,'b-');  grid;   xlabel('�������� t');   ylabel('��Ӧ�� f' ); title('��Ⱥ��Ӧ�ȵı仯');
%�Ա�λ�˲�������ֵ�����Ž�
ss   = [1,2,3,4,5,6];
Pose = [Alfa, Beta, Gama, tt(1), tt(2), tt(3)]';
figure, plot( ss, Pose,'b-o', ss, Optimal_Pose, 'r-*'); grid;   xlabel('����');   ylabel('����ȡֵ' ); title('��ֵ�����Ž�Ա�');

%%����������
M2 = Optimal_RR * ( DA(:,3:5)' - repmat( Optimal_tt, 1, num ) );
for i=1:num,
    M2(1,i) = M2(1,i) / M2(3,i);
    M2(2,i) = M2(2,i) / M2(3,i);
    M2(3,i)  = 0.0;
end
M2 = M2';
%%�����ͶӰ
figure, plot(DA(:,1), DA(:,2), '*', M2(:,1), M2(:,2),'o');  legend('���', '��ͶӰ��');grid;   xlabel('u');   ylabel('v'); title('��Ⱥ�Ż��������ͶӰ���ֲ�');