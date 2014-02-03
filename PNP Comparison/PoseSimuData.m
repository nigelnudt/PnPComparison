clc;
clear all;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1), plot( DA(:,1), DA(:,2), '*');   grid;   xlabel('u');   ylabel('v');   title('像点分布');
figure(2), plot3( DA(:,3),DA(:,5),DA(:,4), 'o');   grid;   xlabel('xw');   ylabel('zw');    zlabel('yw');  title('物点分布');
%%检验仿真误差
M1 = RR* ( DA(:,3:5)' - repmat( tt, 1, num ) );
for i=1:num,
    M1(1,i) = M1(1,i) / M1(3,i);
    M1(2,i) = M1(2,i) / M1(3,i);
    M1(3,i) = 0.0;
end
M1 = M1';
%%物点重投影
figure(3), plot( M1(:,1), M1(:,2),'o'); grid;   xlabel('u');   ylabel('v'); title('物点重投影');
figure(4), plot(DA(:,1), DA(:,2), '*', M1(:,1), M1(:,2),'o');  legend('像点', '重投影点');grid;   xlabel('u');   ylabel('v'); title('重投影误差 2D');
B  = DA(:,1:2) - M1(:,1:2);
figure(5), plot( B(:,1), '-o'); grid;   xlabel('Point Index');   ylabel('u'); title('重投影误差 Err_u');
figure(6), plot( B(:,2), '-o'); grid;   xlabel('Point Index');   ylabel('v'); title('重投影误差 Err_v');
%%计算物点到光心的距离
Dt = DA(:,3:5)' - repmat( tt, 1, num );
for i=1:num,
    XX(i) = norm(Dt(:,i));
end
%%%2013-09-03，验证相似变换恢复的最小二乘公式


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%获取蚁群算法的求解结果，分析误差
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
%%输出适应度
Optimal_Fitness = min( Fitness )
figure, plot(Fitness,'b-');  grid;   xlabel('迭代次数 t');   ylabel('适应度 f' ); title('蚁群适应度的变化');
%对比位姿参数的真值与最优解
ss   = [1,2,3,4,5,6];
Pose = [Alfa, Beta, Gama, tt(1), tt(2), tt(3)]';
figure, plot( ss, Pose,'b-o', ss, Optimal_Pose, 'r-*'); grid;   xlabel('参数');   ylabel('参数取值' ); title('真值和最优解对比');

%%检验仿真误差
M2 = Optimal_RR * ( DA(:,3:5)' - repmat( Optimal_tt, 1, num ) );
for i=1:num,
    M2(1,i) = M2(1,i) / M2(3,i);
    M2(2,i) = M2(2,i) / M2(3,i);
    M2(3,i)  = 0.0;
end
M2 = M2';
%%物点重投影
figure, plot(DA(:,1), DA(:,2), '*', M2(:,1), M2(:,2),'o');  legend('像点', '重投影点');grid;   xlabel('u');   ylabel('v'); title('蚁群优化结果的重投影误差分布');