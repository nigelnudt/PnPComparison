clc;
clear all;
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        < N 点最小二乘算法 >
% 思  想：已知 N 个点的三维“测量坐标”和 “世界坐标”，求解相似变换 R 和 T
% 设计者：张志龙
% 时  间：2013-09-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        第 1 步：读取真值数据
%--------------------------------------------------------------------------
cFileName = 'D:\计算机视觉测量\张志龙老师合作\N点位姿数据 2013-09-13\Pose_PointPattern2.dat';
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        第 2 步：构造3D点云   (经验证正确)
%--------------------------------------------------------------------------
DataIn(:,4:6) = DA(:,3:5);
DataIn(:,1:3) = ( RR * ( DA(:,3:5)' - repmat( tt, 1, num ) ) )';
figure, plot3( DataIn(:,1), DataIn(:,2), DataIn(:,3), 'o', DataIn(:,4), DataIn(:,5), DataIn(:,6), '*');    grid;    legend('测量系', '参考系');   
axis([-0.5 3 -1 2.5 0 3.5])
xlabel('xw');    ylabel('zw');    zlabel('yw');    title('物点分布');
for index1=1:9
    for index2=index1+1:10
        distance(index1,index2)=norm(DataIn(index1,1:3)-DataIn(index2,1:3))-norm(DataIn(index1,4:6)-DataIn(index2,4:6));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        第 3 步：构造单位方向矢量
%        经验证，为单位方向矢量
%--------------------------------------------------------------------------
Temp = DataIn - repmat( DataIn(1,:), 10, 1 );
DA3  = Temp( 2:10, : );
[ lin, col ] = size( DA3 );
for i=1:lin,
    ss1 = norm( DA3(i, 1:3) );    DA3(i,1:3) = DA3(i, 1:3) / ss1;
    ss2 = norm( DA3(i, 4:6) );    DA3(i,4:6) = DA3(i, 4:6) / ss2;   
end  
Index = lin + 1;
for i = 1:lin-1,
    for j = i+1:lin,
        aa = cross( DA3(i,1:3), DA3(j,1:3) );    aa = aa / norm( aa );
        bb = cross( DA3(i,4:6), DA3(j,4:6) );    bb = bb / norm( bb );
        DA3(Index,:) = [aa bb];
        Index = Index + 1;
    end
end
%%
%%新的旋转由小杨来写，已经证明可行
[lin, col] = size( DA3 ); 
B_yang =  zeros( 3*lin, 1 );
A_yang = zeros( 3*lin, 9 );
for i=1:lin,
    A_yang(3*i-2,:) =[DA3(i,1),DA3(i,2),DA3(i,3),0,0,0,0,0,0];
    A_yang(3*i-1,:) =[0,0,0,DA3(i,1),DA3(i,2),DA3(i,3),0,0,0];
    A_yang(3*i,:) =[0,0,0,0,0,0,DA3(i,1),DA3(i,2),DA3(i,3)];
    B_yang(3*i-2:3*i,1)=DA3(i,4:6)';
end
x_yang = inv(A_yang'*A_yang)*A_yang'*B_yang;
y_yang = sqrt(3) * x_yang / norm( x_yang );
Rotate_yang = [ y_yang(1) y_yang(2) y_yang(3) 
           y_yang(4) y_yang(5) y_yang(6);
           y_yang(7) y_yang(8) y_yang(9)];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        第 4 步：求解旋转矩阵
%--------------------------------------------------------------------------
[lin, col] = size( DA3 ); 
B =  ones( lin, 1 );
A = zeros( lin, 9 );
for i=1:lin,
    A(i,:) = [ DA3(i,4)* DA3(i,1),  DA3(i,5)* DA3(i,1),  DA3(i,6)* DA3(i,1), DA3(i,4)* DA3(i,2),  DA3(i,5)* DA3(i,2),  DA3(i,6)* DA3(i,2), DA3(i,4)* DA3(i,3),  DA3(i,5)* DA3(i,3),  DA3(i,6)* DA3(i,3) ];
end
x = inv(A'*A)*A'*B;
y = sqrt(3) * x / norm( x );
Rotate = [ y(1) y(2) y(3) 
           y(4) y(5) y(6);
           y(7) y(8) y(9)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        第 5 步：求解平移矢量
%--------------------------------------------------------------------------
DA4( :,1:3 ) = DataIn( :,1:3 );
DA4( :,4:6 ) = ( Rotate * DataIn( :,4:6 )' )'; 