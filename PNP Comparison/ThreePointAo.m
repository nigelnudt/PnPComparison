clc;
clear all;
close all
%--------------------------------------------------------------------------
%这段代码寻找旋转矩阵 Rotate 和平移矢量 Translate, 将世界坐标系转换到相机坐标系
%设 计 者：张志龙
%设计时间：2013-09-10
%--------------------------------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3D点对数据
Pnt = [1.2914         0    1.2990    0.0021    0.0865    3.5139;
       1.1951         0    1.1596    0.0199   -0.0600    3.4306;
       1.1927         0    1.1711    0.0106   -0.0536    3.4338];   
figure,  legend( '世界坐标系', '相机坐标系' );
subplot(2,2,1,'V6'),  plot3( Pnt(:,1), Pnt(:,2), Pnt(:,3), 'o', Pnt(:,4), Pnt(:,5), Pnt(:,6), '*');   grid;
xlabel('x');   ylabel('y');    zlabel('z');   title('原始数据');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%第一步：计算 3 点确定平面的法矢量 
%--------------------------------------------------------------------------   
Vector(:,1:6) = Pnt(:,1:6) - repmat( Pnt(1,1:6), 3, 1 );
%根据物点的世界坐标值，拟合法矢量 n
n  = cross( Vector(2,1:3), Vector(3,1:3) );
n  = n / norm( n );
%根据物点的相机坐标值拟合法矢量 n1
n1 = cross( Vector(2,4:6), Vector(3,4:6) );
n1 = n1 / norm( n1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%第二步：根据法矢量确定旋转矩阵,旋转两个点集，使与oxy平面平行
%--------------------------------------------------------------------------

%--------------------------------------------
% 1.0 矢量 n
tmp = norm( n(1:2) );
if tmp == 0,
   tmp = 0.00000000000001;
end
%%%%角度Alfa
Theta = asin( abs( n(2) ) / tmp );
Alfa = 0;
if n(1) >= 0
    if n(2) >= 0
        Alfa = Theta;
    else
        Alfa = 2*pi - Theta;
    end
else
    if n(2) >= 0
        Alfa = pi - Theta;
    else
        Alfa = pi + Theta;
    end
end
%%%%角度Beta
Beta = acos( n(3) / norm( n ) );
%%%%计算旋转矩阵
Rzn = [  cos(-Alfa)   -sin(-Alfa)    0;
         sin(-Alfa)    cos(-Alfa)    0;
                0              0     1];
Ryn = [  cos(-Beta)     0     sin(-Beta);
                0       1             0 ;
        -sin(-Beta)     0     cos(-Beta)];
Rn = Ryn * Rzn;
Dat(:,1:3) = ( Rn * Pnt(:,1:3)' )';

%-------------------------------------------
%%%% 2.0 矢量n1
tmp = norm( n1(1:2) ); 
if tmp == 0,
   tmp = 0.0000001;
end
Theta = asin( abs( n1(2) ) / tmp );
%%%%角度Alfa
Alfa1 = 0;
if n1(1) >= 0
    if n1(2) >= 0
        Alfa1 = Theta;
    else
        Alfa1 = 2*pi - Theta;
    end
else
    if n1(2) >= 0
        Alfa1 = pi - Theta;
    else
        Alfa1 = pi + Theta;
    end
end
%%%%角度Beta
Beta1 = acos( n1(3) / norm( n1 ) );
%%%%计算旋转矩阵
Rzn1 = [  cos(-Alfa1)   -sin(-Alfa1)    0;  
          sin(-Alfa1)    cos(-Alfa1)    0;  
                  0              0      1];
Ryn1 = [  cos(-Beta1)   0   sin(-Beta1);
                   0    1            0 ;
         -sin(-Beta1)   0   cos(-Beta1)];
Rn1 = Ryn1 * Rzn1;
Dat(:,4:6) = ( Rn1 * Pnt(:,4:6)' )';
%%显示修正之后的数据
subplot(2,2,2,'V6'), plot3( Dat(:,1), Dat(:,2), Dat(:,3), 'o', Dat(:,4), Dat(:,5), Dat(:,6), '*');   grid;
xlabel('x');   ylabel('y');    zlabel('z');   title('共平面验证');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%第三步：估计两个点集之间的旋转量( 最小二乘法 )，确定旋转矩阵，并加以纠正（ 绕z轴旋转 ）
%--------------------------------------------------------------------------
A = [ Dat(1,1),    Dat(1,2),     1,     0;
      Dat(1,2),   -Dat(1,1),     0,     1;
      Dat(2,1),    Dat(2,2),     1,     0;
      Dat(2,2),   -Dat(2,1),     0,     1;
      Dat(3,1),    Dat(3,2),     1,     0;
      Dat(3,2),   -Dat(3,1),     0,     1];
B = [ Dat(1,4);   %xc1
      Dat(1,5);   %yc1
      Dat(2,4);   %xc2
      Dat(2,5);   %yc2
      Dat(3,4);   %xc3
      Dat(3,5)];  %yc3
xx = inv( A'* A ) * A'* B;        %最小二乘
xx = xx / norm( xx(1:2) );        %归一化约束
Rz2 = [  xx(1)    xx(2)    0;
        -xx(2)    xx(1)    0;
            0        0     1];
Dat2(:,4:6) = Dat(:,4:6);
Dat2(:,1:3) = ( Rz2 * Dat(:,1:3)' )';
subplot(2,2,3,'V6'), plot3( Dat2(:,1), Dat2(:,2), Dat2(:,3), 'o', Dat2(:,4), Dat2(:,5), Dat2(:,6), '*');   grid;
xlabel('x');   ylabel('y');    zlabel('z');   title('共平面、消除旋转的验证');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%第四步：对世界坐标系实施综合的坐标转换
%--------------------------------------------------------------------------
Rotate = Rzn1'* Ryn1'* Rz2 * Ryn * Rzn;
Dat3(:,1:3) = ( Rotate * Pnt(:,1:3)' )';
Dat3(:,4:6) = Pnt(:,4:6);
%估计平移量
Translate = [0 0 0];
for i=1:3,
    Translate = Translate + ( Dat3(i,4:6) - Dat3(i,1:3) );
end
Translate = Translate / 3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%验证：对最终求得的相似变换进行验证
DataOut(:,1:3) = ( Rotate * Pnt(:,1:3)' + repmat( Translate',1,3) )';
DataOut(:,4:6) = Pnt(:,4:6);
subplot(2,2,4,'V6'), plot3( DataOut(:,1), DataOut(:,2), DataOut(:,3), 'o', DataOut(:,4), DataOut(:,5), DataOut(:,6), '*');   grid;   
xlabel('x');   ylabel('y');    zlabel('z');   title('相似变换的结果');