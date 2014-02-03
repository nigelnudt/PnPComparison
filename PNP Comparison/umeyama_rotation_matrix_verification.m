clear
close all
syms a11 a12 a13 a21 a22 a23 a31 a32 a33 b11 b12 b13 b21 b22 b23 b31 b32 b33 
S1=[1,0,0;0,1,0;0,0,-1];
S2=[1,0,0;0,-1,0;0,0,1];
matrix_a=[a11 a12 a13; a21 a22 a23 ;a31 a32 a33 ];
matrix_b=[b11 b12 b13; b21 b22 b23 ;b31 b32 b33 ];
matrix_as=matrix_a*S1
matrix_asb=matrix_as*matrix_b';
matrix_abs=matrix_a*matrix_b'*S2
matrix_diff=matrix_asb-matrix_abs;
tic