function [P,Q]=CS_lapla_matrix_Neumann_bc(n,h)
P=gallery('tridiag',n+2,1/10,1,1/10);
Q=-6/5*gallery('tridiag',n+2);
Q(1,1)=-2681/480;
Q(n+2,n+2)=-2681/480;
Q(1,2)=23/3;
Q(n+2,n+1)=23/3;
Q(1,3)=-113/40;
Q(n+2,n)=-113/40;
Q(1,4)=13/15;
Q(n+2,n-1)=13/15;
Q(1,5)=-59/480;
Q(n+2,n-2)=-59/480;
Q=-Q/h^2;


