clear all
close all
N=160
F=rand(N+2,N+2);
U=laplacefft(F,'NEUMANN');
mesh(U)