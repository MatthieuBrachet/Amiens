%
%test frequences 3D
%
K=1:3;
one=ones(3,1);
X2=kron(one,K);
X2
pause
for i=1:3
X3(:,:,i)=X2;
end
pause
Y2=X2';
pause
for i=1:3
Y3(:,:,i)=Y2;
end
Z2=X2;
for i=1:3
    Z3(i,:,:)=X2
end


