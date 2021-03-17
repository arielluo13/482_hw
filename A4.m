clear all; close all;
[images, labels] = mnist_parse('train-images-idx3-ubyte', 'train-labels-idx1-ubyte');
images = double(reshape(images,784,60000));
[m,n]=size(images);
mn=mean(images,2);
im=images-repmat(mn,1,n);
[U,S,V]=svd(im/sqrt(n-1),'econ');
%% Projection onto 3 V-modes
for label=0:9
    label_indices = find(labels == label);
    plot3(V(label_indices, 2), V(label_indices, 3), V(label_indices, 5),...
        'o', 'DisplayName', sprintf('%i',label), 'Linewidth', 2)
    hold on
end
xlabel('V-Mode 2'), ylabel('V-Mode 3'), zlabel('V-Mode 5')
title('Projection onto V-modes 2, 3, 5')
legend
set(gca,'Fontsize', 14)
%%
z1=S*V';
ind=find(labels==0);
zero=z1(:,ind);
ind=find(labels==1);
one=z1(:,ind);
ind=find(labels==2);
two=z1(:,ind);

n0=size(zero,2);
n1=size(one,2);


lab0=zero(1:10,1:n0);
lab1=one(1:10,1:n0);
lab2=two(1:10,1:n0);
m1=mean(lab1,2);
m0=mean(lab0,2);
m2=mean(lab2,2);
sw=0;
for k = 1:n0
    sw=sw+(lab0(:,k)-m0)*(lab0(:,k)-m0)';
end
for k = 1:n0
    sw=sw+(lab1(:,k)-m1)*(lab1(:,k)-m1)';
end
sb=(m0-m1)*(m0-m1)';
[V2,D]=eig(sb,sw);
[lambda, ind] = max(abs(diag(D)));
w = V2(:,ind);
w = w/norm(w,2);
v0=w'*lab0;
v1=w'*lab1;
if mean(v0)>mean(v1)
    w = -w;
    v0 = -v0;
    v1 = -v1;
end

%%
figure(2)
plot(v0,zeros(size(v0)),'ob','Linewidth',2)
hold on
plot(v1,ones(size(v1)),'dr','Linewidth',2)
title('0 and 1 values')
ylim([0 1.2])

%%
sort0 = sort(v0);
sort1 = sort(v1);
t0=length(sort0);
t1=1;
while sort0(t0) > sort1(t1)
    t0=t0-1;
    t1=t1+1;
end
threshold=(sort0(t0)+sort1(t1))/2;

%%
figure(3)
subplot(1,2,1)
histogram(sort0,30); hold on, plot([threshold threshold], [0 20],'r')
set(gca,'Xlim',[-10 5],'Ylim',[0 20],'Fontsize',14)
title('Digit 0')
subplot(1,2,2)
histogram(sort1,30); hold on, plot([threshold threshold], [0 20],'r')
set(gca,'Xlim',[-5 7],'Ylim',[0 20],'Fontsize',14)
title('Digit 1')

%%
[testim, testlb] = mnist_parse('t10k-images-idx3-ubyte', 't10k-labels-idx1-ubyte');
testim = double(reshape(testim,784,10000));

%%
ppro=U(:,1:10)'*testim;
%%
ind=find((testlb==1)|(testlb==0));
test=testim(:,ind);
testd=ppro(:,ind);
testl=testlb(ind);
%%
pval=w'*testd;
resVec=(pval>threshold);
err = abs(resVec - testl);
errNum = sum(err);
sucRate = 1 - errNum/2115;
%% 
k = 1;
figure(4)
for j = 1:2115
    if resVec(j) ~= testl(j)
        S = reshape(test(:,j),28,28);
        subplot(1,2,k)
        imshow(S)
        k = k+1;
    end
end
%%
sw3=sw;
for k = 1:n0
    sw3=sw3+(lab2(:,k)-m2)*(lab2(:,k)-m2)';
end
m12=mean([lab0 lab1 lab2],2);
sb3=(m0-m12)*(m0-m12)'+(m1-m12)*(m1-m12)'+(m2-m12)*(m2-m12)';
[V3, D3] = eig(sb3,sw3);
[lambda, ind] = max(abs(diag(D3)));
w3 = V3(:,ind);
w3 = w3/norm(w3,2);
v0_3 = w3'*lab0;
v1_3 = w3'*lab1;
v2_3 = w3'*lab2;

%%
figure(5)
plot(v0_3,zeros(size(v0_3)),'ob','Linewidth',2)
hold on
plot(v1_3,ones(size(v1_3)),'dr','Linewidth',2)
hold on
plot(v2_3,2*ones(size(v2_3)),'sk','Linewidth',2)
title('0, 1 and 2 values')
ylim([0 2.2])

%%
sort0_3 = sort(v0_3);
sort1_3 = sort(v1_3);
sort2_3 = sort(v2_3);
t0_3=length(sort0_3);
t1_3=1;
while sort0_3(t0_3) > sort1_3(t1_3)
    t0_3=t0_3-1;
    t1_3=t1_3+1;
end
threshold1=(sort0_3(t0_3)+sort1_3(t1_3))/2;
t1_3=length(sort1_3);
t2_3=1;
while sort1_3(t1_3) > sort2_3(t2_3)
    t1_3=t1_3-1;
    t2_3=t2_3+1;
end
threshold2=(sort1_3(t1_3)+sort2_3(t2_3))/2;
threshold_3=(threshold1+threshold2)/2;
%%
figure(6)
subplot(1,3,1)
histogram(sort0_3,30); hold on, plot([threshold_3 threshold_3], [0 20],'r')
set(gca,'Xlim',[-10 5],'Ylim',[0 20],'Fontsize',14)
title('Digit 0')
subplot(1,3,2)
histogram(sort1_3,30); hold on, plot([threshold_3 threshold_3], [0 20],'r')
set(gca,'Xlim',[-5 7],'Ylim',[0 20],'Fontsize',14)
title('Digit 1')
subplot(1,3,3)
histogram(sort2_3,30); hold on, plot([threshold_3 threshold_3], [0 20],'r')
set(gca,'Xlim',[-5 7],'Ylim',[0 20],'Fontsize',14)
title('Digit 2')
%%
ppro_3=U(:,1:10)'*testim;
%%
ind=find((testlb==1)|(testlb==0)|(testlb==2));
test3=testim(:,ind);
testd3=ppro_3(:,ind);
testl3=testlb(ind);
pval=w3'*testd3;
resVec3=(pval>threshold_3);
%%
k = 1;
figure(6)
for j = 1:3147
    if resVec3(j) ~= testl3(j)
        S = reshape(test3(:,j),28,28);
        subplot(1,3,k)
        imshow(S)
        k = k+1;
    end
end
