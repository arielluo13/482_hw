clear all; close all;clc;
load cam1_4.mat;
load cam2_4.mat;
load cam3_4.mat;
%%
for j=1:300
    frame = vidFrames1_4(:,:,:,j);
    imshow(frame);
    axis image;
end
%%
image(vidFrames3_4(:,:,:,5));
axis image;
%%
ind=[];
lastF = vidFrames1_4(:,:,:,1);
[y,x,c,t]=size(vidFrames1_4);
%cam1_1
for f = 2:t
    diff=[];
    index=[];
    for i = 300:600
        for j = 200:400
            cl1 = lastF(j,i,:);
            cl2 = vidFrames1_4(j,i,:,f);
            if(all(cl2(:)>240))
                d = sum(cl2(:)-cl1(:));
                diff=[diff;d];
                index=[index;[i,j]];
            end
        end
    end
    s=find(diff==max(diff(:,1)));
    index=index(s,:);
    ind=[ind;[mean(index(:,1)),mean(index(:,2))]];
    lastF=vidFrames1_4(:,:,:,f);
end
ind=ind';
x1=ind(1,:);
y1=ind(2,:);

%%
 % cam2_1
[y,x,c,t]=size(vidFrames2_4);
ind=[];
lastF = vidFrames2_4(:,:,:,1);
for f = 2:t
    diff=[];
    index=[];
    for i = 250:500
        for j = 150:400
            cl1 = lastF(j,i,:);
            cl2 = vidFrames2_4(j,i,:,f);
            if(all(cl2(:)>235))
                d = sum(cl2(:)-cl1(:));
                diff=[diff;d];
                index=[index;[i,j]];
            end
        end
    end
    s=find(diff==max(diff(:,1)));
    index=index(s,:);
    ind=[ind;[mean(index(:,1)),mean(index(:,2))]];
    lastF=vidFrames2_4(:,:,:,f);
end
ind=ind';
x2=ind(1,:);
y2=ind(2,:);
%%
% cam3_1
[y,x,c,t]=size(vidFrames3_4);
lastF = vidFrames3_4(:,:,:,1);
ind=[];
 for f = 2:t
    diff=[];
    index=[];
    for i = 300:500
        for j = 150:350
            cl1 = lastF(j,i,:);
            cl2 = vidFrames3_4(j,i,:,f);
            if(all(cl2(:)>220))
                d = sum(cl2(:)-cl1(:));
                diff=[diff;d];
                index=[index;[i,j]];
            end
        end
    end
    s=find(diff==max(diff(:,1)));
    index=index(s,:);
    ind=[ind;[mean(index(:,1)),mean(index(:,2))]];
    lastF=vidFrames3_4(:,:,:,f);
 end
ind=ind';
x3=ind(1,:);
y3=ind(2,:);


 %%
 minlen=min(length(x1),length(x2));
 minlen=min(minlen,length(x3));
 result = [x1(:,1:minlen);y1(:,1:minlen);x2(:,1:minlen);y2(:,1:minlen);x3(:,1:minlen);y3(:,1:minlen)];
 

 [m,n]=size(result);
 mn=mean(result,2);
 result=result-repmat(mn,1,n);
 [U,S,V] = svd(result/sqrt(n-1));
 lambda=diag(S).^2;
 Y=U'*result;

xp=1:minlen;
plot(xp, Y(1,:), 'LineWidth', 1);
