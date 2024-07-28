% --------- xingbo.m --------


clc;clear;close all;

%A母线行波
load xingbo.mat; % 加载故障发生模块数据.mat


m=ans.data;

ua=m(3501:3900,1)-m(1501:1900,1);
ia=m(3501:3900,4)-m(1501:1900,4);
ub=m(3501:3900,2)-m(1501:1900,2);
ib=m(3501:3900,5)-m(1501:1900,5);
uc=m(3501:3900,3)-m(1501:1900,3);
ic=m(3501:3900,6)-m(1501:1900,6);

Q=1/3*[ 1 1 1
        1 -1 0 
        1 0 -1;];
im=Q*[ ia ib ic]';
um1=Q(2,:)*[ ua ub uc]';%a模分量
im1=Q(2,:)*[ ia ib ic]'; % Clarke变换
Lm1=0.8984e-3;
Cm1=12.94e-9;
Zcm1=sqrt(Lm1/Cm1); % ???迹
uf=(um1+im1*Zcm1);%正向行波a模分量
ur=(um1-im1*Zcm1); %反向行波a模分量

im2=im(3,:);
im3=(ib'-ic')/3;
im0=im(1,:);

if (sum(abs(im0)))/400 >= 0.05
    a=((sum(abs(im1)-abs(im2)))/400<1);
    b=(abs(sum(im3))/400<=1);
    if (a==1)&&(b==1)
         disp('A相接地故障');
    end
    %if im1-im2
    a=(sum((abs(im1)-abs(im3)))/400<1);
    b=(sum(abs(im2))/400<=1);
    if (a==1)&&(b==1)
        disp('B相接地故障');
    end
    a=(sum((abs(im2)-abs(im3)))/400<1);
    b=(sum(abs(im1))/400<=1);
    if (a==1)&&(b==1)
        disp('C相接地故障');
    end
    a1=sum(abs(im1)-abs(im2))/400;
    b1=sum(abs(im1)-abs(im3))/400;
    if (a1>0.5)&&(b1>0.5)
        disp('AB相接地故障');
    end
    a1=sum(abs(im3)-abs(im1))/400;
    b1=sum(abs(im3)-abs(im2))/400;
    if (a1>0.5)&&(b1>0.5)
        disp('BC相接地故障');
    end
    a1=sum(abs(im2)-abs(im1))/400;
    b1=sum(abs(im2)-abs(im3))/400;
    if (a1>0.5)&&(b1>0.5)
        disp('AC相接地故障');
    end
end

if abs(sum(im0))/400 <= 0.05
    a2=abs(im1)-abs(2*im2);
    b2=abs(im1)-abs(2*im3);
    a2=sum(a2)/400;
    b2=sum(b2)/400;
    if (abs(a2)<1)&&(abs(b2)<1)
        disp('AB相间短路故障');
    end
    a2=abs(im3)-abs(2*im2);
    b2=abs(im3)-abs(2*im1);
    a2=sum(a2)/400;
    b2=sum(b2)/400;
    if (abs(a2)<1)&&(abs(b2)<1)
        disp('BC相间短路故障');
    end
    a2=abs(im2)-abs(2*im1);
    b2=abs(im2)-abs(2*im3);
    a2=sum(a2)/400;
    b2=sum(b2)/400;
    if (abs(a2)<1)&&(abs(b2)<1)
        disp('AC相间短路故障');
    end
end



uf1=uf';
ur1=ur';
t1=0:10:3990;
t=t1';
plot(t,uf1,'r',  t,ur1,'b--');
xlabel('t/us');ylabel('u/V');
legend('正向行波','反向行波','location','northwest');

[c1,l1]=wavedec(uf,8,'db8');%多尺度一维小波分解函数。其中C为分解结构变量，L为个分解结构以及原始信号长度变量，X为原始信号，N为分解层度，’wname’为小波类型。
db1=wrcoef('d',c1,l1,'db8',1);%一维小波系数单支重构函数。其中，当’type’ = ’a’时重构近似序列，当’type’ = ’d’时重构高频系数；C,L,’wname’,N含义同上
figure
plot(db1)
x=0;
c2=db1;
for i=1:1:(400-2)
    if (c2(1,i+1)>=c2(1,i)) && (c2(1,i+1)>=c2(1,i+2))
        x=x+1;
        ca(1,x)=i+1;
        ca(2,x)=c2(1,i+1);
    end
end
ca(2,:)=ca(2,:)/max(ca(2,:));
[m,n]=size(ca);
x=0;
for i=1:1:(n-2)
    if (ca(2,i+1)>=ca(2,i)) && (ca(2,i+1)>=ca(2,i+2))
        x=x+1;
        cb(1,x)=ca(1,i+1);
        cb(2,x)=ca(2,i+1);
    end
end
x=0;
[m,n]=size(cb);
for i=1:1:(n-1)
    
    if (cb(2,i)>=0.1)%设定阈值，确定某一区域内的极大值
        x=x+1;
        
            t2(1,x)=cb(1,i);
            t2(2,x)=cb(2,i);
   
    end
end

t3=unique(t2(1,:));
v=1/sqrt(Lm1*Cm1);%波速
L=(t3(1,2)-t3(1,1))*(10^(-5))*v/(2);%定位的距离v
fprintf('故障定位距离为：%6.5f\n',L)
        
%B母线行波
load xingbo1.mat;

m=ans.data;

ua=m(3501:3900,1)-m(1501:1900,1);
ia=m(3501:3900,4)-m(1501:1900,4);
ub=m(3501:3900,2)-m(1501:1900,2);
ib=m(3501:3900,5)-m(1501:1900,5);
uc=m(3501:3900,3)-m(1501:1900,3);
ic=m(3501:3900,6)-m(1501:1900,6);

Q=1/3*[ 1 1 1
        1 -1 0 
        1 0 -1;];


um11=Q(2,:)*[ ua ub uc]';%a模分量
im11=Q(2,:)*[ ia ib ic]'; % ????Clarke?任?????????????????

Lm1=0.8984e-3;
Cm1=12.94e-9;
Zcm1=sqrt(Lm1/Cm1); % ???迹
uf1=(um11+im11*Zcm1);%正向行波a模分量
ur1=(um11-im11*Zcm1); %反向行波a模分量


    
    
    
