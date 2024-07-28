% --------- xingbo.m --------


clc;clear;close all;

%Aĸ���в�
load xingbo.mat; % ���ع��Ϸ���ģ������.mat


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
um1=Q(2,:)*[ ua ub uc]';%aģ����
im1=Q(2,:)*[ ia ib ic]'; % Clarke�任
Lm1=0.8984e-3;
Cm1=12.94e-9;
Zcm1=sqrt(Lm1/Cm1); % ???��
uf=(um1+im1*Zcm1);%�����в�aģ����
ur=(um1-im1*Zcm1); %�����в�aģ����

im2=im(3,:);
im3=(ib'-ic')/3;
im0=im(1,:);

if (sum(abs(im0)))/400 >= 0.05
    a=((sum(abs(im1)-abs(im2)))/400<1);
    b=(abs(sum(im3))/400<=1);
    if (a==1)&&(b==1)
         disp('A��ӵع���');
    end
    %if im1-im2
    a=(sum((abs(im1)-abs(im3)))/400<1);
    b=(sum(abs(im2))/400<=1);
    if (a==1)&&(b==1)
        disp('B��ӵع���');
    end
    a=(sum((abs(im2)-abs(im3)))/400<1);
    b=(sum(abs(im1))/400<=1);
    if (a==1)&&(b==1)
        disp('C��ӵع���');
    end
    a1=sum(abs(im1)-abs(im2))/400;
    b1=sum(abs(im1)-abs(im3))/400;
    if (a1>0.5)&&(b1>0.5)
        disp('AB��ӵع���');
    end
    a1=sum(abs(im3)-abs(im1))/400;
    b1=sum(abs(im3)-abs(im2))/400;
    if (a1>0.5)&&(b1>0.5)
        disp('BC��ӵع���');
    end
    a1=sum(abs(im2)-abs(im1))/400;
    b1=sum(abs(im2)-abs(im3))/400;
    if (a1>0.5)&&(b1>0.5)
        disp('AC��ӵع���');
    end
end

if abs(sum(im0))/400 <= 0.05
    a2=abs(im1)-abs(2*im2);
    b2=abs(im1)-abs(2*im3);
    a2=sum(a2)/400;
    b2=sum(b2)/400;
    if (abs(a2)<1)&&(abs(b2)<1)
        disp('AB����·����');
    end
    a2=abs(im3)-abs(2*im2);
    b2=abs(im3)-abs(2*im1);
    a2=sum(a2)/400;
    b2=sum(b2)/400;
    if (abs(a2)<1)&&(abs(b2)<1)
        disp('BC����·����');
    end
    a2=abs(im2)-abs(2*im1);
    b2=abs(im2)-abs(2*im3);
    a2=sum(a2)/400;
    b2=sum(b2)/400;
    if (abs(a2)<1)&&(abs(b2)<1)
        disp('AC����·����');
    end
end



uf1=uf';
ur1=ur';
t1=0:10:3990;
t=t1';
plot(t,uf1,'r',  t,ur1,'b--');
xlabel('t/us');ylabel('u/V');
legend('�����в�','�����в�','location','northwest');

[c1,l1]=wavedec(uf,8,'db8');%��߶�һάС���ֽ⺯��������CΪ�ֽ�ṹ������LΪ���ֽ�ṹ�Լ�ԭʼ�źų��ȱ�����XΪԭʼ�źţ�NΪ�ֽ��ȣ���wname��ΪС�����͡�
db1=wrcoef('d',c1,l1,'db8',1);%һάС��ϵ����֧�ع����������У�����type�� = ��a��ʱ�ع��������У�����type�� = ��d��ʱ�ع���Ƶϵ����C,L,��wname��,N����ͬ��
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
    
    if (cb(2,i)>=0.1)%�趨��ֵ��ȷ��ĳһ�����ڵļ���ֵ
        x=x+1;
        
            t2(1,x)=cb(1,i);
            t2(2,x)=cb(2,i);
   
    end
end

t3=unique(t2(1,:));
v=1/sqrt(Lm1*Cm1);%����
L=(t3(1,2)-t3(1,1))*(10^(-5))*v/(2);%��λ�ľ���v
fprintf('���϶�λ����Ϊ��%6.5f\n',L)
        
%Bĸ���в�
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


um11=Q(2,:)*[ ua ub uc]';%aģ����
im11=Q(2,:)*[ ia ib ic]'; % ????Clarke?��?????????????????

Lm1=0.8984e-3;
Cm1=12.94e-9;
Zcm1=sqrt(Lm1/Cm1); % ???��
uf1=(um11+im11*Zcm1);%�����в�aģ����
ur1=(um11-im11*Zcm1); %�����в�aģ����


    
    
    
