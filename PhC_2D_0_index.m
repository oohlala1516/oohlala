clear; clc; epssys=1.0e-6; %设定一个最小量，避免系统截断误差或除0错误
global NG G f Nkpoints
Nkpoints=50; %每个方向上取的点数，
NrSquare=10;
NG =(2*NrSquare+1)^2; % NG is the number of the G value
G = zeros(NG,2);

h=pi/180;
%use radius as independant variable to confine angle
a=1;
Rc=0.2*a;  %input the radius
th1=acos(-(2*Rc^2/a^2-1));% th2=acos((2*Rc^2/a^2-1));
%theta=pi/3:h:pi/2;
theta=th1:h:pi/2;
% theta=th1+(length(theta)-3)*h
%theta=pi/2-10*h:h:pi/2+10*h;
length(theta)

Lset=[];
Iset=[];
Inext=[];
ITdiff2=[];
ITdiff3=[];
ITdiff4=[];
ITdiff5=[];
void=[];
s=1;
for tt = theta
    a=1; a1=a*[1 0]; a2=a*[cos(tt) sin(tt)];
    b1=2*pi/a*[1 -cot(tt)]; b2=2*pi/a*[0 csc(tt)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %晶格的参数
    epsa = 12.5; %inner
    epsb = 1; %outer
    %Pf = Ac/Au 填充率
    Au =a^2*sin(tt); %二维格子原胞面积
    Ac = pi*(Rc)^2; %介质柱横截面积
    Pf=Ac/Au;
    if Pf>1
        error('Pf>1')
    end

    %construct the G list
    i = 1;
    for l = -NrSquare:NrSquare
        for m = -NrSquare:NrSquare
            G(i,:)=l*b1+m*b2;
            i = i+1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %生成k空间中的f(Gi-Gj)的值，i,j 从1到NG。
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f=zeros(NG,NG);
    for i=1:NG
        for j=1:NG
            Gij=norm(G(i,:)-G(j,:));
            if (Gij < epssys)
                f(i,j)=(1/epsa)*Pf+(1/epsb)*(1-Pf);
                %f(i,j)=1/(epsa*Pf+epsb*(1-Pf));
            else
                f(i,j)=(1/epsa-1/epsb)*Pf*2*besselj(1,Gij*Rc)/(Gij*Rc); %
            end;
        end;
    end;
    
    %高对称度顶点的坐标
    I=[0 0];
    X=[0 pi/a*csc(tt)];
    P=[-pi/a pi/a*cot(tt)];
    M=[pi/a*(csc(tt)-cot(tt)-tan(tt))*cot(tt) pi/a*csc(tt)];
    N=[pi/a*csc(tt)*tan(tt/2) pi/a*csc(tt)];
    T=[2*pi/a*csc(tt)*sin(tt/2)*cos(tt/2) 2*pi/a*csc(tt)*sin(tt/2)*sin(tt/2)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %路径I-T-N-I-X-M-I-P-M;减少计算量，未计算不经过I点的路径
    [IT_TE_eig]=getEigValue(I,T,1); %/(2*pi/a)
    [NI_TE_eig]=getEigValue(N,I,1);
    [IX_TE_eig]=getEigValue(I,X,1);
    [MI_TE_eig]=getEigValue(M,I,1);
    [IP_TE_eig]=getEigValue(I,P,1);
    
    %计算L值
    a=(max([IT_TE_eig(2,2),IT_TE_eig(2,3),IT_TE_eig(2,4)])-min([IT_TE_eig(2,2),IT_TE_eig(2,3),IT_TE_eig(2,4)]));
    b=(max([NI_TE_eig(Nkpoints-1,2),NI_TE_eig(Nkpoints-1,3),NI_TE_eig(Nkpoints-1,4)])-min([NI_TE_eig(Nkpoints-1,2),NI_TE_eig(Nkpoints-1,3),NI_TE_eig(Nkpoints-1,4)]));
    c=(max([IX_TE_eig(2,2),IX_TE_eig(2,3),IX_TE_eig(2,4)])-min([IX_TE_eig(2,2),IX_TE_eig(2,3),IX_TE_eig(2,4)]));
    d=(max([MI_TE_eig(Nkpoints-1,2),MI_TE_eig(Nkpoints-1,3),MI_TE_eig(Nkpoints-1,4)])-min([MI_TE_eig(Nkpoints-1,2),MI_TE_eig(Nkpoints-1,3),MI_TE_eig(Nkpoints-1,4)]));
    e=(max([IP_TE_eig(2,2),IP_TE_eig(2,3),IP_TE_eig(2,4)])-min([IP_TE_eig(2,2),IP_TE_eig(2,3),IP_TE_eig(2,4)]));
    g=(max([IT_TE_eig(1,2),IT_TE_eig(1,3),IT_TE_eig(1,4)])-min([IT_TE_eig(1,2),IT_TE_eig(1,3),IT_TE_eig(1,4)]));
    
    L=a+b+c+d+e-5*g;
    Lset=[Lset,L];
    
    h=a+b+c+d+e;
    Inext=[Inext,h];
    Iset=[Iset,g];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     d2=IT_TE_eig(1,2)-IT_TE_eig(2,2);
%     d3=IT_TE_eig(1,3)-IT_TE_eig(2,3);
%     d4=IT_TE_eig(1,4)-IT_TE_eig(2,4); 
%     d5=max([IT_TE_eig(2,2),IT_TE_eig(2,3),IT_TE_eig(2,4)])-min([IT_TE_eig(2,2),IT_TE_eig(2,3),IT_TE_eig(2,4)])-max([IT_TE_eig(1,2),IT_TE_eig(1,3),IT_TE_eig(1,4)])-min([IT_TE_eig(1,2),IT_TE_eig(1,3),IT_TE_eig(1,4)])
    
    %v=max([IT_TE_eig(1,2),IT_TE_eig(1,3),IT_TE_eig(1,4)])-min([IT_TE_eig(1,2),IT_TE_eig(1,3),IT_TE_eig(1,4)]);
%     v=IT_TE_eig(1,4)-IT_TE_eig(1,2);
%     
%     ITdiff2=[ITdiff2,d2];
%     ITdiff3=[ITdiff3,d3];
%     ITdiff4=[ITdiff4,d4];
%     ITdiff5=[ITdiff5,d5];
%     void=[void,v];
    
    s %进程指示器
    s=s+1;
    G = zeros(NG,2);
end
Angle=theta*180/pi;
%%%%%%%%%%角度化%%%%%%%%%%%%%%%%%%
figure(2)
title('定义的L值')
scatter(Angle,Lset)
xlabel('Angle/\circ');
ylabel('L')
[m,p]=max(Lset);
target=[theta(p),Angle(p)];
target

figure(3)
title('gamma点处带隙宽度')
scatter(Angle,Iset)
figure(4)
title('gamma点两侧5个点带隙宽度的和')
scatter(Angle,Inext)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% plot(Angle,ITdiff2,'r',Angle,ITdiff3,'b',Angle,ITdiff4,'g',Angle,ITdiff5,'k')
% figure(3)
% scatter(Angle,void)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%获取本征值
function [eigValue]=getEigValue(k_begin,k_end,mode)
%mode :0 for TH 1 for TE
global NG G f Nkpoints
THETA=zeros(NG,NG); %待解的TE波矩阵
stepsize=0:1/(Nkpoints-1):1; %每个方向上的步长
TX_TE_eig = zeros(Nkpoints,NG); 
for n=1:Nkpoints %scan the 10 points along the direction
    %fprintf(['\n k-point:',int2str(n),'of',int2str(Nkpoints),'.\n']);
    step = stepsize(n)*(k_end-k_begin)+k_begin;  % get the k
    for i=1:(NG-1)   % G
        for j=(i+1):NG % G'
            kGi = step+G(i,:); %k+G
            kGj = step+G(j,:); %K+G'
            switch mode
                case 0  %TE mode
                    THETA(i,j)=f(i,j)*dot(kGi,kGj); %(K+G)(K+G')f(G-G')
                case 1 %TH mode
                     THETA(i,j)=f(i,j)*norm(kGi)*norm(kGj);
            end
            THETA(j,i)=conj(THETA(i,j)); 
        end
    end    
    for i=1:NG
        kGi = step+G(i,:);
        THETA(i,i)=f(i,i)*norm(kGi)*norm(kGi); 
    end
    eigValue(n,:)=sort(sqrt(eig(THETA))).';
end
end
