clear; clc; epssys=1.0e-6; %设定一个最小量，避免系统截断误差或除0错误
global NG G f Nkpoints

a=1; a1=a*[1 0]; a2=a*[1/2 3^(0.5)/2];
b1=2*pi/a*[1 -3^(0.5)/3];b2=2*pi/a*[0 2*3^(0.5)/3];
Nkpoints=10; %每个方向上取的点数，
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%定义晶格的参数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsa = 3.46; %inner &1.38
epsb = 1.5110; %outer %3.46
Pf = 0.1298; %Pf = Ac/Au 填充率，可根据需要自行设定 原来0.1257
Au =a^2*3^(0.5)/2; %二维格子原胞面积
Rc = (Pf *Au/pi)^(1/2); %介质柱截面半径
Ac = pi*(Rc)^2; %介质柱横截面积

%construct the G list
NrSquare = 10;
NG =(2*NrSquare+1)^2; % NG is the number of the G value
G = zeros(NG,2);
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
else
f(i,j)=(1/epsa-1/epsb)*Pf*2*besselj(1,Gij*Rc)/(Gij*Rc);
end;
end;
end;
T=(2*pi/a)*[epssys 0];
M=(2*pi/a)*[2/3 0];
X=(2*pi/a)*[1/2 3^(0.5)/6];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
[MT_TE_eig]=getEigValue(M,T,0);
[TX_TE_eig]=getEigValue(T,X,0);
[XM_TE_eig]=getEigValue(X,M,0);


kaxis = 0;
MTaxis = kaxis:norm(M-T)/(Nkpoints-1):(kaxis+norm(M-T));
kaxis = kaxis + norm(M-T);
TXaxis = kaxis:norm(T-X)/(Nkpoints-1):(kaxis+norm(T-X));
kaxis = kaxis + norm(T-X);
XMaxis = kaxis:norm(X-M)/(Nkpoints-1):(kaxis+norm(X-M));
kaxis = kaxis + norm(M-T);

Ntraject = 3;
figure (1)
hold on;
Nk=Nkpoints;
for k=1:4 %NG
for i=1:Nkpoints
EigFreq_TE(i+0*Nk) = MT_TE_eig(i,k)/(2*pi/a);
EigFreq_TE(i+1*Nk) = TX_TE_eig(i,k)/(2*pi/a);
EigFreq_TE(i+2*Nk) = XM_TE_eig(i,k)/(2*pi/a);

end
linehandle=plot(MTaxis(1:Nk),EigFreq_TE(1+0*Nk:1*Nk),'r',...
TXaxis(1:Nk),EigFreq_TE(1+1*Nk:2*Nk),'r',...
XMaxis(1:Nk),EigFreq_TE(1+2*Nk:3*Nk),'r',...
'LineWidth',1 );
set (linehandle, 'linesmoothing', 'on');
end



grid on;
xlabel('K-Space');
ylabel('Frequency(\omegaa/2\piC)');
axis([0 XMaxis(Nkpoints) 0 0.8]);
set(gca,'XTick',[MTaxis(1), MTaxis(Nkpoints),...
TXaxis(Nkpoints),XMaxis(Nkpoints)]);
tmixlabel = strvcat('M','T','X','M');
set(gca,'XTickLabel',tmixlabel);


%获取本征值

function [eigValue]=getEigValue(k_begin,k_end,mode)
%mode :0 for TE 1 for TH
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

