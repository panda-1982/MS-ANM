clear;
index=1;
rmse_fi_1=0;
rmse_fi_2=0;
mse_fi_1=zeros(1,51);
mse_fi_2=zeros(1,51);
m=9; % 阵元数
p=3; % 信号数
%俯仰角为入射方向与圆阵轴线之间的夹角,方位角为x轴到入射线在圆阵上投影的夹角
st1=90; % 俯仰角
dt1=-20; % 方位角
st2=90;
dt2=-5;
st3=90;
dt3 = 120;

st=[st1;st2;st3];
dt=[dt1;dt2;dt3];
degrad=pi/180;
lamda=0.3;
k0=2*pi/lamda;
radius=2*lamda;
M1=19;
T=361;
% 圆半径为半波长
A=zeros(m,T);
A1=zeros(m,T);
k=(0:m-1)';
for t=0:(T-1)
    A(:,t+1)=exp(i*2*pi*(radius/lamda)*(sin(st1*degrad)*cos((t-1)*degrad)*cos(2*pi*k/m)+sin(st1*degrad)*sin((t-1)*degrad)*sin(2*pi*k/m)));
end

g=zeros(m,360);
for t=1:m
    g(t,:)=ifftshift(ifft(A(t,:),360));
end
G=g(:,181-M1:181+M1);
A=zeros(m,p);
k=(0:m-1)';
for t=1:p
  A(:,t)=exp(i*2*pi*(radius/lamda)*(sin(st(t)*degrad)*cos(dt(t)*degrad)*cos(2*pi*k/m)+sin(st(t)*degrad)*sin(dt(t)*degrad)*sin(2*pi*k/m)));
end
load sig3
load crb
snr = 0;
sn=[snr;snr;snr];
nn=size(s1,1);
tt=1:nn;
S=[s1(tt),s2(tt),s3(tt)]';

Ps=S*S'/nn;
ps=diag(Ps);
refp=2*10.^(sn/10);
tmp=sqrt(refp./ps);
SS=diag(tmp)*S;
iter = 300;
mse_GL = zeros(1,10);
mse_MST1 = zeros(1,10);
mse_SPA1 = zeros(1,10);
mse_music = zeros(1,10);
rmse_GL =0;
rmse_music =0;
rmse_MST1 =0;
rmse_SPA1 =0;
delta_index = 1;
num_snapshot = 200;
est_freq = zeros(3,100);
for snr= -10:2:20
    for num = 1:iter
        Xt=A*(SS(:,1:num_snapshot));
        X = awgn(Xt,snr,'measured');
        [freq, coef] = MS_MUSIC(X,G,p);
        rmse_music=rmse_music+(sum(abs(freq-dt).^2));
        [freq, coef] = MST_GL_ANM(X,G,p);
        rmse_GL=rmse_GL+(sum(abs(freq-dt).^2));
        [freq, coef] = MST_ANM_LD(X,G,p);
        rmse_MST1=rmse_MST1+(sum(abs(freq-dt).^2));
        [freq, coef] = MS_SPA(X, G, p);
        rmse_SPA1=rmse_SPA1+(sum(abs(freq-dt).^2));
    end
    mse_music(delta_index) = sqrt(rmse_music/iter);
    mse_GL(delta_index) = sqrt(rmse_GL/iter);
    mse_MST1(delta_index) = sqrt(rmse_MST1/iter);
    mse_SPA1(delta_index) = sqrt(rmse_SPA1/iter);
    rmse_music =0;
    rmse_GL =0;
    rmse_MST1 =0;
    rmse_SPA1 =0;
    delta_index = delta_index + 1;
end
semilogy(-10:2:20,mse_music, '-*');
hold on;
semilogy(-10:2:20,mse_GL,'->');
hold on;
semilogy(-10:2:20,mse_MST1,'-*');
hold on;
semilogy(-10:2:20,mse_SPA1,'-p');
grid on;
semilogy(-10:2:20,sum(crb,1),'-k');
legend('MUSIC','GLE-ANM','Low Dimensional MS-ANM','SPICE-GL','CRB')
ylabel('RMSE (\circ)')
xlabel('SNR(dB)')
