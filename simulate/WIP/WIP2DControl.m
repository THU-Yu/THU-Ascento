clear;
mw=2;%��������
r=0.02;%���ְ뾶 m
Iw=2e-4;%����ת������
mp=3;%��������
h=0.04;%���峵��� m
u=0.2;%Ħ��ϵ��
Ips=5e-3;%����ˮƽת������
Ipz=0.049  * 0.000001;%��ֱת������
D=0.025;%����� m
g=9.8;% m/s^2 �������ٶ�
%%
a = Ips - ((mp * h * r)^2)/(Iw+(mp+mw)*r*r);
b = -1 - (mp * h * r)/(Iw+(mp+mw)*r*r);
c = -mp*g*h;
A = [0,1;-c/a,0];
B = [0;b/a];
C = [1,0];
D = [];
p_uncontrol = eig(A)
Qk = ctrb(A,B)
Qc = obsv(A,C)
rk = rank(Qk)
rc = rank(Qc)
%%
%p = [-2.01+2.23i,-2.01-2.23i];%[-41.5895, -133.9830];%[-4.3893+51.6636i,-4.3893-51.6636i];%[-2.01+2.23i,-2.01-2.23i];
%K = acker(A,B,p);
K = lqr(A,B,eye(2,2),eye(1,1));
sys = ss(A-B*K,B,C,D);
x0=[0.1,0]';
t = 0:0.1:5;
[y,t,x]=initial(sys,x0,t);
hold on;
plot(t,y,'g');gtext("ƫ��theta");
xlabel("ʱ��/s");
ylabel("�Ƕ�/rad");
hold off;