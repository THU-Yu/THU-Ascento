clear;
mw=2;%车轮质量
r=0.02;%车轮半径 m
Iw=2e-4;%车轮转动惯量
mp=3;%车体质量
h=0.04;%车体车轴距 m
u=0.2;%摩擦系数
Ips=5e-3;%车体水平转动惯量
Ipz=0.049  * 0.000001;%竖直转动惯量
D=0.025;%车轴距 m
g=9.8;% m/s^2 重力加速度
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
plot(t,y,'g');gtext("偏角theta");
xlabel("时间/s");
ylabel("角度/rad");
hold off;