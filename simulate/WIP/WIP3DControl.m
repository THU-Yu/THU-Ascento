clear;
mw=1;%��������
R=0.02;%���ְ뾶 m
Iw=1e-4;%����ת������
mp=6;%��������
h=0.04;%���峵��� m
u=0.2;%Ħ��ϵ��
Ips=5e-3;%������ֱת������
Ipz=2e-4;%ˮƽת������
D=0.05;%����� m
g=9.8;% m/s^2 �������ٶ�
%%
fm=(mp^2*R^4-2*(mp*h^2+Ips)*(mw*R^2+mp*R^2/2+Iw));
a1=-u/(2*mw*R^2+Iw+(2*R^2*Ipz)/(D^2));
a2=u*(mp*h^2+Ips)/fm;
a3=mp^2*h^2*R^2*g/fm;
a4=-mp*h*u/fm;
a5=-mp*g*h*(mw*R^2+mp*R^2/2+Iw)/fm;
b1=-R/(2*mw*D*R^2+D*Iw+2*R^2*Ipz/D);
b2=-R*(mp*h^2+Ips+mp*h*R)/fm;
b3=(mp*h*R+mw*R^2+mp*R^2/2+Iw)/fm;
%s-s����
A=[0 1 0 0 0 0;0 a1 0 0 0 0;0 0 0 1 0 0;0 0 0 a2 a3 0;0 0 0 0 0 1;0 1 0 a4 a5 0];
B=[0 0;b1 -b1;0 0;b2 b2;0 0;b3 b3];
C=eye(6,6);%[1 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0];
D=[];
p_uncontrol = eig(A)
Qk = ctrb(A,B)
Qc = obsv(A,C)
rk = rank(Qk)
rc = rank(Qc)
%%
%p = [-2.01+2.23i,-2.01-2.23i,-100,-101,-102,-103];
%K = place(A,B,p);
K = lqr(A,B,eye(6,6),eye(2,2))
%K = [-0.5,-0.02,-0.001,-1.5,-1.1412,-0.7071;0.5,0.02,-0.001,-1.5,-1.1412,-0.7071]
sys = ss(A-B*K,B,C,D);
x0=[0.1745,0,0,0,0.1,0]';
t = 0:0.1:10;
[y,t,x]=initial(sys,x0,t);
hold on;
plot(t,y(:,5),'g');gtext("ƫ��theta");
plot(t,y(:,2),'r');gtext("�ٶ�v");
plot(t,y(:,1),'b');gtext("ƫ��delta");
xlabel("ʱ��/s");
hold off;