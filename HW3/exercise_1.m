

%Constants set for all Methods
Cm=0.01; % Membrane Capcitance uF/cm^2
dt=0.04; % Time Step ms
t=0:dt:25; %Time Array ms
I=0.1; %External Current Applied
ENa=55.17; % mv Na reversal potential
EK=-72.14; % mv K reversal potential
El=-49.42; % mv Leakage reversal potential
gbarNa=1.2; % mS/cm^2 Na conductance
gbarK=0.36; % mS/cm^2 K conductance
gbarl=0.003 % mS/cm^2 Leakage conductance
V(1)=-60; % Initial Membrane voltage
m(1)=am(V(1))/(am(V(1))+bm(V(1))); % Initial m-value
n(1)=an(V(1))/(an(V(1))+bn(V(1))); % Initial n-value
h(1)=ah(V(1))/(ah(V(1))+bh(V(1))); % Initial h-value






%% Runge-Kutta Method
V(1)=-60; % Initial Membrane voltage
m(1)=am(V(1))/(am(V(1))+bm(V(1))); % Initial m-value
n(1)=an(V(1))/(an(V(1))+bn(V(1))); % Initial n-value
h(1)=ah(V(1))/(ah(V(1))+bh(V(1))); % Initial h-value

for i=1:length(t)-1 % Loop through each step until time is
finished

 %4 step method of Runge-Kutta
 K1=dt*HH(i,[V(i); n(i); m(i); h(i)]);
 k1=K1(1,1);n1=K1(2,1);m1=K1(3,1);h1=K1(4,1);% obtain 4 k
%%variables (V,m,n,h) from HH function

K2=dt*HH(i+(0.5*dt),[V(i)+(0.5*k1);n(i)+(0.5*n1);m(i)+(0.5*m1);h(i
)+(0.5*h1)]);
 k2=K2(1,1);n2=K2(2,1);m2=K2(3,1);h2=K2(4,1);

K3=dt*HH(i+(0.5*dt),[V(i)+(0.5*k2);n(i)+(0.5*n2);m(i)+(0.5*m2);h(i
)+(0.5*h2)]);
 k3=K3(1,1);n3=K3(2,1);m3=K3(3,1);h3=K3(4,1);
 K4=dt*HH(i+dt,[V(i)+k3;n(i)+n3;m(i)+m3;h(i)+h3]);
 k4=K4(1,1);n4=K4(2,1);m4=K4(3,1);h4=K4(4,1);

 %create next step for each variable
 V(i+1)=V(i)+1/6*(k1+2*k2+2*k3+k4);
 n(i+1)=n(i)+1/6*(n1+2*n2+2*n3+n4);
 m(i+1)=m(i)+1/6*(m1+2*m2+2*m3+m4);
 h(i+1)=h(i)+1/6*(h1+2*h2+2*h3+h4);
end
%set variables for graphing later
RK=V;
RKm=m;
RKn=n;
RKh=h;