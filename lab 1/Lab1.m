%Will Johsnon Lab 0 
%Spring simulation using RK4


g = -9.81; %This is the gravitational constant meters per second squared
k = 10; %This is the spring constant Newtens per meter
m_W = 0.2; %This is the mass of the weight in kilograms
u_L = 1; %This is the unweighted length of the spring in meters
init_displacement = 0.3; %This is the initial displacement for when the spring is hanging in meters
weight_disp = (m_W*g)/k; %this is the the amount of displacement caused by the weight
initial_springLen = u_L  + init_displacement + weight_disp; %This is where the spring start length wise
s = @(curr_len) (curr_len - u_L); %This is the function to calculate total displacment
grav_F = m_W*g; %gravitational force
rest_F = @(currlen) (-k*s(currlen)); %restorational force
a = @(F1,F2) ((F1+F2)/m_W); %function for acceleration

%these are simulation constants
t = 10; %total time in seconds
dt = 0.001; %change in time when calculating new point
num_iter = t/dt; %number of iteration

%derivtive functions
dPdt = @(time,pos,vel) (vel); %caclulates velocity based on change in displacment
dVdt = @(time,pos,vel) (a(grav_F,rest_F(pos))); %caclulates accel based on change in veloc

%Euler Method
disp = zeros(1,num_iter); 
veloc = zeros(1,num_iter);
disp(1) = initial_springLen;
veloc(1) = 0; 
for i = 2:num_iter
    rate_V = a(grav_F,rest_F(disp(i-1)));
    rate_P = veloc(i-1);
    amt_V = rate_V*dt;
    amt_P = rate_P*dt;
    veloc(i) = veloc(i-1) + amt_V;
    disp(i) = disp(i-1)+amt_P;
end

time = 0:dt:10;
plot(time(1:num_iter),disp);


%RK4 Method
disp = zeros(1,num_iter); 
veloc = zeros(1,num_iter);
disp(1) = initial_springLen;
veloc(1) = 0; 


for j = 2:num_iter
    del_P1 = dPdt(j,disp(j),veloc(j-1))*dt;
    del_V1 = dVdt(j,disp(j-1))*dt;
    del_P2 = dPdt(j+0.5*dt,disp(j)+0.5*del_P1,veloc(j-1)+0.5*del_V1)*dt;
    del_V2 = dVdt(j+0.5*dt,disp(j)+0.5*del_P1,veloc(j-1)+0.5*del_V1)*dt;
    del_P3 = dPdt(j+0.5*dt,disp(j)+0.5*del_P2,veloc(j-1)+0.5*del_V2)*dt;
    del_V3 = dVdt(j+0.5*dt,disp(j)+0.5*del_P2,veloc(j-1)+0.5*del_V2)*dt;
    del_P4 = dPdt(j+dt,disp(j)+del_P3,veloc(j-1)+del_V3)*dt;
    del_V4 = dVdt(j+dt,disp(j)+del_P3,veloc(j-1)+del_V3)*dt;
    change_disp = (del_P1 + 2*del_P2+2*del_P3 + del_P4)/6;
    change_veloc = (del_V1 + 2*del_V2+2*del_V3 + del_V4)/6;
    veloc(j) = veloc(j-1)+change_veloc;
    disp(j) = disp(j-1)+change_disp;
end
figure();
plot(time(1:num_iter),disp);
    
