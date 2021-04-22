% Matt Bass
%
% mnbass21_cs346_lab1.m
% CS346 -- Computational Modeling and Simulation
% Spring, 2021



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ ME:
 
% (ABOUT SCRIPT)
% This script is to run a Predator-Prey Model for where both the 

    
    
%(HOW TO RUN):

    
% These variables the basic running of the finite differential equations
% simulation
simLength = 10; %(int) This is the number of months to tun the sim
dt = 0.001; %(float) Time step of the sim in units of months 
numIterations = simLength/dt; %how many steps are in the simulation
velocity = 0;
mass = 0.2;
acceleration_due_to_gravity = 9.81;
weight = mass * acceleration_due_to_gravity;
spring_constant = 10;
unweighted_length = 1;
weight_displacement = weight / spring_constant;
displacement = 0.3;
length = unweighted_length + weight_displacement + displacement;

restoring_spring_force = -spring_constant * (length - unweighted_length);
total_force = weight + restoring_spring_force;
acceleration = total_force / mass;

dv = acceleration;
dL = velocity;


t = 0; % starting time is 0

% initilize lists to hold populations to plot and the time at each step
tArray = ones(1,numIterations)*t;
velocityArray = ones(1,numIterations)*velocity;
lengthArray = ones(1,numIterations)*length;

% simulation loop
for i = 2:numIterations
    t = i * dt;
    velocity = velocity + dv * dt;
    length = length + dL * dt;
    restoring_spring_force = -spring_constant * (length - unweighted_length);
    total_force = weight + restoring_spring_force;
    acceleration = total_force / mass;
    dL = velocity;
    dv = acceleration;
    %assign to arrays
    tArray(i) = t;
    velocityArray(i) = velocity;
    lengthArray(i) = length;
end
figure
plot(tArray,lengthArray)

% See write up for more
 
    
% Change the variables above in the simulation for different results 
% Then run the script with the run button in MatLab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


