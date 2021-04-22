gravity = 9.81; % meters per second squared
mass = 0.2; % kilograms
spring_constant = 10; % newtons per meter
unweighted_length = 1; % meters
initial_displacement = 0.3; % meters
weight_displacement = -(mass*gravity)/spring_constant; % meters
length = unweighted_length + init_displacement + ...
    weight_displacement; % meters
gravitational_force = mass * gravity; % newtons
restoring_force = -spring_constant * initial_displacement; % newtons

t = 10; % seconds
dt = .001; % seconds

time_array = 1:1:(1/dt)*t;
length_array = zeros(1, (1/dt)*t);
displacement_array = zeros(1, (1/dt)*t);
velocity_array = zeros(1, (1/dt)*t);
acceleration_array = zeros(1, (1/dt)*t);

acceleration_array(1, 1) = (gravitational_force + restoring_force) / ... 
    mass; % meters per second squared 

for time_index = 1:1:(1/dt)*t
        if time_index == 1
            velocity_array(1, time_index) = 0;
            length_array(1, time_index) = length;
            acceleration_array(1, time_index) = (gravitational_force + ... 
                restoring_force) / mass; % meters per second squared 
        else
            velocity_array(1, time_index) = ... 
                velocity_array(1, time_index - 1) + ... 
                acceleration_array(1, time_index - 1)  * ... 
                dt; % meters per second
            length_array(1, time_index) = ... 
                length_array(1, time_index - 1) + ...
                velocity_array(1, time_index - 1) * dt; % meters
            displacement_array(1, time_index) = ...
                length_array(1, time_index) - unweighted_length; % meters
            restoring_force = -spring_constant * ...
                displacement_array(1, time_index); % newtons
            acceleration_array(1, time_index) = ...
                (gravitational_force + restoring_force) / ... 
                mass; % meters per second squared 
        end
end

plot(time_array, length_array) 
   