% Jackson Prettyman, Matt Bass, Will Johnson
% May 10, 2021
% CA model of COVID infection accross a population

% figure dimensions
m = 10; %number of rows 
n = 10; %number of cols
x = 1:m;
y = 1:n; 

% creating grid
cell_value = zeros(m,n);
cell_value_extend = zeros(m+2,n+2);
observation_period = {1:num_iterations};
observation_period_extended = {1:num_iterations};

% initial population values
total_popultion = m*n; % total number of people
Susceptible = 0; % total number of susceptible people
Infected = 0; % total number of infected people
Recovered = 0; % total number of recovered people

% constants
recovery_rate = .5; % rate at which infected recovers 
transmission_constant = .50; %chance a susceptible becomes infected when ...
                             %when interacting with 1 infected
                                             
% initialization probabilities
probability_initially_susceptible = 0.9; %probability to start susceptible 
probability_initially_infected = 0.1; %probability to start infected
probability_initially_phase_1 = 0.5; %probability if infected, start in phase 1

%neighbor function for use in infection simulation
%input: coordinates of grid extended
%output odds of infection given number of infected around susceptible
neighbor_sum = @(x,y,cell_value_extended) (1-(1-transmission_constant)^ ...
    ((cell_value_extended(x+1,y) > 0 && cell_value_extended(x+1,y) < 3)+...
    (cell_value_extended(x-1,y) > 0 && cell_value_extended(x-1,y) < 3)+...
    (cell_value_extended(x,y+1) > 0 && cell_value_extended(x,y+1) < 3)+...
    (cell_value_extended(x,y-1) > 0 && cell_value_extended(x,y-1) < 3)));


%simulation constants
days = 14; % number of days
num_iterations = days; % number of iterations 

% initializing the grid and grid list
cell_value(x,y) = cell_value(x,y) + ...
    (rand>probability_initially_susceptible);
cell_value_extended(2:m+1,2:n+1) = cell_value;
observation_period{1} = cell_value;
observation_period_extended{1} = cell_value_extended;

for i = 2:num_iterations
cell_value = observation+period{i-1};
cell_value_extended = observation_period_extended{i-1};
%Alright so this is the set up for the loop to simulate the problem i am
%having is trying to figure out how to make sure the neighbor function only
%applies to the 0 values. I am pretty sure this is the use of conditionals
%but am not sure and not confident. One thought would be to find all index
%where values are 0 and use that to call on neighbor function in a loop.
%Check the conditional i have below for checking recovery rate once we
%have a functioning one I will remove all magic numbers
if cell_value(x,y) == 8
    cell_value(x,y) = 0;
end

cell_value(x,y) = cell_value(x,y)+(rand>neighbor_sum(x,y,cell_value_extended));
cell_value_extended(2:m+1,2:n+1) = cell_value;
observation_period{i} = cell_value;
observation_period_extended{i} = cell_value_extended;
end

