% Jackson Prettyman, Matt Bass, Will Johnson
% May 10, 2021
% CA model of COVID infection accross a population

% figure dimensions
m = 10; %number of rows 
n = 10; %number of cols
x = 1:m;
y = 1:n; 

% creating grid
cell_value = zeros(m, n);
cell_value_extend = zeros(m + 2,n + 2);
observation_period = {1:number_iterations};
observation_period_extend = {1:number_iterations};

% simulation constants
days = 14; % number of days
number_iterations = days; % number of iterations 

% initial population values
total_popultion = m * n; % total number of people
initial_susceptible = 0; % total number of susceptible people
initial_infected = 0; % total number of infected people
initial_recovered = 0; % total number of recovered people

% population arrays
susceptible_array = zeros(1, number_iterations);
susceptible_array(1) = initial_susceptible;
infected_array = zeros(1, number_iterations);
infected_array(1) = initial_infected;
recovered_array = zeros(1, number_iterations);
recovered_array(1) = initial_recovered;

% constants
length_of_infection = 2; % length in days that infection lasts
recovery_rate = 1 / length_of_infection; % rate at which infected recovers 
transmission_constant = .50; % chance a susceptible becomes infected ...
                             % when interacting with 1 infected
                                             
% initialization probabilities
probability_initially_susceptible = 0.9; % probability to start susceptible 
probability_initially_infected = 0.1; % probability to start infected
probability_initially_phase_1 = 0.5; % probability if infected 
                                     % starts in phase 1

% neighbor function for use in infection simulation
% input: coordinates of grid extended
% output odds of infection given number of infected around susceptible
neighbor_sum = @(x, y, cell_value_ex) ...
    (1 - (1 - transmission_constant)^(...
    (cell_value_ex(x + 1, y) > 0 && cell_value_ex(x + 1, y) < 3) ...
    + (cell_value_ex(x - 1, y) > 0 && cell_value_ex(x - 1, y) < 3) ...
    + (cell_value_ex(x, y + 1) > 0 && cell_value_ex(x, y + 1) < 3) ...
    + (cell_value_ex(x, y - 1) > 0 && cell_value_ex(x ,y - 1) < 3)));

% initializing the grid and grid list
cell_value(x, y) = cell_value(x, y) + (rand > ...
                   probability_initially_susceptible);
cell_value_extend(2:m + 1, 2:n + 1) = cell_value;
observation_period{1} = cell_value;
observation_period_extend{1} = cell_value_extend;

for i = 2:number_iterations
cell_value = observation + period{i - 1};
cell_value_extend = observation_period_extend{i - 1};
%Alright so this is the set up for the loop to simulate the problem i am
%having is trying to figure out how to make sure the neighbor function only
%applies to the 0 values. I am pretty sure this is the use of conditionals
%but am not sure and not confident. One thought would be to find all index
%where values are 0 and use that to call on neighbor function in a loop.
%Check the conditional i have below for checking recovery rate once we
%have a functioning one I will remove all magic numbers
if cell_value(x, y) == 8
    cell_value(x, y) = 0;
end

cell_value(x, y) = cell_value(x, y)+(rand > neighbor_sum(x, y, ...
                   cell_value_extend));
cell_value_extend(2:m + 1,2:n + 1) = cell_value;
observation_period{i} = cell_value;
observation_period_extend{i} = cell_value_extend;
end

