% Jackson Prettyman, Matt Bass, Will Johnson
% May 10, 2021
% CA model of COVID infection accross a population

% figure dimensions
rows= 10; %number of rows 
cols = 10; %number of cols
x = 1:rows;
y = 1:cols; 

% simulation constants
days = 14; % number of days
number_iterations = days; % number of iterations 

% creating grid
cell_value = zeros(rows, cols);
cell_value_extend = zeros(rows + 2, cols + 2);
observation_period = {1:number_iterations};
observation_period_extend = {1:number_iterations};

% constants
length_of_infection = 2; % length in days that infection lasts
length_of_immunity = 5
recovery_rate = 1 / length_of_infection; % rate at which infected recovers 
transmission_constant = .25; % chance a susceptible becomes infected ...
                             % when interacting with 1 infected

% neighbor function for use in infection simulation
% input: coordinates of grid extended
% output: odds of infection given number of infected around susceptible
probability_of_infection = @(x, y, cell_value_ex) ...
    ((transmission_constant) ...
    * ((cell_value_ex(x + 1, y) > 0 && cell_value_ex(x + 1, y) < ...
    length_of_infection + 1) ...
    + (cell_value_ex(x - 1, y) > 0 && cell_value_ex(x - 1, y) < ...
    length_of_infection + 1) ...
    + (cell_value_ex(x, y + 1) > 0 && cell_value_ex(x, y + 1) < ...
    length_of_infection + 1) ...
    + (cell_value_ex(x, y - 1) > 0 && cell_value_ex(x ,y - 1) < ...
    length_of_infection + 1)));

% initial population values
total_population = rows * cols; % total number of people

% initialization probabilities
probability_initially_susceptible = 0.9; % probability to start susceptible 

% initializing the grid and grid list
initial_infected = 0;
for j = x
    for k = y
        if rand > probability_initially_susceptible
            cell_value(j, k) = cell_value(j, k) + 1;
            initial_infected = initial_infected + 1;
        end
    end
end

cell_value_extend(2:rows + 1, 2:cols + 1) = cell_value;
observation_period{1} = cell_value;
observation_period_extend{1} = cell_value_extend;

% population arrays
susceptible_array = zeros(1, number_iterations);
susceptible_array(1) = total_population - initial_infected;
infected_array = zeros(1, number_iterations);
infected_array(1) = initial_infected;
recovered_array = zeros(1, number_iterations);

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
    if cell_value(x, y) == length_of_infection + length_of_immunity + 1 
        cell_value(x, y) = 0;
    end

    cell_value(x, y) = cell_value(x, y)+(rand > neighbor_sum(x, y, ...
                       cell_value_extend));
    cell_value_extend(2:rows + 1,2:cols + 1) = cell_value;
    observation_period{i} = cell_value;
    observation_period_extend{i} = cell_value_extend;
end

