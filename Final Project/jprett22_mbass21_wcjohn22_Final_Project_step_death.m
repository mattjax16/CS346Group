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
observation_period = {};
observation_period_extend = {};

% constants
length_of_infection = 2; % length in days that infection lasts
length_of_immunity = 5;
transmission_constant = .25; % chance a susceptible becomes infected ...
                             % when interacting with 1 infected
chance_of_death = .05;

                             
%cell states
susceptible_state = 0; 
dead_state = 9; 

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
     cell_value = observation_period{i - 1};
     cell_value_extend = observation_period_extend{i - 1};
     for row = 2:rows + 1
         for col = 2:cols + 1
             if cell_value(row - 1, col - 1) == 0
                cell_value(row - 1, col - 1) = ...
                    cell_value(row - 1, col - 1) + (rand < ...
                    probability_of_infection(row, col, cell_value_extend));
             elseif cell_value(row - 1, col - 1) > 0 && ...
                 cell_value(row - 1, col - 1) < 8
                cell_value(row - 1, col - 1) = ...
                    cell_value(row - 1, col - 1) + 1;
                if cell_value(row - 1, col - 1) > 0 && ...
                        cell_value(row - 1, col - 1)<3 ...
                        && rand<= chance_of_death
                    cell_value(row - 1, col - 1) = dead_state;
                end
                if cell_value(row - 1, col - 1) == 8
                    cell_value(row - 1, col - 1) = 0;  
                end
             end
         end 
     end
     cell_value_extend(2:rows + 1, 2:cols + 1) = cell_value;
     observation_period{i} = cell_value;
     observation_period_extend{i} = cell_value_extend;
end






%Plotting the CA for each time step

upper_cv_bound = 4;
lower_cv_bound = 0;

%set the color map
set(groot,'DefaultFigureColormap',jet(4));
figure;
for i=1:1:number_iterations
    
    ca_model = observation_period{i};
   
    %Get max and mic cell value for the entire CA simulation
    
    imagesc(ca_model);
    caxis([lower_cv_bound,upper_cv_bound]);
    title(sprintf('Day: %d',i));
    cbh = colorbar();
    cbh.Ticks = linspace(0, 4, 1); 
    hold;
    axis equal; axis tight; axis xy;
    
    %wait to go onto next image
    fprintf('Waiting for any key to be pressed\n');
    w = waitforbuttonpress;
end