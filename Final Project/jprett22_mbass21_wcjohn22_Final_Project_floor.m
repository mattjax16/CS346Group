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
%cell states
%the different states being used based on their cell value. 
susceptible_state = 0;
recovered_state = 3;
infected_state = 1;

% neighbor function for use in infection simulation
% input: coordinates of grid extended
% output: odds of infection given number of infected around susceptible
probability_of_infection = @(x, y, cell_value_ex) ...
    ((transmission_constant) ...
    * ((cell_value_ex(x + 1, y) > susceptible_state ...
    && cell_value_ex(x + 1, y) < ...
    recovered_state) ...
    + (cell_value_ex(x - 1, y) > susceptible_state ...
    && cell_value_ex(x - 1, y) < ...
    recovered_state) ...
    + (cell_value_ex(x, y + 1) > susceptible_state ...
    && cell_value_ex(x, y + 1) < ...
    recovered_state) ...
    + (cell_value_ex(x, y - 1) > susceptible_state ...
    && cell_value_ex(x ,y - 1) < ...
    recovered_state)));

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

for i = 2:number_iterations    
     cell_value = observation_period{i - 1};
     cell_value_extend = observation_period_extend{i - 1};
     for row = 2:rows + 1
         for col = 2:cols + 1
             if cell_value(row - 1, col - 1) == susceptible_state
                cell_value(row - 1, col - 1) = ...
                    cell_value(row - 1, col - 1) + (rand < ...
                    probability_of_infection(row, col, cell_value_extend));
             elseif cell_value(row - 1, col - 1) > susceptible_state && ...
                 cell_value(row - 1, col - 1) < recovered_state+ ...
                 length_of_immunity
                cell_value(row - 1, col - 1) = ...
                    cell_value(row - 1, col - 1) + 1;
                if cell_value(row - 1, col - 1) == recovered_state+...
                        length_of_immunity
                    cell_value(row - 1, col - 1) = susceptible_state;  
                end
             end
             %the following conditionals check values to implament graphs
             %for population tracking
             if cell_value(row - 1, col - 1) == susceptible_state
                 susceptible_array(i) = susceptible_array(i) + 1;
             elseif cell_value(row - 1, col - 1) > susceptible_state && ...
                     cell_value(row - 1, col - 1) < recovered_state
                 infected_array(i) = infected_array(i)+1;
             end
         end 
     end
     cell_value_extend(2:rows + 1, 2:cols + 1) = cell_value;
     observation_period{i} = cell_value;
     observation_period_extend{i} = cell_value_extend;
end






%Plotting the CA for each time step

upper_cv_bound = 8;
lower_cv_bound = 0;

%set the color map
set(groot,'DefaultFigureColormap',jet(64));
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
t = 1:number_iterations;
figure();
plot(t,infected_array, t, susceptible_array);
title("Number of infected and susceptable over the simulation");
xlabel("Days");
ylabel("Number of Infected People");
legend("Infected", "Susceptible");