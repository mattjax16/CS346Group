% Matt Bass, Jackson , Will
%
% mnbass21_jbpret22_wcjohn22_hw2_2_CS346.m.m
% CS346 -- Computational Modeling and Simulation
% Spring, 2021



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ ME:
 
% (ABOUT SCRIPT)
% This script is to run a Predator-Prey Model for where both the 
% predator and the prey are hunted by humans. (In this model the 
% predator will be sharks the prey will be tuna and they will both be 
% fished by humans). 
    
% This script will plot the populations of the shark and tuna 
% over the time rang of the simulation. The rate of change of these 
% two populations come from Lotka-Volterra Models differential
% equations:
%           (dY/dt) = K_Y(Y) - K_PY(YP) - K_YH(YH)
%    
%           (dP/dt) = (b)K_PY(YP) - K_P(P) - K_PH(PH)
    
    
%(HOW TO RUN):
% change the variables below to run the simulation and produce a graph
% plotting the tuna population (Y) and shark population (P) over the
% time of the simulation
    
% These variables the basic running of the finite differential equations
% simulation
numMonths = 12; %(int) This is the number of months to tun the sim
dt = 0.0001; %(float) Time step of the sim in units of months 
           % See write up for more
    
% Starting population variables that change
shark_population = 20; %(float, variable P in equation)   
tuna_population = 106;  %(float, variable Y in equation)
    
human_population = 20; %(float, variable H in equation)
fishing_rate = 0.01; %(float between 0 and 1, reresents the rate-constants 
                     % K_PH and K_YH)
	
    
tuna_birth_fraction = 2; %(K_Y)
tuna_death_proportionality_constant = 0.02; %(K_PY float on scale
                                                % on tenths)
	
shark_birth_fraction = 0.01; %((b)K_PY float smaller than K_PY)
shark_death_proportionality_constant = 1.06; %(K_P)


% boolean variables to determine which plots are shown
% if both plots are chosen they will be in a sub-plot
show_population_over_time_plot = true;
show_population_cycle_phase_plane = true;
    
% Change the variables above in the simulation for different results 
% Then run the script with the run button in MatLab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%Initilized the starting values fot rest of Sim

% Tuna 
tuna_births = tuna_birth_fraction*tuna_population; % K_Y(Y)
tuna_deaths = tuna_population*((tuna_death_proportionality_constant * ...
    shark_population) +(fishing_rate * ...
    human_population));% (K_PY(YP) + K_YH(YH)) (negetive pulled out)

% Shark
% (b)K_PY(YP)
shark_births = shark_birth_fraction * tuna_population * shark_population;

% (K_P(P) + K_PH(PH)) (negetive pulled out)
shark_deaths = shark_population*(shark_death_proportionality_constant + ...
    fishing_rate * human_population);


% Create number of iterations for the sim
numIterations = numMonths/dt;

t = 0; % starting time is 0

% initilize lists to hold populations to plot and the time at each step
tArray = ones(1,numIterations)*t;
sharkArray = ones(1,numIterations)*shark_population;
tunaArray = ones(1,numIterations)*tuna_population;

% simulation loop
for i = 2:numIterations
    t = i * dt;
    tuna_population = tuna_population + (tuna_births - tuna_deaths) * dt;
    shark_population = shark_population + (shark_births - shark_deaths)...
                        * dt;

    tuna_births = tuna_birth_fraction*tuna_population;
    tuna_deaths = tuna_population*((tuna_death_proportionality_constant * ...
                 shark_population) +(fishing_rate * human_population));

    shark_births = shark_birth_fraction * tuna_population * ...
        shark_population;

    shark_deaths = shark_population*(shark_death_proportionality_constant...
        + fishing_rate * human_population);
    
    %assign to arrays
    tArray(i) = t;
    sharkArray(i) = shark_population;
    tunaArray(i) = tuna_population;
end


% See if any graphs need to be created
if show_population_over_time_plot == true ||...
        show_population_cycle_phase_plane == true
    figure
end


%Create graph of populations over simulation
if show_population_over_time_plot == true
    %check what to see if subplot is needed
    if show_population_cycle_phase_plane == true
        subplot(1,2,1)
        sgtitle(["Predator Prey Model of Sharks hunting Tuna with both "+...
        "fished by Humans over " + {numMonths} + " Months","","dt= "+{dt} + ...
        "   H= " + {human_population} + "   K_H_P & K_Y_P= " + {fishing_rate}...
        + "   K_P= " + {shark_death_proportionality_constant} + "   K_Y= " + ...
        {tuna_birth_fraction} + "   K_P_Y= " + {tuna_death_proportionality_constant}...
        + "   (b)K_P_Y= " + {shark_birth_fraction}])
    
    end
    
    plot(tArray,sharkArray, 'r')
    hold
    plot(tArray,tunaArray, 'b')
    
    legend(["Sharks (predator P)", "Tuna (prey Y)" ],'FontSize',9)
    xlabel("Time (Months)")
    ylabel("# in Population")
    
    if show_population_cycle_phase_plane == false
        title(["Predator Prey Population Phase Cycle Model of Sharks" + ...
        "hunting Tuna with both "+...
        "fished by Humans", "H= " + {human_population} + ...
        "   K_H_P & K_Y_P= " + {fishing_rate} + "   K_P= " + ...
        {shark_death_proportionality_constant} + "   K_Y= " + ...
        {tuna_birth_fraction} + "   K_P_Y= " + ...
        {tuna_death_proportionality_constant}...
        + "   (b)K_P_Y= " + {shark_birth_fraction}])
    end
   
end


% Population Cycle Phase Plane plot

if show_population_cycle_phase_plane == true
    
    %see if subplot is need
    if show_population_over_time_plot == true
         subplot(1,2,2)
    end
    
    % Make the vector field for the graph
    % based off of:
    %           (dY/dt) = K_Y(Y) - K_PY(YP) - K_YH(YH)
    %    
    %           (dP/dt) = (b)K_PY(YP) - K_P(P) - K_PH(PH)
    
    max_pop = max(max(sharkArray,tunaArray))*1.05;
    [vect_tuna,vect_shark]=meshgrid(0:10:max_pop,...
        0:10:max_pop*1.05);

    U = tuna_birth_fraction*vect_tuna - (tuna_death_proportionality_constant... 
        .*vect_tuna.*vect_shark+fishing_rate.*human_population.*vect_tuna);

    V = shark_birth_fraction.*vect_tuna.*vect_shark-(vect_shark.*...
        (shark_death_proportionality_constant+fishing_rate*...
        human_population));

    L=sqrt(U.^2+V.^2);
    quiver(vect_tuna,vect_shark,U./L,V./L,.5,'k')
    hold on
    scatter(tunaArray,sharkArray)
    xlabel("Tuna(Y)")
    ylabel("Shark(P)")
    
    %see if title is needed 
    if show_population_over_time_plot == false
        title(["Predator Prey Model of Sharks hunting Tuna with both "+...
        "fished by Humans over " + {numMonths} + " Months","","dt= "+{dt} + ...
        "   H= " + {human_population} + "   K_H_P & K_Y_P= " + {fishing_rate}...
         + "   K_P= " + {shark_death_proportionality_constant} + "   K_Y= " + ...
        {tuna_birth_fraction} + "   K_P_Y= " + {tuna_death_proportionality_constant}...
        + "   (b)K_P_Y= " + {shark_birth_fraction}])
    end
    % Finding the Equilibrium point (see write up for equations)
    equilibrium_tuna_pop = (tuna_birth_fraction - (fishing_rate * ...
                human_population))/tuna_death_proportionality_constant;

    equilibrium_shark_pop = (shark_death_proportionality_constant + ...
        (fishing_rate * human_population))/shark_birth_fraction;



    scatter(equilibrium_shark_pop,equilibrium_tuna_pop, 'blue', 'Filled')
    legend(["Population Vector Field","Population Phaser Cycle initial (" + ...
        {tunaArray(1)} + " Y , " + {sharkArray(1)} + " P)",...
        "Non 0's Equilibrium Point (" + {equilibrium_shark_pop} + " Y , " + ...
        {equilibrium_tuna_pop} + " P)" ],'FontSize',9, 'Location' ,...
        'southeast')

    axis equal tight
end