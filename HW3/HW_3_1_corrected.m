%Constants
capitance = 0.1; % microfarads per centimeter squared
applied_current = 0; %nanoampheres

sodium_displacement = 50; % displacement from the equilibrium  ...
                          % potential in milllivolts
potassium_displacement = -77; % displacement from the equilibrium ...
                              % potential in milllivolts
leakage_displacement = -54.4; % displacement from the equilibrium ...
                               % potential in milllivolts

sodium_conductance = 120; % maximum conductance in millisiemens per ...
                          % centimeters squared
potassium_conductance = 36; % maximum conductance in millisiemens per ...
                            % centimeters squared
leakage_conductance = 0.3; % maximum conductance in millisiemens per ...
                           % centimeters squared
depolorization = 0; 
hyperpolorization = 1;

time = 3; % milliseconds
delta_time = .001; % milliseconds
number_iterations = time/delta_time;
time_array = zeros(1, number_iterations); % time array

an = @(voltage) (0.01 * (voltage + 55))/(1 - exp(-(voltage + 55) / 10));
     % n-gate opening rate constant in meters per second
am = @(voltage) (0.1 * (voltage + 40))/(1 - exp(-(voltage + 40) / 10));
     % m-gate opening rate constant in meters per second
ah = @(voltage) (0.07 * exp(-(voltage + 65) / 20));
     % h-gate opening rate constant in meters per second
bn = @(voltage) (0.125 * exp(-(voltage + 65) / 80));
     % n-gate closing rate constant in meters per second
bm = @(voltage) (4 * exp(-(voltage + 65) / 18));
     % m-gate closing rate constant in meters per second
bh = @(voltage) (1 / (exp(-(voltage + 35) / 10) + 1));
     % h-gate closing rate constant in meters per second

sodium_channel_current = @(voltage, n, m, h) (sodium_conductance * m^3 * ...
    h * (voltage - sodium_displacement)); %nanoampheres
potassium_channel_current = @(voltage, n, m, h) (potassium_conductance * ...
    n^4 * (voltage - potassium_displacement)); %nanoampheres
leakage_current = @(voltage, n, m, h) (leakage_conductance * ...
	(voltage - leakage_displacement)); %nanoampheres

dvdt = @(voltage, I, n, m, h) (((1 / capitance) * ...
     (I - (sodium_channel_current(voltage, n, m, h) + ...
     leakage_current(voltage, n, m, h)))));
dvdt2 = @(voltage, I, n, m, h) (((1 / capitance) * ...
     (I - (potassium_channel_current(voltage, n, m, h) + ...
     leakage_current(voltage, n, m, h)))));
dvdt3 = @(voltage, I, n, m, h) (((1 / capitance) * ...
     (I - (leakage_current(voltage, n, m, h)))));
 

dndt = @(voltage, n, m, h) ((an(voltage) * (1 - n) - bn(voltage) * n));
dmdt = @(voltage, n, m, h) ((am(voltage) * (1 - m) - bm(voltage) * m));
dhdt = @(voltage, n, m, h) ((ah(voltage) * (1 - h) - bh(voltage) * h));
   
% initial values
V = zeros(1, number_iterations); % initialized voltage array
V(1) = -65; % initial voltage in millivolts
n = zeros(1, number_iterations); % initialized n array
n(1) = an(V(1)) / (an(V(1)) + bn(V(1))); % initial n-value
m = zeros(1, number_iterations); % initialized m array
m(1) = am(V(1)) / (am(V(1)) + bm(V(1))); % initial m-value
h = zeros(1, number_iterations); % initialized h array
h(1) = ah(V(1)) / (ah(V(1)) + bh(V(1))); % initial h-value
       
for i = 2:number_iterations
    time_array(i) = (i - 1) * delta_time;
    if i == 500;
         applied_current = 15;
     end
    if i == 1001
         applied_current = 0;
    end
    if V(i-1) < 49.3 && depolorization == 0

        k1 = delta_time * dvdt(V(i - 1), applied_current, n(i - 1), m(i - 1), h(i - 1));
        n1 = delta_time * dndt(V(i - 1), n(i - 1), m(i - 1), h(i - 1));
        m1 = delta_time * dmdt(V(i - 1), n(i - 1), m(i - 1), h(i - 1)); 
        h1 = delta_time * dhdt(V(i - 1), n(i - 1), m(i - 1), h(i - 1));

        k2 = delta_time * dvdt(V(i - 1) + (0.5 * k1), applied_current, n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        n2 = delta_time * dndt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        m2 = delta_time * dmdt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        h2 = delta_time * dhdt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));

        k3 = delta_time * dvdt(V(i - 1) + (0.5 * k2), applied_current, n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        n3 = delta_time * dndt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        m3 = delta_time * dmdt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        h3 = delta_time * dhdt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));

        k4 = delta_time * dvdt(V(i - 1) + k3, applied_current, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        n4 = delta_time * dndt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        m4 = delta_time * dmdt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        h4 = delta_time * dhdt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
    elseif V(i-1) > 49.3 && depolorization == 0 && hyperpolorization == 1
        
        depolorization = 1;
        k1 = delta_time * dvdt2(V(i - 1), applied_current, n(i - 1), m(i - 1), h(i - 1));
        n1 = delta_time * dndt(V(i - 1), n(i - 1), m(i - 1), h(i - 1));
        m1 = delta_time * dmdt(V(i - 1), n(i - 1), m(i - 1), h(i - 1)); 
        h1 = delta_time * dhdt(V(i - 1), n(i - 1), m(i - 1), h(i - 1));

        k2 = delta_time * dvdt2(V(i - 1) + (0.5 * k1), applied_current, n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        n2 = delta_time * dndt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        m2 = delta_time * dmdt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        h2 = delta_time * dhdt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));

        k3 = delta_time * dvdt2(V(i - 1) + (0.5 * k2), applied_current, n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        n3 = delta_time * dndt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        m3 = delta_time * dmdt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        h3 = delta_time * dhdt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));

        k4 = delta_time * dvdt2(V(i - 1) + k3, applied_current, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        n4 = delta_time * dndt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        m4 = delta_time * dmdt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        h4 = delta_time * dhdt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
    elseif V(i-1) < 49.3 && depolorization == 1 && hyperpolorization == 1 && V(i-1) >= -65
        
        k1 = delta_time * dvdt2(V(i - 1), applied_current, n(i - 1), m(i - 1), h(i - 1));
        n1 = delta_time * dndt(V(i - 1), n(i - 1), m(i - 1), h(i - 1));
        m1 = delta_time * dmdt(V(i - 1), n(i - 1), m(i - 1), h(i - 1)); 
        h1 = delta_time * dhdt(V(i - 1), n(i - 1), m(i - 1), h(i - 1));

        k2 = delta_time * dvdt2(V(i - 1) + (0.5 * k1), applied_current, n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        n2 = delta_time * dndt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        m2 = delta_time * dmdt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        h2 = delta_time * dhdt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));

        k3 = delta_time * dvdt2(V(i - 1) + (0.5 * k2), applied_current, n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        n3 = delta_time * dndt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        m3 = delta_time * dmdt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        h3 = delta_time * dhdt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));

        k4 = delta_time * dvdt2(V(i - 1) + k3, applied_current, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        n4 = delta_time * dndt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        m4 = delta_time * dmdt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        h4 = delta_time * dhdt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
    elseif V(i-1)<-65 && depolorization ==1 && hyperpolorization ==1 
        hyperpolorization = 0;
        k1 = delta_time * dvdt2(V(i - 1), applied_current, n(i - 1), m(i - 1), h(i - 1));
        n1 = delta_time * dndt(V(i - 1), n(i - 1), m(i - 1), h(i - 1));
        m1 = delta_time * dmdt(V(i - 1), n(i - 1), m(i - 1), h(i - 1)); 
        h1 = delta_time * dhdt(V(i - 1), n(i - 1), m(i - 1), h(i - 1));

        k2 = delta_time * dvdt2(V(i - 1) + (0.5 * k1), applied_current, n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        n2 = delta_time * dndt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        m2 = delta_time * dmdt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        h2 = delta_time * dhdt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));

        k3 = delta_time * dvdt2(V(i - 1) + (0.5 * k2), applied_current, n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        n3 = delta_time * dndt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        m3 = delta_time * dmdt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        h3 = delta_time * dhdt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));

        k4 = delta_time * dvdt2(V(i - 1) + k3, applied_current, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        n4 = delta_time * dndt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        m4 = delta_time * dmdt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        h4 = delta_time * dhdt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
    elseif V(i-1)<-65 && depolorization == 1 && hyperpolorization == 0
        
        
        k1 = delta_time * dvdt2(V(i - 1), applied_current, n(i - 1), m(i - 1), h(i - 1));
        n1 = delta_time * dndt(V(i - 1), n(i - 1), m(i - 1), h(i - 1));
        m1 = delta_time * dmdt(V(i - 1), n(i - 1), m(i - 1), h(i - 1)); 
        h1 = delta_time * dhdt(V(i - 1), n(i - 1), m(i - 1), h(i - 1));

        k2 = delta_time * dvdt2(V(i - 1) + (0.5 * k1), applied_current, n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        n2 = delta_time * dndt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        m2 = delta_time * dmdt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));
        h2 = delta_time * dhdt(V(i - 1) + (0.5 * k1), n(i - 1) + ...
            (0.5 * n1), m(i - 1) + (0.5 * m1), h(i - 1) + (0.5 * h1));

        k3 = delta_time * dvdt2(V(i - 1) + (0.5 * k2), applied_current, n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        n3 = delta_time * dndt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        m3 = delta_time * dmdt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));
        h3 = delta_time * dhdt(V(i - 1) + (0.5 * k2), n(i - 1) + ...
            (0.5 * n2), m(i - 1) + (0.5 * m2), h(i - 1) + (0.5 * h2));

        k4 = delta_time * dvdt2(V(i - 1) + k3, applied_current, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        n4 = delta_time * dndt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        m4 = delta_time * dmdt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
        h4 = delta_time * dhdt(V(i - 1) + k3, n(i - 1) + n3, m(i - 1) + m3, ...
            h(i - 1) + h3);
    end

    V(i) = V(i - 1) + 1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    n(i) = n(i - 1) + 1 / 6 * (n1 + 2 * n2 + 2 * n3 + n4);
    m(i) = m(i - 1) + 1 / 6 * (m1 + 2 * m2 + 2 * m3 + m4);
    h(i) = h(i - 1) + 1 / 6 * (h1 + 2 * h2 + 2 * h3 + h4);
end

figure()
plot(time_array, V)