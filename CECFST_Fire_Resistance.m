% CECFST heat transfer analysis + strength reduction calculation model

% input parameters
%Leff = 2.0; % exposed length, by Kodur (2011) and Lie TT (1990) et al.
Leff = 3.22; % pin-pin effective length，by Richard Liew et al.(2012)
L = 0.24; % concrete encasement dimension, in m
N = 50; % mesh number
dx = L / (N-1); % grid size
dt = 0.5; % time increment, in second
time_total = 150 * 60; % in second (150 min in total)

% Updated radius for concrete and steel ring parameters
radius_concrete = 0.065; % Radius of the concrete in meters
thickness_steel = 0.005; % Thickness of the steel ring in meters
radius_total = radius_concrete + thickness_steel; % Total radius including steel ring

% Mechanical Properties at ambient temperature
f_r = 523 * 10^6;  % Measured strength of steel tube in Pa
f_a = 360 * 10^6;  % Measured strength of steel tube in Pa
f_ci = 31.5 * 10^6;  % Measured strength of inner concrete in Pa
f_co = 36.6 * 10^6;  % Measured strength of outer concrete in Pa

E_ci = 22 * (f_ci*10^(-7))^0.3*10^9;  % Calculate E based on EC2
E_co = 22 * (f_co*10^(-7))^0.3*10^9;
E_s = 214*10^9; % test measured value for steel

% Thermal properties 
% heat conduction in W/m·K; lower limit
h = 25; % convective coefficient，W/m^2·K
epsilon = 0.7; % concrete emissivity
sigma = 5.67e-8; % Stefan-Boltzmann constant

% Set initial temp (Singapore average temperature)
T = 31 * ones(N); % initial temp Matrix in Celsius

% Initialise data (create array to store values over time）
TC_1=zeros(1, floor(time_total/dt) + 1);
TC_2=zeros(1, floor(time_total/dt) + 1);
TC_3=zeros(1, floor(time_total/dt) + 1);
TC_4=zeros(1, floor(time_total/dt) + 1);
TC_5=zeros(1, floor(time_total/dt) + 1);
TC_6=zeros(1, floor(time_total/dt) + 1);
TC_7=zeros(1, floor(time_total/dt) + 1);
Furnace_curve = zeros(1, floor(time_total/dt) + 1);  % furnace fire curve
N_OuterConcrete_sum = zeros(1, floor(time_total / dt) + 1);
N_rebar_sum = zeros(1, floor(time_total / dt) + 1);
EI_OuterConcrete = zeros(N, N);
EI_rebar_sum = zeros(1, floor(time_total / dt) + 1);
EI_OuterConcrete_sum = zeros(1, floor(time_total / dt) + 1);  % Array to store the sum of EI for every time step

% Time iterative
figure;
for t = 0:dt:time_total
    t_idx = floor(t/ dt) + 1;  % Define index for storing data (array index must be integer)
    % Furnace fire temperature function (test measured)
    if t/60 <= 20
        T_fire = 31 + 400 * log10(t/60 + 1);
    elseif t/60 <= 90
        T_fire = 560 + 8.9 * (t/60 - 20);
    elseif t/60 <= 150
        T_fire = 1183 + 3 * (t/60 - 90);
    end
    
    Furnace_curve(t_idx) = T_fire;  
    
    T_new = T; % Update the value of temp for successive time step
    for i = 1:N
        for j = 1:N
            x = (i-0.5) * dx; % x-coordinate of the center of the current cell
            y = (j-0.5) * dx; % y-coordinate of the center of the current cell
            if (i - N/2)^2 + (j - N/2)^2 >= (radius_total/L * N)^2 % Outer concrete 
                % Calculate thermal conductivity
                k_local = thermalConductivity(T(i, j));
                % Utility function to calculate specific heat
                c_p = calculateSpecificHeat(T(i, j));
                density = 2400; % in kg/m³
                if i == 1 || i == N || j == 1 || j == N  % Boundary points
                    A = dx^2/2; % Area of an external grid cell
                else
                    A = dx^2; % Area of an internal grid cell                
                end 
                k_T_local = calc_reduction_concrete(T(i, j));
                N_OuterConcrete = A * f_co * k_T_local; % Calculate resistance for the cell
                N_OuterConcrete_sum(t_idx) = N_OuterConcrete_sum(t_idx) + N_OuterConcrete;
                
                % Determine flexural stiffness (EI)
                I_xx = dx^4/12+dx^2*(y-L/2)^2;
                E_local = calculateConcreteModulusReduction(T(i, j))*E_co; 
                EI_OuterConcrete(i,j) = 0.8 * I_xx * E_local; % 0.8 for thermal stress reduction
                EI_OuterConcrete_sum(t_idx) = EI_OuterConcrete_sum(t_idx) + EI_OuterConcrete(i, j);
            elseif (i - N/2)^2 + (j - N/2)^2 <= (radius_concrete/L * N)^2 % Inner concrete
                % Calculate thermal conductivity
                k_local = thermalConductivity(T(i, j));
                % Utility function to calculate specific heat
                c_p = calculateSpecificHeat(T(i, j)); 
                density = 2400; % in kg/m³                
            elseif (i - N/2)^2 + (j - N/2)^2 <= (radius_total/L * N)^2 && (i - N/2)^2 + (j - N/2)^2 >= (radius_concrete/L * N)^2
                % Steel tube
                % Calculate thermal conductivity
                k_local = thermalConductivitySteel(T(i, j));
                % Utility function to calculate specific heat
                c_p = calculateSpecificHeatSteel(T(i, j)); 
                density = 7850; % in kg/m³                
            end         

            % Identify neighbors
            neighbors = 0;
            sum_T_neighbors = 0;
            if i > 1
                sum_T_neighbors = sum_T_neighbors + T(i-1, j);
                neighbors = neighbors + 1;
            end
            if i < N
                sum_T_neighbors = sum_T_neighbors + T(i+1, j);
                neighbors = neighbors + 1;
            end
            if j > 1
                sum_T_neighbors = sum_T_neighbors + T(i, j-1);
                neighbors = neighbors + 1;
            end
            if j < N
                sum_T_neighbors = sum_T_neighbors + T(i, j+1);
                neighbors = neighbors + 1;
            end

            % Conduction for all points
            conduction = k_local * (sum_T_neighbors - neighbors * T(i, j)) / dx^2;

            % Apply additional calculations for boundary points
            if i == 1 || i == N || j == 1 || j == N  % Boundary points
                % convection
                convection = h * (T_fire - T(i, j)) / dx;
                
                % radiation
                radiation = epsilon * sigma * ((T_fire + 273.15)^4 - (T(i, j) + 273.15)^4) / dx;
                
                % Update external points' temperature
                T_new(i, j) = T(i, j) + dt * (conduction + convection + radiation) / (density * c_p); 
            else  % Internal points
                % Update internal points' temperature
                T_new(i, j) = T(i, j) + dt * conduction / (density * c_p); 
            end
        end
    end
    T = T_new; % Update T_i+1

    % Record thermocouple temp
    TC_6(t_idx) = T(1, round(N/2)); % concrete surface 
    TC_3(t_idx) = T(round(N/2), round(N*(0.5*L+radius_total)/L)); % steel tube 
    TC_4(t_idx) = T(round(N/2), round(N*206/240)); % stirrup 
    TC_5(t_idx) = T(round(N*197/240), round(N*197/240)); % rebar 
    TC_7(t_idx) = T(round(N*(L/2+radius_total/sqrt(2))/L), round(N*(L/2+radius_total/sqrt(2))/L)); % steel tube upper right 
 
    temp_min =20; % minimum temperature for colormap
    temp_max = 1500; % max temperature for colormap

    % Calculate for rebar
    A_rebar = 0.000314; % Cross-sectional area of rebar in square meters
    I_rebar = 1.63 * 10^(-6); % Inertia of rebar in meters^4
    k_f_rebar = calc_reduction_rebar(TC_5(t_idx)); % Calculate strength reduction factor for rebar at current temperature
    N_rebar = A_rebar * f_r * k_f_rebar; % Calculate resistance for rebar
    N_rebar_sum(t_idx) = N_rebar_sum(t_idx) + N_rebar;
    k_E_rebar = calculateSteelModulusReduction(TC_5(t_idx)); % Calculate modulus reduction factor for rebar at current temperature
    EI_rebar = 0.8 * I_rebar * E_s * k_E_rebar; % Calculate flexural resistance for rebar, 0.8 for thermal stress reduction
    EI_rebar_sum(t_idx) = EI_rebar_sum(t_idx) + EI_rebar;

    % plot Cartesian step temp distribution
    if mod(t, 600) == 0
        imagesc(T);
        colormap(turbo);
        colorbar;
        clim([temp_min temp_max]); % Set consistent colormap range
        title(sprintf('Temperature Distribution at t = %d minutes', t / 60));
        drawnow;
    end
end

% plot T-t curve
figure;
hold on;
TC_3 = TC_3.';  % Transpose TC_3 from row to column
TC_4 = TC_4.';  % Transpose TC_4 from row to column
TC_5 = TC_5.';  % Transpose TC_5 from row to column
TC_6 = TC_6.';  % Transpose TC_6 from row to column
TC_7 = TC_7.';  % Transpose TC_7 from row to column
plot(0:dt:time_total, TC_3, 'g-', 'LineWidth', 2);
plot(0:dt:time_total, TC_4, 'c-', 'LineWidth', 2);
plot(0:dt:time_total, TC_5, 'm-', 'LineWidth', 2);
plot(0:dt:time_total, TC_6, 'y-', 'LineWidth', 2);
plot(0:dt:time_total, Furnace_curve, 'k-.', 'LineWidth', 2);  % furnace fire 
xlabel('Time (seconds)');
ylabel('Temperature (°C)');
title('Temperature-Time Curves Comparison');
legend('TC3','TC4','TC5','TC6', 'Fire Curve');
grid on;
hold off;

% Adjusting to polar coordinate system (assume 1-layer thin-walled steel tube)
N_concrete = 23; % Adjusted number of divisions for concrete
N_steel = 2; % Number of divisions for steel
N_circle = N_concrete + N_steel; % Total divisions
dr = radius_total / N_circle; % Radial step size

% Initialise 
T_circle = 31 * ones(N_circle, 1); % Initial temperature array/Celsius
N_CFST_m = zeros(1, N_circle); % This will hold the resistance for each ring element m
EI_CFST_m = zeros(1, N_circle); % This will hold the flexural resistance for each ring element m
N_CFST_sum = 0; % This will store the cumulative sum of resistances
EI_CFST_sum = 0; % This will store the cumulative sum of resistances

% data for plotting
A_m_array = zeros(1, N_circle);
N_CFST_sum_array = zeros(1, floor(time_total/dt) + 1);
EI_CFST_sum_array = zeros(1, floor(time_total/dt) + 1);

% Time iteration
figure;
for t = 0:dt:time_total
    t_idx = floor(t / dt) + 1;  % Index for storing data (array index must be integer)
    N_CFST_sum = 0; % Reset the sum at the start of each time step
    EI_CFST_sum = 0; 
    T_new_circle = T_circle; % Make a copy to update temperatures
    for m = 1:N_circle
        r_m = m * dr; % Radial position of current node
        if m <= N_concrete
            % Concrete properties calculations
            k_local = thermalConductivity(T_circle(m));
            c_p_local = calculateSpecificHeat(T_circle(m));
            density = 2400; 
            f_local = f_ci;
            k_c_f = calc_reduction_concrete(T_circle(m));
            % Update CFST strength for concrete only
            A_m = area_m(dr, m);
            A_m_array(m) = A_m;  % Store the area in the array (unnecessary)
            N_CFST_m(m) = k_c_f * A_m * f_local;
            N_CFST_sum = N_CFST_sum + N_CFST_m(m);  % Sum only concrete resistances
            % Update CFST EI for concrete only
            I_m = inertia_m(dr, m);
            k_c_E = calculateConcreteModulusReduction(T_circle(m));
            EI_CFST_m(m) = 0.8 * k_c_E * I_m * E_ci; % thermal stress effect 0.8 as per EC4 Annex G and Zhou & Han
            EI_CFST_sum = EI_CFST_sum + EI_CFST_m(m);  % Sum only concrete flexural resistance
        else
            % Steel properties 
            k_local = thermalConductivitySteel(T_circle(m));
            c_p_local = calculateSpecificHeatSteel(T_circle(m));
            density = 7850;
            f_local = f_a;
            k_a_f = calc_reduction_tube(T_circle(N_circle));
            k_a_E = calculateSteelModulusReduction(T_circle(N_circle));
        end
        
        % Initialize conduction term
        conduction = 0;

        if m == 1 % Center node
            % Only outward radial direction, and approximation as forward difference
            conduction = k_local * (T_circle(2) - T_circle(1)) / (dr^2);
        elseif m <= N_concrete % Concrete intermediate nodes
            % Standard finite difference approximation for second derivative in radial direction
            conduction = k_local * (T_circle(m-1) - 2 * T_circle(m) + T_circle(m+1)) / (dr^2);
            % Adding the 1/r dT/dr term using central difference
            radial_derivative = (T_circle(m+1) - T_circle(m-1)) / (2 * dr);
            conduction = conduction + k_local * radial_derivative / r_m;        
        else % Outermost node, using average of TC_3 and TC_7
            T_new_circle(m) = (TC_3(t/dt + 1) + TC_7(t/dt + 1)) / 2;
            continue;
        end

        % Update temperatures based on conduction only
        T_new_circle(m) = T_circle(m) + (dt * conduction) / (density * c_p_local);
    end

    %Update CFST strength with steel tube
    A_steel = pi * (radius_total^2 - radius_concrete^2);
    N_tube = A_steel * f_local * k_a_f;
    N_CFST_sum = N_CFST_sum + N_tube; 
    N_CFST_sum_array(t_idx) = N_CFST_sum; % Store the sum of resistances at this time step
    I_steel = pi * (radius_total^4 - radius_concrete^4)/4;
    EI_tube = I_steel * E_s * k_a_E;
    EI_CFST_sum = EI_CFST_sum + EI_tube; 
    EI_CFST_sum_array(t_idx) = EI_CFST_sum; % Store the sum of resistances at this time step
    T_circle = T_new_circle; % Update temperature distribution

    % Record center temperature for plotting
    TC_1(t/dt + 1) = T_circle(2);
    TC_2(t/dt + 1) = T_circle(floor(N_circle/2));
    % Record CFST strength for plotting

    % Plot the temperature distribution every 10 minutes
    if mod(t, 600) == 0
        theta = linspace(0, 2*pi, 100); % Angular coordinates
        r_matrix = linspace(0, radius_total, length(T_circle)); % Radial coordinates from center to edge
        r = repmat(r_matrix, length(theta), 1); % Extend radial coordinates across the angular dimension
        [X, Y] = pol2cart(theta', r);
    
        % Use surf to plot the temperature distribution across these coordinates
        surf(X, Y, repmat(T_circle', length(theta), 1), 'EdgeColor', 'none'); % Use temperature data as color
        title(sprintf('Isothermal Lines at t = %d minutes', t / 60));
        %clim([temp_min temp_max]);
        xlabel('X Coordinate');
        ylabel('Y Coordinate');
        zlabel('Temperature (°C)');
        view(2); % View the plot from above
        colormap(turbo);
        colorbar;
        axis equal; % Ensure the axes are scaled equally to maintain circular shape
    end
end

N_CECFST_array = N_CFST_sum_array + N_OuterConcrete_sum + N_rebar_sum;
EI_CECFST_array = EI_CFST_sum_array + EI_OuterConcrete_sum + EI_rebar_sum;
N_cr_array = pi^2 * EI_CECFST_array / Leff^2;
lambda = sqrt(N_CECFST_array ./ N_cr_array); % slenderness ratio as per Eurocode
Phi = 0.5 * (1 +  0.49 * (lambda - 0.2) + lambda.^(2)); % clause 6.1.3.2 as per EC3
chi = 1 ./ (Phi + sqrt(Phi.^2 - lambda.^2));
N_fi_Rd = chi .* N_CECFST_array; % Design fire resistance 

% Plotting N_CECFST-Time Curve
figure;
plot(0:dt:time_total, N_CECFST_array, 'b-', 'LineWidth', 2); % Adding markers
title('N_CECFST_array Plot');
xlabel('Time (seconds)');
ylabel('Normal Force (N)');
grid on;

figure;
plot(0:dt:time_total, N_fi_Rd, 'g-', 'LineWidth', 2); % Adding markers
title('N_fi_Rd Plot');
xlabel('Time (seconds)');
ylabel('Normal Force (N)');
grid on;

% Plotting Temperature-Time Curves
figure;
hold on;
plot(0:dt:time_total, TC_1, 'r-', 'LineWidth', 2); % From circle model
plot(0:dt:time_total, TC_2, 'g-', 'LineWidth', 2); % From circle model
xlabel('Time (seconds)');
ylabel('Temperature (°C)');
title('Temperature-Time Curves Comparison');
legend('Column centre Temp','Tube surface Temp');
grid on;
hold off;

TC_1 = TC_1.';  % Transpose TC_1 from row to column matrix
TC_2 = TC_2.';  % Transpose TC_2 from row to column matrix


% Defined functions 
function k_c = thermalConductivity(T)
    k_c = 1.36 - 0.136 * (T / 100) + 0.0057 * (T / 100).^2;
end

%assume 5% moisture content
%{
function c_p = calculateSpecificHeat(T)
    if T < 100
        c_p = 900;
    elseif T < 115
        c_p = 2600;
    elseif T < 200
        c_p = 2600 - 18.82 * (T - 115);
    elseif T < 400
        c_p = 1000 + (T - 200) / 2;
    else
        c_p = 1100;
    end
end
%}

%assume 10% moisture content
function c_p = calculateSpecificHeat(T)
    if T < 100
        c_p = 900;
    elseif T < 115
        c_p = 5600;
    elseif T < 200
        c_p = 5600 -(5600-1000)/(200-115) * (T - 115);
    elseif T < 400
        c_p = 1000 + (T - 200) / 2;
    else
        c_p = 1100;
    end
end

function c_a = calculateSpecificHeatSteel(T)
    if T <= 600
        c_a = 425 + 7.73*10^-1*T - 1.69*10^-3*T^2 + 2.22*10^-6*T^3;
    elseif T > 600 && T <= 735
        c_a = 666 - 13002 / (T - 738);
    elseif T > 735 && T <= 900
        c_a = 545 + 17820 / (T - 731);
    else
        c_a = 650; % Constant specific heat for temperatures above 900°C
    end
end

function k_a = thermalConductivitySteel(T)
    if T <= 800
        k_a = 54 - 3.33*10^-2*T;
    else
        k_a = 27.3; 
    end
end

function kc_T = calc_reduction_concrete(T)
    % Define the corresponding strength reduction factor for concrete
    temperatures = [20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1500];
    reduction_factors = [1.00, 1.00, 0.95, 0.85, 0.75, 0.6, 0.45, 0.3, 0.15, 0.08, 0.04, 0.01, 0.00];
    % Interpolate the strength reduction factor at the given temperature
    kc_T = interp1(temperatures, reduction_factors, T, 'linear');
end

function kr_T = calc_reduction_rebar(T)
    % Define the temperature and corresponding strength reduction factor for rebar 
    temperatures = [20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1500];
    reduction_factors = [1.00, 0.96, 0.92, 0.81, 0.63, 0.44, 0.26, 0.08, 0.06, 0.05, 0.03, 0.02, 0.00];
    % Interpolate the strength reduction factor at the given temperature
    kr_T = interp1(temperatures, reduction_factors, T, 'linear');
end

function ka_T = calc_reduction_tube(T)
    % Define the temperature and corresponding strength reduction factor for rebar 
    temperatures = [20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1500];
    reduction_factors = [1.00, 1.00, 1.00, 1.00, 1.00, 0.78, 0.47, 0.23, 0.11, 0.06, 0.04, 0.02, 0.00];
    % Interpolate the strength reduction factor at the given temperature
    ka_T = interp1(temperatures, reduction_factors, T, 'linear');
end

% Define the area of each circular ring element
function A = area_m(dr, m)
    r_inner = (m-1) * dr;
    r_outer = m * dr;
    A = pi * (r_outer^2 - r_inner^2);
end

% Define the inertia of each circular ring element
function I = inertia_m(dr, m)
    r_inner = (m-1) * dr;
    r_outer = m * dr;
    I = pi * (r_outer^4 - r_inner^4)/4;
end

% Define concrete modulus as function of temperature (DA Krishna)
function k_c_E = calculateConcreteModulusReduction(T)
        k_c_E = max(0, -0.001282*T + 1.0265);
end

% Define steel modulus as function of temperature (EC2-1-2)
function k_s_E = calculateSteelModulusReduction(T)
    temperatures = [20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1500];
    reduction_factors = [1.00, 1.00, 0.87, 0.72, 0.56, 0.40, 0.24, 0.08, 0.06, 0.05, 0.03, 0.02, 0.00];
    % Interpolate the strength reduction factor at the given temperature
    k_s_E = interp1(temperatures, reduction_factors, T, 'linear');
end
