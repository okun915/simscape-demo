function [E_cell,v_nernst,v_act,v_ohm,v_diff,heat_produced] = potential_cal(I,stack_T)

%potential_cal(I, cell_area, cell_num, T, alpha, ioref, ohm_effi, ilim)
%可变值
cell_area = 100; % [cm^2]单池面积
cell_num = 10; % [1]单池数量
alpha = 0.5; % Charge transfer coefficient
ioref = 8e-05; % [A/cm^2]Exchange current density
ohm_effi = 1; % Ohm overvoltage efficiency
ilim = 10; % [A/cm^2] limit Exchange current density


% 不变的常数
T_std = 298.15; % [K] Standard temperature
p_std = 101325; % [Pa] Standard pressure
std_state_entropy_change = -0.9*10^(-3); % [J/K/mol] standard state entropy change
G_H2O = -237140; % [J/mol] Gibbs free energy of water
E_O2 = 70000; % [J/mol] activation energy
F = 96485.33212; % [C/mol] Faraday constant
R = 8.31446261815324; % [J/K/mol] Universal gas constant
current_density = I/cell_area;

% 对ioref进行温度压力矫正
io = exp(-E_O2/R*(1/stack_T-1/T_std)) * ioref; 

% standard nernst voltage
E_std = G_H2O/(-2*F);
if ge(current_density, io)
    v_nernst = E_std + (stack_T-T_std)*std_state_entropy_change /(2*F); % + R*T/(2*F)* log((a_H2_C * a_O2_A^0.5) / a_H2O_A); %针对的是整个反应所以应该用电池的阴阳极活度来计算，而不应考虑催化层上的活度
else
    v_nernst = 0;
end

% Activation losses from Tafel equation
if ge(current_density, io)
    v_act = R * stack_T / (2 * alpha * F)*log(current_density/io);
else
    v_act = 0;
end

% Resistive voltage loss
%                         % Water content
%                         lambda_acl = membrane_water(RH_acl);
%                         lambda_ccl = membrane_water(RH_ccl);
%                     lambda_membrane = (lambda_acl + lambda_ccl)/2;
%                 % Membrane conductivity
%                 if ge(lambda_membrane, 1)
%                     sigma_30 = 0.005139*lambda_membrane - 0.00326; % [1/(Ohm*cm)]
%                 else 
%                     sigma_30 = 0.005139 - 0.00326; % [1/(Ohm*cm)]
%                 end
%             sigma = sigma_30 * exp(1268*(1/303.15 - 1/T));
%         % Membrane resistance
%         R_ohm = t_membrane / sigma;
    R_ohm = 0.3;
% Resistive voltage loss
v_ohm = R_ohm * current_density * ohm_effi;

% mass transform loss
v_diff = R * stack_T / (2 * alpha * F)*log(ilim/(ilim-current_density));

E_cell = v_nernst + v_act + v_ohm + v_diff;

% Electrical power consumed
heat_produced = cell_num * (E_cell-v_nernst) * current_density * cell_area; % [w]MEA产热速率

end