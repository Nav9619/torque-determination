
clear;
clc;

% initial variables
theta_2_init = 85.31;  % initial input angle
theta_3_init = 270.13;
theta_4_init = 301.43;
omega_2 = 4  % given input angular velocity in degrees per sec 
r1 = 76.14;    % link lengths (which are fixed throughout)
r2 = 43.16;
r3 = 30;
r4 = 66.28;
r5 = 19.40;  % used for position and velocity analysis for CP on door 
r6 = 15;

v2 = r2 * omega_2; % input velocity;
F_out = 0.25*(30*70);  % force due to pressure calculated from given dynamic pressure and door area

n = 0   % iteration counter
theta_2 = theta_2_init;  % for iteration
theta_3 = theta_3_init;
theta_4 = theta_4_init;

omega_3 = 0;

% create matrix to hold kinematic coefficients
theta_3_prime = 1;  % assign random values to these initially
theta_4_prime = 1;
coeff_matrix = [theta_3_prime; theta_4_prime];

% create matrix for change in theta 3 and theta 4
delta_theta_3 = 0;
delta_theta_4 = 0;
delta_matrix = [delta_theta_3;delta_theta_4];

theta_2_array = [];
omega_3_array = [];  % to store the omega 3 values
inputTorquearray = [];


while (theta_2 < 121.31)
  
  n = n+1;
  theta_2_array = horzcat(theta_2_array,theta_2); % store current theta2 value in array
  
  
  % for finding change in theta 3 and theta 4
  f_a = r2*cosd(theta_2)- r3*cosd(theta_3)- r4*cosd(theta_4)-r1;
  f_b = r2*sind(theta_2)- r3*sind(theta_3)- r4*sind(theta_4);
  function_matrix1 = [f_a;f_b];
  
  % for finding kinematic coeffient  
  f1 = r2 * sind(theta_2);
  f2 = -r2 * cosd(theta_2);
  function_matrix2 = [f1;f2];
  
  % create Jacobian matrix
  J11 = r3*sind(theta_3);
  J12 = r4*sind(theta_4);
  J21 = -r3*cosd(theta_3);
  J22 = -r4*cosd(theta_4);
  [J]= [J11 J12;
       J21 J22];
  
  % find change in theta 3 and 4
  [delta_matrix] = -J^(-1) * [function_matrix1];
  delta_theta_3 = delta_matrix(1);
  delta_theta_4 = delta_matrix(2);
  
  % find kinematic coefficients
  [coeff_matrix] = -J^(-1) * [function_matrix2];
  theta_3_prime = coeff_matrix(1);
  theta_4_prime = coeff_matrix(2);
  
  
  omega_3 = theta_3_prime * omega_2; % find new value of omega 3 for this iteration
  
  omega_3_array = horzcat(omega_3_array,omega_3); % store the omega 3 value for this iteration
  
  % position analysis for CP on door
  r7x = (r2*cosd(theta_2)) - (r5*sind(theta_3)) + (r6*cosd(theta_3));
  r7y = (r2*sind(theta_2)) + (r5*cosd(theta_3)) + (r6*sind(theta_3));
  % r_out = sqrt((r7x)^2 +(r7y)^2);
  theta_out = atan(r7y/r7x);
  
  % velocity anallysis for CP on door
  v7x = (-r2*omega_2*sind(theta_2)) - (r5*omega_3*cosd(theta_3)) - (r6*omega_3*sind(theta_3));
  v7y = (r2*omega_2*cosd(theta_2)) - (r5*omega_3*sind(theta_3)) + (r6*omega_3*cosd(theta_3));
  v_out = sqrt((v7x)^2 + (v7y)^2);
  
  
  % Calculate input force and torque
  F_in = ((v_out*cosd(theta_out))/(v2*cosd(theta_2))) * F_out; 
  T_in = F_in * r2;
  
  inputTorquearray = horzcat(inputTorquearray,T_in); % store the input torque value for this iteration 
  
  % find new values of theta 3 and 4 for next iteration
  theta_3 = theta_3 + delta_theta_3;
  theta_4 = theta_4 + delta_theta_4;
  
  theta_2 = theta_2 + 1; % increment theta 2 by 1 degree for next iteration
  
end

% Plot figures 
plot(theta_2_array, omega_3_array);
plot(theta_2_array, inputTorquearray);

