function [ pulse_vector, time_vector ] = root_rcro( kt,Tb, samples_per_bit, r)
%ROOT_RCRO Summary of this function goes here
%   Root raised cosine rolloff pulse for |t| <= k_t x T_b
    stop_time = -kt*Tb;
    start_time= kt*Tb;
    time_vector = linspace(stop_time, start_time, samples_per_bit*2*kt);
    pulse_vector = zeros(1,length(time_vector));
    R = 1/Tb;
    % Case 1
    pulse_vector = (sin(pi*R.*time_vector*(1-r)) + ...
        4*R*r.*time_vector.*cos(pi*R.*time_vector.*(1+r))) ./ ...
        (pi*R.*time_vector.*(1-(4*R*r.*time_vector).^2));
    
    % Case 2
    pulse_vector(find(abs(time_vector) == (Tb/(4*r)))) =  ...
        (r/(sqrt(2))) * ((1+2/pi)*sin(pi/(4*r)) + (1-(2/pi))*cos(pi/(4*r)));
    % Case 3
    pulse_vector(find(time_vector == 0)) = (1 - r + (4*r/pi));
end

