function [ pulse_vector, time_vector ] = rcro_pulse( k_t, bit_period, samples_per_bit, r )
%rcro_pulse Generate Raised Cosine Rolloff pulse vector
    f0 = 1/(bit_period*2);
    f_delta = r*f0;
    
    stop_time = -k_t*bit_period;
    start_time= k_t*bit_period;
    time_vector = linspace(stop_time, start_time, samples_per_bit*2*k_t);
    
    pulse_vector_1 = ((2*f0*sin(2*pi*f0.*time_vector)  ./ ...
                        (2*pi*f0.*time_vector)));
    % Separate the second half of the equation to handle small roundoff
    % errors and 0/0
    pulse_vector_2 = cos(2*pi*f_delta*time_vector);
    pulse_vector_3 = (1-(4*f_delta.*time_vector).^2);
    pulse_vector_2(abs(pulse_vector_2) < 5e-14) = 0;
    pulse_vector_3(abs(pulse_vector_3) < 5e-14) = 0;
    pulse_vector_2 = pulse_vector_2 ./ pulse_vector_3;
    % Apply L'Hopital's rule
    pulse_vector_2(find(isnan(pulse_vector_2))) = ...
                        sin(2*pi*f_delta*time_vector(find(isnan(pulse_vector_2)))) ...
                        * pi ./ ...
                        (16*f_delta*time_vector(find(isnan(pulse_vector_2))));
    % Replace potential NaN with one, mainly at t=0.
    pulse_vector = pulse_vector_1 .* pulse_vector_2;
    pulse_vector(find(isnan(pulse_vector))) = 2*f0;
end

