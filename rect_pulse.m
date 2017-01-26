function [pulse_vector,time_vector ] = rect_pulse( bit_period, ...
                                            samples_per_bit)
%RECT_PULSE Generate rectangular pulse vector for some num_samples and
%           period of time
    total_time = bit_period * 10;
    stop_point = total_time / 2;
	time_vector = linspace(-stop_point, stop_point, 10*samples_per_bit);
    vector_length = length(time_vector);
    pulse_vector = zeros(1, vector_length);
    pulse_vector(4.5*samples_per_bit:5.5*samples_per_bit - 1) = ones(1,samples_per_bit);

end

