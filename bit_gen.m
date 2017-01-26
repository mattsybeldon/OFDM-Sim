function [ bitstream ] = bit_gen( number_of_bits, logic_0_val, logic_1_val )
%BIT_GEN Generates bitstream with specified num_bits and logic values.
    bitstream = round(rand(1,number_of_bits));
    bitstream = bitstream .* ones(1,number_of_bits).* logic_1_val;
    bitstream(find(bitstream == 0)) = logic_0_val;
end

