function simulate_mqam_ber(constellation_size)
    % Main function to run M-QAM BER simulation
    % Input: constellation_size - Size of QAM constellation (4, 16, 64, etc.)
    
    % Initialize parameters
    params = initialize_parameters(constellation_size);
    
    % Generate constellation
    [constellation_levels, normalization_factor] = ...
        generate_constellation(params.constellation_size);
    
    % Setup Gray coding
    gray_to_binary_map = setup_gray_coding(params.bits_per_symbol);
    
    % Initialize result arrays
    number_of_bit_errors = zeros(1, length(params.eb_n0_db));
    
    % Main simulation loop
    for snr_index = 1:length(params.eb_n0_db)
        number_of_bit_errors(snr_index) = ...
            simulate_single_snr(params, snr_index, ...
                              constellation_levels, ...
                              normalization_factor, ...
                              gray_to_binary_map);
    end
    
    % Calculate and plot results
    plot_results(params, number_of_bit_errors);
end

function params = initialize_parameters(constellation_size)
    % Initialize simulation parameters
    params = struct();
    params.number_of_symbols = 10^5;
    params.constellation_size = constellation_size;
    params.bits_per_symbol = log2(params.constellation_size);
    
    % Verify valid constellation size
    if mod(log2(constellation_size), 1) ~= 0 || constellation_size < 4
        error('Constellation size must be a power of 2 and >= 4');
    end
    
    params.eb_n0_db = 0:15;
    params.es_n0_db = params.eb_n0_db + ...
        10 * log10(params.bits_per_symbol);
end

function [constellation_levels, normalization_factor] = ...
        generate_constellation(constellation_size)
    % Generate M-QAM constellation points
    sqrt_M = sqrt(constellation_size);
    base_levels = -(sqrt_M - 1);
    step_size = 2;
    end_level = (sqrt_M - 1);
    
    constellation_levels = base_levels:step_size:end_level;
    
    % Calculate normalization factor to maintain unit average symbol energy
    % For M-QAM with equally spaced points, average power is: 2(M-1)/3
    average_power = 2*(constellation_size-1)/3;
    normalization_factor = 1/sqrt(average_power);
end

function gray_to_binary_map = setup_gray_coding(bits_per_half_symbol)
    % Setup Gray coding mapping for each dimension
    binary_indices = 0:(2^(bits_per_half_symbol)-1);
    gray_map_indices = bitxor(binary_indices, floor(binary_indices/2));
    [~, gray_to_binary_map] = sort(gray_map_indices);
end

function [real_bits, imag_bits, real_decimal, imag_decimal] = ...
        generate_random_data(params)
    % Generate and process random binary data
    input_bits = rand(1, params.number_of_symbols * ...
                        params.bits_per_symbol) > 0.5;
    input_bits_matrix = reshape(input_bits, ...
                              params.bits_per_symbol, ...
                              params.number_of_symbols).';
    
    % Setup binary to decimal conversion
    power_vector = (params.bits_per_symbol/2 - 1):-1:0;
    weight_vector = 2.^power_vector;
    bin_to_dec_weights = ones(params.number_of_symbols, 1) * ...
                        weight_vector;
    
    % Process real and imaginary parts
    real_bits = input_bits_matrix(:, 1:params.bits_per_symbol/2);
    imag_bits = input_bits_matrix(:, ...
                params.bits_per_symbol/2 + 1:params.bits_per_symbol);
    
    real_decimal = sum(real_bits .* bin_to_dec_weights, 2);
    imag_decimal = sum(imag_bits .* bin_to_dec_weights, 2);
end

function symbols = map_to_symbols(real_decimal, imag_decimal, ...
                                constellation_levels, ...
                                normalization_factor)
    % Map decimals to constellation points
    real_gray_decimal = bitxor(real_decimal, ...
                              floor(real_decimal/2));
    imag_gray_decimal = bitxor(imag_decimal, ...
                              floor(imag_decimal/2));
    
    real_symbols = constellation_levels(real_gray_decimal + 1);
    imag_symbols = constellation_levels(imag_gray_decimal + 1);
    
    symbols = normalization_factor * ...
             (real_symbols + 1i * imag_symbols);
end

function received_signal = add_noise(symbols, params, snr_index)
    % Add Gaussian noise to symbols
    noise_real = randn(1, params.number_of_symbols);
    noise_imag = randn(1, params.number_of_symbols);
    noise = 1/sqrt(2) * (noise_real + 1i * noise_imag);
    
    attenuation_factor = 10^(-params.es_n0_db(snr_index)/20);
    received_signal = symbols + attenuation_factor * noise;
end

function [detected_real_bits, detected_imag_bits] = ...
        detect_symbols(received_signal, params, ...
                      constellation_levels, ...
                      normalization_factor, ...
                      gray_to_binary_map)
    % Detect symbols and convert back to bits
    received_real = real(received_signal)/normalization_factor;
    received_imag = imag(received_signal)/normalization_factor;
    
    % Symbol detection using minimum distance detection
    [~, real_indices] = min(abs(received_real - constellation_levels'), [], 1);
    [~, imag_indices] = min(abs(received_imag - constellation_levels'), [], 1);
    
    detected_real = constellation_levels(real_indices);
    detected_imag = constellation_levels(imag_indices);
    
    % Convert detected symbols to binary
    [detected_real_bits, detected_imag_bits] = ...
        symbols_to_bits(detected_real, detected_imag, ...
                       params, gray_to_binary_map);
end

function [detected_real_bits, detected_imag_bits] = ...
        symbols_to_bits(detected_real, detected_imag, ...
                       params, gray_to_binary_map)
    % Convert detected symbols back to bits
    sqrt_M = sqrt(params.constellation_size);
    scale_factor = (sqrt_M - 1);
    
    detected_real_decimal = ...
        gray_to_binary_map(floor((detected_real + scale_factor)/2 + 1)) - 1;
    detected_imag_decimal = ...
        gray_to_binary_map(floor((detected_imag + scale_factor)/2 + 1)) - 1;
    
    detected_real_bits = reshape_detected_bits(...
        detected_real_decimal, ...
        params.bits_per_symbol/2, ...
        params.number_of_symbols);
    detected_imag_bits = reshape_detected_bits(...
        detected_imag_decimal, ...
        params.bits_per_symbol/2, ...
        params.number_of_symbols);
end

function reshaped_bits = reshape_detected_bits(decimal_values, ...
                                             bits_per_half_symbol, ...
                                             number_of_symbols)
    % Reshape detected bits into proper format
    bits_str = dec2bin(decimal_values, bits_per_half_symbol);
    bits_str = bits_str.';
    bits_str = bits_str(1:end).';
    
    reshaped_bits = reshape(str2num(bits_str).', ...
                           bits_per_half_symbol, ...
                           number_of_symbols).';
end

function bit_errors = count_bit_errors(real_bits, imag_bits, ...
                                     detected_real_bits, ...
                                     detected_imag_bits)
    % Count total bit errors
    real_bit_errors = size(find(real_bits - detected_real_bits), 1);
    imag_bit_errors = size(find(imag_bits - detected_imag_bits), 1);
    bit_errors = real_bit_errors + imag_bit_errors;
end

function number_of_bit_errors = simulate_single_snr(params, ...
                                                  snr_index, ...
                                                  constellation_levels, ...
                                                  normalization_factor, ...
                                                  gray_to_binary_map)
    % Simulate for a single SNR value
    [real_bits, imag_bits, real_decimal, imag_decimal] = ...
        generate_random_data(params);
    
    symbols = map_to_symbols(real_decimal, imag_decimal, ...
                           constellation_levels, ...
                           normalization_factor);
    
    received_signal = add_noise(symbols, params, snr_index);
    
    [detected_real_bits, detected_imag_bits] = ...
        detect_symbols(received_signal, params, ...
                      constellation_levels, ...
                      normalization_factor, ...
                      gray_to_binary_map);
    
    number_of_bit_errors = count_bit_errors(...
        real_bits, imag_bits, ...
        detected_real_bits, detected_imag_bits);
end

function plot_results(params, number_of_bit_errors)
    % Calculate and plot BER results
    total_bits = params.number_of_symbols * params.bits_per_symbol;
    simulated_ber = number_of_bit_errors/total_bits;
    
    % Calculate theoretical BER for M-QAM
    M = params.constellation_size;
    k = log2(M);
    eb_n0_linear = 10.^(params.eb_n0_db/10);
    theoretical_ber = 4 * (1 - 1/sqrt(M)) * ...
        qfunc(sqrt(3*k*eb_n0_linear/(M-1))) / k;
    
    % Create plot
    figure;
    semilogy(params.eb_n0_db, theoretical_ber, 'bs-', 'LineWidth', 2);
    hold on;
    semilogy(params.eb_n0_db, simulated_ber, 'mx-', 'LineWidth', 2);
    
    % Configure plot
    axis([0 15 10^-5 1]);
    grid on;
    legend('Theoretical', 'Simulation');
    xlabel('Eb/N0, dB');
    ylabel('Bit Error Rate');
    title(['Bit Error Probability Curve for ' ...
           num2str(M) '-QAM Modulation']);
end