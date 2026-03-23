function data_matrix1 = load_data(filename)

    data = load(filename);
    variables = fieldnames(data);
    signals = struct();

    for i = 1:length(variables)
        if contains(variables{i}, 'cf_ad9361_A_voltage')
            signals.(variables{i}) = data.(variables{i});
        end
    end

    signal_names = fieldnames(signals);
    num_signals = length(signal_names);

    if num_signals > 0
        len = length(signals.(signal_names{1}));
        data_matrix = zeros(len, num_signals);
        for i = 1:num_signals
            data_matrix(:, i) = signals.(signal_names{i});
        end
    else
        error('未找到符合条件的信号数据。');
    end

    if size(data_matrix, 2) < 8
        error('信号数据不足8列，无法构造4个复数通道。');
    end

    channel1 = complex(data_matrix(:,1), data_matrix(:,2));
    channel2 = complex(data_matrix(:,3), data_matrix(:,4));
    channel3 = complex(data_matrix(:,5), data_matrix(:,6));
    channel4 = complex(data_matrix(:,7), data_matrix(:,8));
    data_matrix1 = [channel1, channel2, channel3, channel4];
    
end
