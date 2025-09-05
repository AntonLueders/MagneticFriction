close all;
clear all;

files = dir('*.mat');
load(files(1).name);

timeStamps = zeros(1, length(data));
order_parameter = zeros(length(files), length(data));

for j = 1:length(files)
    % Load the .mat file
    load(files(j).name);

    % Store the results of the different time steps in arrays
    for i = 1:length(data)
    
        timeStamps(j,i) = data{1,i}.time;
        array = data{1,i}.array;

        % Calculates the order parameter
        temp_order_parameter = 0;
        for row = 2:7
            for col=1:size(array, 2)
                 temp_order_parameter = temp_order_parameter ...
                     + cos(array(row,col) - array(row-1,col));
            end 
        end
        order_parameter(j,i) = temp_order_parameter ...
            / (size(array, 1)*(size(array, 2) - 1));
        
    end

end

% Calculates the average order parameter
order_parameter_avg = [];
timeLimit = 2.5;
for j = 1:length(files)
    order_parameter_avg = [order_parameter_avg sum(order_parameter(j,timeStamps(j,:) > timeLimit))...
        / length(order_parameter(j,timeStamps(j,:) > timeLimit))];
end

order_parameter_avg