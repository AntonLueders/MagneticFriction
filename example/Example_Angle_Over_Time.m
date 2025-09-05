close all;
clear all;

files = dir('*.mat');
load(files(1).name);

timeStamps = zeros(1, length(data));
angle = zeros(length(files), length(data));

for j = 1:length(files)
    % Load the .mat file
    load(files(j).name);

    % Store the results of the different time steps in arrays
    for i = 1:length(data)
    
        timeStamps(j,i) = data{1,i}.time;
        array = data{1,i}.array;
        angle(j,i) = data{1,i}.array(4,4);

    end

end

% Plots the figures
for j = 1:length(files)
    figure(j)
    plot(timeStamps(j,:), angle(j,:));
end