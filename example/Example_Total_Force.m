close all;
clear all;

files = dir('*.mat');
load(files(1).name);

timeStamps = zeros(1, length(data));
x_force=zeros(length(files), length(data));
y_force=zeros(length(files), length(data));
z_force=zeros(length(files), length(data));

for j = 1:length(files)
    % Load the .mat file
    load(files(j).name);

    % Store the results of the different time steps in arrays
    for i = 1:length(data)
    
        timeStamps(j,i) = data{1,i}.time;
        
        x_force(j,i)=data{1,i}.force(1);
        y_force(j,i)=data{1,i}.force(2);
        z_force(j,i)=data{1,i}.force(3);
    end

end

% Calculates the average friction
Force_avg = [];
timeLimit = 2.5;
for j = 1:length(files)
    Force_avg = [Force_avg sum(x_force(j,timeStamps(j,:) > timeLimit))...
        / length(x_force(j,timeStamps(j,:) > timeLimit)) ...
        + 0.09 * sum(z_force(j,timeStamps(j,:) > timeLimit))...
        / length(z_force(j,timeStamps(j,:) > timeLimit)) - 0.34];
end

Force_avg