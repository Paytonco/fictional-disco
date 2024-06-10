% Load the .fig file
figFile = 'M_NoStress.fig';  % Specify the name of the .fig file to be loaded
open(figFile);  % Open the .fig file in MATLAB

% Get the current figure handle
figHandle = gcf;  % Get the handle of the current figure

% Find all surf objects in the figure
surfObjects = findall(figHandle, 'Type', 'surface');  % Find all surface objects in the figure

% Assuming there is only one surf object in the figure
if length(surfObjects) == 1  % Check if there is exactly one surface object
    surfObj = surfObjects;  % Assign the surface object to surfObj
else
    error('Multiple surf objects found. Please ensure there is only one surf plot in the figure.');  % Display an error message if there are multiple surface objects
end

% Get the Z data from the surf object
zData = get(surfObj, 'ZData');  % Retrieve the Z data (height values) from the surface object

% Flatten the Z data matrix to a vector
zVector = zData(:);  % Flatten the Z data matrix into a single column vector

% Calculate the maximum and median values
maxZ = max(zVector);  % Calculate the maximum value of the Z data
medianZ = median(zVector);  % Calculate the median value of the Z data

% Display the results
fprintf('Maximum Z value: %.4f\n', maxZ);  % Display the maximum Z value with four decimal places
fprintf('Median Z value: %.4f\n', medianZ);  % Display the median Z value with four decimal places
