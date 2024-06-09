% Load the .fig file
figFile = 'M_NoStress.fig';
open(figFile);

% Get the current figure handle
figHandle = gcf;

% Find all surf objects in the figure
surfObjects = findall(figHandle, 'Type', 'surface');

% Assuming there is only one surf object in the figure
if length(surfObjects) == 1
    surfObj = surfObjects;
else
    error('Multiple surf objects found. Please ensure there is only one surf plot in the figure.');
end

% Get the Z data from the surf object
zData = get(surfObj, 'ZData');

% Flatten the Z data matrix to a vector
zVector = zData(:);

% Calculate the maximum and median values
maxZ = max(zVector);
medianZ = median(zVector);

% Display the results
fprintf('Maximum Z value: %.4f\n', maxZ);
fprintf('Median Z value: %.4f\n', medianZ);
