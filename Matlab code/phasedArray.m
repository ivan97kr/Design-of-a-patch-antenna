
% MATLAB Code from Sensor Array Analyzer App

% Generated by MATLAB 9.11 and Phased Array System Toolbox 4.6

% Generated on 14-Jan-2022 16:17:52

% Create a Uniform Linear Array Object
Array = phased.ULA('NumElements',5,...
'ArrayAxis','x');
% The multiplication factor for lambda units to meter conversion
Array.ElementSpacing = 0.658*0.111;
sll = 41.58;
Array.Taper = chebwin(5, sll);

% Create an isotropic antenna element
 %Elem = phased.IsotropicAntennaElement;
 %Elem.FrequencyRange = [0 2.7e9];
 Elem = MyPCB1_40db;
 %Elem = pf_final;
Array.Element = Elem;

% Assign Frequencies and Propagation Speed
Frequency = 2.7e9;
PropagationSpeed = physconst('Lightspeed');

% Assign Steering Angles
SteeringAngles = [45;0];
% Assign Phase shift quantization bits
PhaseShiftBits = 0;

% Create Figure


% figure;
% pattern(MyPCB1_40db,Frequency)


% Calculate Steering Weights

Freq3D = 2.7e9;
% Find the weights
w = zeros(getNumElements(Array), length(Frequency));
SteerVector = phased.SteeringVector('SensorArray', Array,...
 'PropagationSpeed', PropagationSpeed, 'NumPhaseShifterBits', PhaseShiftBits(1));
for idx = 1:length(Frequency)
    w(:, idx) = step(SteerVector, Frequency(idx), SteeringAngles(:, idx));
end

% Plot Array Geometry
% figure;
% viewArray(Array,'ShowNormal',false,...
%   'ShowTaper',false,'ShowIndex','None',...
%   'ShowLocalCoordinates',true,'ShowAnnotation',false,...
%   'Orientation',[0;0;0]);


% Plot 3d graph
format = 'polar';
figure;
pattern(Array, Freq3D , 'PropagationSpeed', PropagationSpeed,...
 'Type','directivity', 'CoordinateSystem', format,'weights', w(:,1),...
  'ShowArray',false,'ShowLocalCoordinates',true,...
  'ShowColorbar',true,'Orientation',[0;0;0]);

% Find the weights
w = zeros(getNumElements(Array), length(Frequency));
for idx = 1:length(Frequency)
    SteerVector = phased.SteeringVector('SensorArray', Array,...
      'PropagationSpeed', PropagationSpeed, ...
      'NumPhaseShifterBits', PhaseShiftBits(idx));
    w(:, idx) = step(SteerVector, Frequency(idx), SteeringAngles(:, idx));
end

% Plot 2d azimuth graph
format = 'polar';
cutAngle = 0;
figure;
pattern(Array, Frequency, -180:180, cutAngle, 'PropagationSpeed', PropagationSpeed,...
 'Type', 'directivity', 'CoordinateSystem', format ,'weights', w);

% Find the weights
w = zeros(getNumElements(Array), length(Frequency));
for idx = 1:length(Frequency)
    SteerVector = phased.SteeringVector('SensorArray', Array,...
       'PropagationSpeed', PropagationSpeed,...
       'NumPhaseShifterBits', PhaseShiftBits(idx));
    w(:, idx) = step(SteerVector, Frequency(idx), SteeringAngles(:, idx));
end


% Plot 2d azimuth graph
format = 'rectangular';
cutAngle = 0;
figure;
pattern(Array, Frequency, -180:180, cutAngle, 'PropagationSpeed', PropagationSpeed,...
 'Type', 'directivity', 'CoordinateSystem', format ,'weights', w);


steeringvector = phased.SteeringVector('SensorArray',Array);
sv = steeringvector(Frequency,SteeringAngles);