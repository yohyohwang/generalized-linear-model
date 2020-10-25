% Script hipp_glm.m
% MATLAB shell to fit a GLM model of the relation between
% spiking and the rat's position, history, and network effects,
% and visualize the spatial component of this model.
%
% The code is initialized with an overly simple GLM model construction.
% Please improve it!

% load the rat trajectory and spiking data;
load('hipp_data.mat');

% fit a GLM model to the x and y positions.  Adjust this line to improve
% model fit.
[b,dev,stats] = glmfit([xN yN],spikes,'poisson');

%visualize your model
% construct a grid of positions to plot the model against...
figure;
[x_new,y_new]=meshgrid(-1:.1:1);

% compute lambda for each point on this grid using the GLM model
lambda = exp( b(1) + b(2)*x_new + b(3)*y_new );
lambda(find(x_new.^2+y_new.^2>1))=nan;

%plot lambda as a function of position over this grid
h_mesh = mesh(x_new,y_new,lambda,'AlphaData',0);
hold on;
plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));

