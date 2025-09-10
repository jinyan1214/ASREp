clear 
close all
% Append paths to plot functions
addpath('../../Postprocessing_in_FEM')
% Load mesh
load('./input/burd_mat_20250901T162943.mat');
% Load tensile strain
el_strains = readmatrix('./output/principal_strain_elastic.txt');
ep_strains = readmatrix('./output/principal_strain_elastoplastic.txt');
% Load nodal disp
el_u = readmatrix('./output/nodal_displacement_elastic.txt');
ep_u = readmatrix('./output/nodal_displacement_elastoplastic.txt');
% Plot elastic strain
figure
PlotFieldonDefoMesh(wholeNodesXYZ, wholeElem2n, 100,...
    [el_u(1:3:end),el_u(2:3:end),el_u(3:3:end)],el_strains(:,1)*1000000,...
    'Principal tensile strains($\mu\varepsilon$)')
daspect([1 1 1])
saveas(gcf, './output/strain_plot_elastic.png')
saveas(gcf, './output/strain_plot_elastic.fig')
% Plot elastoplastic strain
figure
PlotFieldonDefoMesh(wholeNodesXYZ, wholeElem2n, 100,...
    [ep_u(1:3:end),ep_u(2:3:end),ep_u(3:3:end)],ep_strains(:,1)*1000000,...
    'Principal tensile strains($\mu\varepsilon$)')
daspect([1 1 1])
saveas(gcf, './output/strain_plot_elastoplastic.png')
saveas(gcf, './output/strain_plot_elastoplastic.fig')


