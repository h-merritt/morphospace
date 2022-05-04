%% morphospace analysis

%clear all
%close all
%clc

% load data
structd = load('schaefer400_structfunc_data_fixgs.mat')
nnode = 400;
nsub = 92;
datam = zeros(nnode,nnode,nsub);

addpath("/Users/haily/Documents/projects/social/startup+code")
%subjm is 1x92 double array of subject labels
subjm = cell(nsub,1);
for isub = 1:size(info,1)
    subjm{isub} = structd.sf_hcp_data(isub).sub_name;
end
subjm = str2double(subjm);
subjm = subjm';
% create list of subs in structd
s_nsub=95;
struct_subs = zeros(1,s_nsub);
for i = 1:nsub
    struct_subs(:,i) = str2num(structd.sf_hcp_data(i).sub_name);
end
for is = 1:length(subjm)
    index = find(struct_subs==subjm(is));
    datam(:,:,is) = structd.sf_hcp_data(index).struct.avg;
end

%% compute measures
% initialize matrices to store measures
% average clustering coefficient
cluster_mean = zeros(nsub,1);
cluster_var = zeros(nsub,1);
% average shortest path length
path = zeros(nsub,1);
% density
dens = zeros(nsub,1);
% degree
degree_mean = zeros(nsub,1);
degree_var = zeros(nsub,1);
% treeness
%t = zeros(nsub);
% orderability
%o = zeros(nsub);
% efficiency
%e = zeros(nsub,1);

% loop through subs and compute all measures 
for i = 1:nsub
    cluster_mean(i) = mean(clustering_coef_bu(datam(:,:,i)));
    cluster_var(i) = var(clustering_coef_bu(datam(:,:,i)));
    path(i) = charpath(datam(:,:,i));
    dens(i) = density_und(datam(:,:,i));
    degree_mean(i) = mean(degrees_und(datam(:,:,i)));
    degree_var(i) = var(degrees_und(datam(:,:,i)));
    %e(i) = efficiency_bin(data(:,:,i));
    disp(i)
end

% get degrees for figures
for i = 1:nsub
    degree = degrees_und(data(:,:,i));
    density = density_und(data(:,:,i));
    data_nrm = weight_conversion(data, 'normalize');
    cluster_mean = mean(clustering_coef_bu(data_nrm(:,:,i)));
    cluster_var = var(clustering_coef_bu(data_nrm(:,:,i)));
end

degrees = zeros(nnode,nsub);
cluster = zeros(nnode,nsub);
for i = 1:nsub
    degrees(:,i) = degrees_und(data(:,:,i));
    data_nrm = weight_conversion(data, 'normalize');
    cluster(:,i) = clustering_coef_bu(data_nrm(:,:,i));
end
imagesc(degrees)
imagesc(max(degrees))
degree_var = var(degrees);
plot(degree_var)
imagesc(cluster)
imagesc(max(cluster))
cluster_v = var(cluster);
plot(cluster_v)
scatter(degree_var,cluster_v)
[r,p] = corrcoef(degree_var,cluster_v);
% corr coef = -0.1149, p = 0.2675

%% calculate variability of each
%vars = [var(acc_mean),var(acc_var),var(aspl),var(d),var(degree_mean),var(degree_var),var(e)];
vars = [var(cluster_mean),var(cluster_var),var(path),var(dens),var(degree_mean),var(e)];
var_acc_mean = var(acc_mean);
var_acc_var = var(acc_var);
var_aspl = var(aspl);
var_d = var(d);
var_degree_mean = var(degree_mean);
var_degree_var = var(degree_var);
var_e = var(e);
% plot variability of measures
plot(vars)

% all measures in one matrix
morpho = zeros(92,6);
morpho(:,1) = degree_mean;
morpho(:,2) = degree_var;
morpho(:,3) = path;
morpho(:,4) = cluster_mean;
morpho(:,5) = cluster_var;
morpho(:,6) = dens;
save('morpho.mat','morpho');

%% compare correlations
measures = [cluster_mean,cluster_var,aspl,d,degree_mean,degree_var,e];
corrs = zeros(size(measures,2));
for m = 1:size(measures,2)
    for n = 1:size(measures,2) 
        c = corrcoef(measures(:,m),measures(:,n));
        corrs(m,n) = c(2);
    end
end

imagesc(corrs)

% get group-averaged data
group = zeros(ncortex,ncortex);
for i = 1:ncortex
    for j = 1:ncortex
        group(i,j) = mean(data_un(i,j,:));
    end
end

%% plot network in anatomical space

% density
dens = 0.01;

% threshold network
SCthr = threshold_proportional(data_un(:,:,1),dens);

% get coor
coor = struct.sf_hcp_data(1).parc.coords;

% get edge list
[ex,ey,ez] = adjacency_plot_und(SCthr,coor);

% make plot
f = figure(...
    'units','inches',...
    'position',[5,5,4,4]);
plot3(ex,ey,ez);
hold on;
scatter3(coor(:,1),coor(:,2),coor(:,3),'filled');
axis image;


%% behavioral
% extract behavioral measures for the HCP100UR subjects
pheno=readtable('restricted_HCP_betzel.txt');
[~,idx,~] = intersect(pheno.Subject,subjm);
tblp = pheno(idx,:);
tblb = behav(idx,:);
vars_morpho = [tblb.PercStress_Unadj,tblb.ASR_Anxd_Raw,tblb.DSM_Depr_Raw,tblb.CogFluidComp_AgeAdj,tblp.Age_in_Yrs,tblp.SSAGA_Educ];
save('phenoVars.mat','vars_morpho');