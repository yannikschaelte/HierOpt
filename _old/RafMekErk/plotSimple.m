function fh = plotSimple(parameters, varargin)
% plotMultiStarts plots the result of the multi-start 
% optimization stored in parameters.
%
% USAGE:
% fh = plotMultiStarts(parameters)
% fh = plotMultiStarts(parameters,fh)
% fh = plotMultiStarts(parameters,fh,options)
%
% plotMultiStarts() uses the following PestoPlottingOptions members:
%  * PestoPlottingOptions::add_points
%  * PestoPlottingOptions::title
%  * PestoPlottingOptions::draw_bounds
%
% Parameters:
% parameters: parameter struct containing information about parameters
%   and log-posterior.
% varargin:
% fh: handle of figure in which profile likelihood is plotted. If no
%   figure handle is provided, a new figure is opened.
% options: options of plotting as instance of PestoPlottingOptions
%
% Return values:
% fh: figure handle
%
% History: 
% * 2012/05/31 Jan Hasenauer
% * 2016/10/07 Daniel Weindl

%% CHECK AND ASSIGN INPUTS

% Check, if parameters has all necessary fieds
% parameters = checkSanityOfStructs(parameters, 'parameters');

% Open figure
if length(varargin) >= 1 && ~isempty(varargin{1}) && isvalid(varargin{1})
    fh = figure(varargin{1});
else
    fh = figure('Name','plotMultiStarts');
end

% Options
defaultOptions = PestoPlottingOptions();
defaultOptions.add_points.par = [];
defaultOptions.add_points.logPost = [];
defaultOptions.add_points.col = [0,0.8,0];
defaultOptions.add_points.ls = '-';
defaultOptions.add_points.lw = 1;
defaultOptions.add_points.m = 'd';
defaultOptions.add_points.ms = 8;
defaultOptions.add_points.name = 'add. point';

if length(varargin) >= 2
    options = varargin{2};
    options = handlePlottingOptionArgument(options);
    options = setdefault(options, defaultOptions);
    options.add_points = setdefault(options.add_points, defaultOptions.add_points);
else
    options = defaultOptions;
end

%% SORT RESULTS
% [parameters] = sortMultiStarts(parameters);

%% CLUSTERING
n_starts = length(parameters.MS.logPost);
if (n_starts > 1)
    clust = cluster(linkage(pdist(parameters.MS.logPost)),'cutoff',0.1,'criterion','distance');
else
    clust = 1;
end
uclust = unique(clust);
    
for iclust = 1:length(uclust)
    sizecluster(iclust) = sum(clust == uclust(iclust));
end

%% ASSIGN COLORS
Col = colormap(gray(n_starts+ceil(n_starts/3)));
Col = Col.^(1/3);
Col(1,:) = [1,0,0];

% sort clusters
for iclust = 1:length(uclust)
    Jclust(iclust) = max(parameters.MS.logPost(find(clust == uclust(iclust))));
end
Jclust(isnan(Jclust)) = -Inf;
[~,idx] = sort(Jclust,'descend');
uclust = uclust(idx);
sizecluster = sizecluster(idx);

if(sizecluster(1)>1)
    ColClust = [1,0,0;flipud(parula(max(sum(sizecluster>1)-1,0)))];
else
    ColClust = flipud(parula(sum(sizecluster>1)));
end

for iclust = 1:length(uclust)
    if(sizecluster(iclust)>1)
    Col(clust == uclust(iclust),:) = repmat(ColClust(sum(sizecluster(1:iclust)>1),:),[sizecluster(iclust),1]);
    end
end

%% PLOT OBJECTIVES
% subplot(2,2,1);
n_finished_starts = 0;
for j = 1 : n_starts
    if ~isnan(parameters.MS.logPost(j))
        n_finished_starts = j;
    else
        break;
    end
end

plot(1:n_finished_starts,parameters.MS.logPost(1:n_finished_starts),'-','color',0.9*[1,1,1],'linewidth',2);
hold on;
for j = n_finished_starts:-1:1
    plot(j,parameters.MS.logPost(j),'o','color',Col(j,:),'linewidth',2);
    hold on;
end
if ~isempty(options.add_points.logPost)
    if length(options.add_points.logPost) == 1
        if (n_starts == 1)
            addPointsX = [0.85 1.15];
            addPointsY = options.add_points.logPost * [1 1];
        else
            addPointsX = 1 : n_starts;
            addPointsY = options.add_points.logPost * ones(size(1:n_starts));
        end
        plot(addPointsX,addPointsY,options.add_points.ls,'color',...
            options.add_points.col(1,:),'linewidth',options.add_points.lw); 
        hold on;
    else
        plot(1:length(options.add_points.logPost),options.add_points.logPost,...
            options.add_points.ls,'color',options.add_points.col(1,:),...
            'linewidth',options.add_points.lw,'marker',options.add_points.m,...
            'markersize',options.add_points.ms); hold on;        
    end
end
hold off;
xlim([1-0.2,n_starts+0.2]);
xlabel('start');
ylabel('log-likelihood');
if options.title
    title('all estimates');
end

ylim([-150,Inf]);