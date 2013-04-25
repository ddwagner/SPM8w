function varargout = art(varargin)
% see help2/art_readme.txt for most recent tutorial.
%
%art - module for automatic and manual detection and removal of outliers.
%This utility asks for SPM.mat file and a text motion parameters file.
%It then displays four graphs. The top graph is the global brain
%activation mean as a function of time.
%The second is a z-normalized (stdv away from mean) global brain activation
%as a function of time.
%The Third shows the linear motion parameters (X,Y,Z) in mm as a function of time
%and the fourth shows the rotational motion parameters (roll,pitch,yaw) in radians as
%a function of time.
%using default threshold values for each of the bottom three graphs we
%define outliers as points that exceed the threshold in at least one of
%these graphs. The thresholds are shown as horizontal black lines in each of the graphs.
%Points which are identified as outliers, are indicated by a vertical red
%line in the graph that corresponds to the outlying parameter(s). For
%example, the if the absolute value of the Y motion parameter for time t=17 is above
%the motion threshold, it is identified as an outlier and indicated by a
%red vertical line at t=17 in the third graph. The union of all
%outliers is indicated by vertical lines on the top graph. The list of
%outliers is also displayed in the editable text box below the graphs.
%The current values of the thresholds are displayed by the side of the
%corresponding graphs. These values may can be changed by the user either
%by pressing the up/down buttons, which increment/decrement the current
%value by 10%, or by specifying a new value in the text box.
%In Addition, the user can manually add or remove points from the list of
%outliers by editting the list. Note that the list is only updated once the
%curser points outside the text box (i.e. click the mouse somewhere outside
%the text box). Since any changes made by the user are overridden once the
%thresholds are updated, it is recommended to do any manual changes as the
%last step before saving.
%Pressing the save button lets the user choose wheter to save the
%motion statistics (.mat or .txt) the list of outliers (.mat or .txt),
%or save the graphs (.jpg, .eps or matlab .fig).
%
% ----------------------------------------------------------------------
% - Switched out range function to remove stat toolbox dependency - DDW March/13
% - First attempt to spm8'ify the code. It's a little hacky hacky but 
%   so far it works DDW.March 2012
% - got joe moran's version and made some small changes to make it
%   matlab6.5 compatbile. Also added a subject selection bit and automated
%   some stuff. DDW.March 2008
% - added "show correlations" and "show spectrum"
% - minor GUI changes to support Windows and open a large graph in a separate
%   window. Also fixed starnge motion params filename bug. Shay Mozes, 5/2/2007
% - fixed bug in display of motion outlier on the graph, Shay Mozes 4/30/2007
% - superimpose task conditions on the z-graph, Shay Mozes, 4/24/2007
% - added support for SPM5, Shay Mozes, 4/2007
% - now supporting FSL .par format, 4/9/2007
% - new GUI and features Shay Mozes, 2006 + Mar. 2007
% from art_global.m, by Paul Mazaika, April 2004.
% from artdetect4.m, by Jeff Cooper, Nov. 2002
% from artdetect3.m, by Sue Whitfield
% artdetect.m  Sue Whitfield 2000

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @art_OpeningFcn, ...
    'gui_OutputFcn',  @art_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before art is made visible.
% This is where the logic begins
function art_OpeningFcn(hObject, eventdata, handles, varargin)

% -----------------------
% Initialize, begin loop
% -----------------------
warning('OFF','all')
pfig    = [];
spm_ver = 8;
global num_sess
p                = varargin{1};
num_sess         = varargin{2}; %get numofruns from batcher - DDW

%Set figure title
thefig = gcf;
set(thefig,'Name',sprintf('Art Outlier Detection | Subject: %s Run %d',p.subj,num_sess));

% ------------------------
% Default values for outliers
% ------------------------
%  Deviations over 2*std are outliers.
z_thresh          = 2; 

%movement thresholds - shay
mvmt_thresh       = 1;    %outliers have linear momevent norm > mm_thresh
rotat_thresh      = 0.05; %outliers have angular movement > rad_thresh in at least one orientation - change this???
mvmt_diff_thresh  = 0.5;
rotat_diff_thresh = 0.05;

% ------------------------
% Collect files
% ------------------------
global_type_flag = 1;

%global_type_flag = spm_input('Which global mean to use?', 1, 'm', ...
%    'Regular | Every Voxel | User Mask | Auto ( Generates ArtifactMask )',...
%    [1 2 3 4], 1);

P              = [];
M              = [];
mv_data        = [];
output_path    = p.func;
drop_flag = 0;
P         = spm_select('FPList','FUNCTIONAL',['^',p.boldtok,num2str(num_sess),'.*\.nii']);
M         = load(fullfile(p.func,['rp_',p.rptok,num2str(num_sess),'.txt']));
mv_data   = vertcat(mv_data,M); %Redundant with M - DDW

% -------------------------
% get file identifiers and Global values
% -------------------------
fprintf('%-4s: ','Mapping files...')
VY        = spm_vol(P);
fprintf('%3s\n','...done')
nscans = size(VY,1);

% ------------------------
% Compute Global variate
%--------------------------
%GM     = 100;
g      = zeros(nscans,1);
fprintf('%-4s: %3s','Calculating globals...',' ')
% Regular global mean
for i  = 1:nscans
    g(i) = spm_global(VY(i));
end
fprintf('...done\n')


% ------------------------
% Compute default out indices by z-score
% ------------------------
gsigma = std(g);
gmean = mean(g);
pctmap = 100*gsigma/gmean;
z_thresh = 0.1*round(z_thresh*10); % Round to nearest 0.1 Z-score value

%update text fields
set(handles.data_stdv,'String',num2str(gsigma));
set(handles.zthresh,'String',num2str(z_thresh));

%HACK TO AUTOSET MVMT THRESHOLDS TO DIFFERENCE -DDW
set(handles.mvthresh,'String',num2str(mvmt_thresh));
set(handles.rtthresh,'String',num2str(rotat_thresh));
%set(handles.mvthresh,'String',num2str(mvmt_diff_thresh));
%set(handles.rtthresh,'String',num2str(rotat_diff_thresh));

%columns 8:13 store the difference series
mv_data(2:end,8:13)  = abs(mv_data(2:end,1:6) - mv_data(1:end-1,1:6));

%HACK TO AUTOSET MOVEMENT DATA TO DIFFERENCE - DDW
%mv_data(:,1:6) = mv_data(:,8:13);

for i=1:size(mv_data,1)
    %7th column holds euclidean norms of movement
    mv_data(i,7) = norm(mv_data(i,1:3));
    %14-15th column stores the sums/norms of linear and angular motion
    %    mv_data(i,14) = norm(mv_data(i,8:10));
    %    mv_data(i,15) = norm(mv_data(i,11:13));
end

%save application data for use in callbacks
setappdata(handles.zthresh,'g',g);
setappdata(handles.mvthresh,'mv_data',mv_data);
setappdata(handles.mvthresh,'altval',num2str(mvmt_diff_thresh));
setappdata(handles.rtthresh,'altval',num2str(rotat_diff_thresh));
%setappdata(handles.savefile,'data',P);
setappdata(handles.savefile,'path',output_path);

%plot global mean
axes(handles.globalMean);
%%%===DDW EDIT=== (removed dependency on stat toolbox)
%rng = range(g);
rng = max(g)-min(g);
%%%===DDW EDIT===
plot(g);

ylabel(['Range = ' num2str(rng)], 'FontSize', 8);
if ( global_type_flag == 1 ) title('Global Mean - Regular SPM'); end
if ( global_type_flag == 2 ) title('Global Mean - Every Voxel'); end
if ( global_type_flag == 3 ) title('Global Mean - User Defined Mask'); end
if ( global_type_flag == 4 ) title('Global Mean - Generated ArtifactMask'); end

%plot in stddev
z_Callback(hObject, eventdata, handles,1.0);

%plot movement
mv_Callback(hObject, eventdata, handles,1.0);

%plot rotation
rt_Callback(hObject, eventdata, handles,1.0);

%calculate all outliers and plot
calc_all(hObject, eventdata, handles)

%calculate and print statistics of movement

mv_stats = [mean(abs(mv_data)); std(abs(mv_data)); max(abs(mv_data)) ];
%global_stats = [gmean, gsigma];
setappdata(handles.mvthresh,'mv_stats',mv_stats);

fprintf('\n\nGlobal statistics -  mean: %7.4f stdv: %7.4f',gmean,gsigma);
fprintf('\n\nStatistics of movement data:\n\n');
fprintf('%5s%10s%10s%10s%11s%10s%9s%10s\n',' ','x','y','z','norm',' pitch','roll','yaw');
fprintf('%7s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n','mean ',mv_stats(1,1:3),mv_stats(1,7),mv_stats(1,4:6));
fprintf('%7s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n','stdv ',mv_stats(2,1:3),mv_stats(2,7),mv_stats(2,4:6));
fprintf('%7s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n\n\n','max ',mv_stats(3,1:3),mv_stats(3,7),mv_stats(3,4:6));

% Choose default command line output for art
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


%function that retreives all different outliers and outputs
%a single list of outliers
function calc_all(hObject, eventdata, handles)

%get data
idx = getappdata(handles.zthresh,'zoutliers');
idx = [idx , getappdata(handles.mvthresh,'mv_x_outliers')];
idx = [idx , getappdata(handles.mvthresh,'mv_y_outliers')];
idx = [idx , getappdata(handles.mvthresh,'mv_z_outliers')];
idx = [idx , getappdata(handles.rtthresh,'rt_p_outliers')];
idx = [idx , getappdata(handles.rtthresh,'rt_r_outliers')];
idx = [idx , getappdata(handles.rtthresh,'rt_y_outliers')];

if (get(handles.norms,'Value') == get(handles.norms,'Max'))
    idx = [idx , getappdata(handles.rtthresh,'rt_norm_outliers')];
    idx = [idx , getappdata(handles.mvthresh,'mv_norm_outliers')];
end

idx = unique(idx);

%update data
set(handles.all_outliers, 'String', int2str(idx));

%plot
all_outliers_Callback(hObject, eventdata, handles);

function all_outliers_Callback(hObject, eventdata, handles)
% Plots the outliers in the global graph based on the current
% outliers in the all_outliers edit field
%get data
g = getappdata(handles.zthresh,'g');
%    idx = round(str2num(get(handles.all_outliers,'String')));
tmps = get(handles.all_outliers,'String');
if length(tmps) > 0
    idx = round(str2num(tmps(1,:)));
    [nstrings,stam] = size(tmps);
    for i=2:nstrings
        idx = [idx , round(str2num(tmps(i,:)))];
    end
    set(handles.all_outliers, 'String', int2str(idx));
else
    idx = [];
end

%plot global mean
axes(handles.globalMean);
cla;
%%%===DDW EDIT=== (removed dependency on stat toolbox)
%rng = range(g);
rng = max(g)-min(g);
%%%===DDW EDIT===
plot(g);
l=ylabel(['Range = ' num2str(rng)], 'FontSize', 8);

y_lim = get(gca, 'YLim');
for i = 1:length(idx)
    line((idx(i)*ones(1, 2)), y_lim, 'Color', 'red');
end

% --- Executes on button press in z_up.
function z_up_Callback(hObject, eventdata, handles)
z_Callback(hObject, eventdata, handles, 1.05);
calc_all(hObject, eventdata, handles)

function z_figure(hObject, eventdata, handles)
incr=1.0;
h = figure;
%get data
z_thresh = str2num(get(handles.zthresh,'String'));
g = getappdata(handles.zthresh,'g');

%calc new outliers
z_thresh = z_thresh*incr;
%%%===DDW EDIT===
%out_idx = (find(abs(zscore(g)) > z_thresh))';
out_idx = (find(abs(statsw_zscore(g)) > z_thresh))';
%%%===DDW EDIT===

%update text
set(handles.zthresh,'String',num2str(z_thresh));
setappdata(handles.zthresh,'zoutliers',out_idx);

%update plot
%%%===DDW EDIT===
%plot((zscore(g)));
plot((statsw_zscore(g)));
%%%===DDW EDIT===
l=ylabel('stdv away \newlinefrom mean');
thresh_x = 1:length(g);
thresh_y = z_thresh*ones(1,length(g));
line(thresh_x, thresh_y, 'Color', 'black');
line(thresh_x, -1*thresh_y, 'Color', 'black');

axes_lim = get(gca, 'YLim');
axes_height = axes_lim;
for i = 1:length(out_idx)
    line((out_idx(i)*ones(1, length(axes_height))), axes_height, 'Color', 'black');
end

if (get(handles.showDesign,'Value') == get(handles.showDesign,'Max'))
    [SPM,design] = get_design(handles);
    hold on
    colors = {'go','ro','co','mo','yo'};
    for i=1:size(design,2)
        plot(1:size(design,1) , design(:,i),colors{mod(i,5)+1},'MarkerSize',4);        
    end
    hold off
end

%function for getting SPM design matrix information
function [SPM,design] = get_design(handles)
    SPM = getappdata(handles.showDesign,'SPM');
    tmpfile = spm_select(1,'^.*\.mat$',['Select design matrix:']);
    load(tmpfile);
    setappdata(handles.showDesign,'SPM',SPM);
    sessions = getappdata(handles.showDesign,'sessions');
    if (isempty(sessions))
        tmpsess = inputdlg('What session(s) to use? (e.g. 1 or [1,2])');
        sessions = eval(char(tmpsess));
        setappdata(handles.showDesign,'sessions',sessions);
    end
    rows = [];
    cols = [];
    for s = sessions
        rows = [rows SPM.Sess(s).row];
        cols = [cols SPM.Sess(s).col];
    end
  
    design = SPM.xX.X(rows,cols);
    
% base callback for (re)plotting outliers graph
function z_Callback(hObject, eventdata, handles, incr)

%get data
z_thresh = str2num(get(handles.zthresh,'String'));
g = getappdata(handles.zthresh,'g');

%calc new outliers
z_thresh = z_thresh*incr;
%%%===DDW EDIT===
%out_idx = (find(abs(zscore(g)) > z_thresh))';
out_idx = (find(abs(statsw_zscore(g)) > z_thresh))';
%%%===DDW EDIT===
%update text
set(handles.zthresh,'String',num2str(z_thresh));
setappdata(handles.zthresh,'zoutliers',out_idx);

%update plot
axes(handles.zvalue);
cla;
%%%===DDW EDIT===
plot((zscore(g)));
%plot((statsw_zscore(g)));
%%%===DDW EDIT===
l=ylabel('stdv away \newlinefrom mean');
set(l,'VerticalAlignment','Bottom');
set(gca,'XTickLabel',[]);

%plot threshold lines
thresh_x = 1:length(g);
thresh_y = z_thresh*ones(1,length(g));
line(thresh_x, thresh_y, 'Color', 'black');
line(thresh_x, -1*thresh_y, 'Color', 'black');

%plot outliers
axes_lim = get(gca, 'YLim');
axes_height = axes_lim;
for i = 1:length(out_idx)
    line((out_idx(i)*ones(1, length(axes_height))), axes_height, 'Color', 'black');
end

%show design
if (get(handles.showDesign,'Value') == get(handles.showDesign,'Max'))
    [SPM,design] = get_design(handles);

    hold on
    colors = {'go','ro','co','mo','yo'};
    for i=1:size(design,2)
        plot(1:size(design,1) , design(:,i),colors{mod(i,5)+1},'MarkerSize',4);        
    end
    hold off

end

% --- Executes on button press in z_down.
function z_down_Callback(hObject, eventdata, handles)

z_Callback(hObject, eventdata, handles,0.95);
calc_all(hObject, eventdata, handles)

function zthresh_Callback(hObject, eventdata, handles)
% hObject    handle to zthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
z_Callback(hObject, eventdata, handles,1.0);
calc_all(hObject, eventdata, handles)


% base callback for (re)plotting movement graphs
function mv_Callback(hObject, eventdata, handles, incr)

%get data
mvmt_thresh = str2num(get(handles.mvthresh,'String'));
mv_data = getappdata(handles.mvthresh,'mv_data');

%calc new outliers
mvmt_thresh = mvmt_thresh*incr;

out_mvmt_idx = (find(abs(mv_data(:,1:3)) > mvmt_thresh))';
out_mvmt_idx_X=[];
out_mvmt_idx_Y=[];
out_mvmt_idx_Z=[];
for i =1:length(out_mvmt_idx)
    if out_mvmt_idx(i) <= length(mv_data)
        out_mvmt_idx_X=[out_mvmt_idx_X; out_mvmt_idx(i)];
    elseif out_mvmt_idx(i) > length(mv_data) & out_mvmt_idx(i) < 2*length(mv_data)
        out_mvmt_idx_Y=[out_mvmt_idx_Y; out_mvmt_idx(i)-length(mv_data)];
    else
        out_mvmt_idx_Z = [out_mvmt_idx_Z; out_mvmt_idx(i)-2*length(mv_data)];
    end
end
out_mvmt_idx_X = out_mvmt_idx_X';
out_mvmt_idx_Y = out_mvmt_idx_Y';
out_mvmt_idx_Z = out_mvmt_idx_Z';

%find norm outliers
if (get(handles.norms,'Value') == get(handles.norms,'Max'))
    normv = zeros(length(mv_data),1);
    for i=1:length(mv_data)
        normv(i) = norm(mv_data(i,1:3));
    end
    out_mvmt_idx_norm = find(normv>mvmt_thresh);
    out_mvmt_idx_norm = out_mvmt_idx_norm';
    setappdata(handles.mvthresh,'mv_norm_outliers',out_mvmt_idx_norm);
end


%update text
set(handles.mvthresh,'String',num2str(mvmt_thresh));
setappdata(handles.mvthresh,'mv_x_outliers',out_mvmt_idx_X);
setappdata(handles.mvthresh,'mv_y_outliers',out_mvmt_idx_Y);
setappdata(handles.mvthresh,'mv_z_outliers',out_mvmt_idx_Z);

axes(handles.mvmtGraph);
cla;
if (get(handles.norms,'Value') == get(handles.norms,'Max'))
    plot(normv);
else
    plot(mv_data(:,1:3));
end
l=ylabel('movement \newline   [mm]');
set(l,'VerticalAlignment','Bottom');
set(gca,'XTickLabel',[]);
if (get(handles.norms,'Value') == get(handles.norms,'Max'))
    l=legend('norm','Location','East');
else
%    l=legend('x mvmt', 'y mvmt', 'z mvmt','Location','East');
    l=legend('x', 'y', 'z','Location','East');
end
%set(l,'Position',[116.4 15.4 12 5.1]);
h = gca;
set(h,'Ygrid','on');

thresh_mv_x = 1:length(mv_data);
thresh_mv_y = mvmt_thresh*ones(1,length(mv_data));
line(thresh_mv_x, thresh_mv_y, 'Color', 'black');
line(thresh_mv_x, -1*thresh_mv_y, 'Color', 'black');

axes_lim = get(gca, 'YLim');
axes_height = axes_lim;
for i = 1:length(out_mvmt_idx_X)
    line((out_mvmt_idx_X(i)*ones(1, length(axes_height))), axes_height, 'Color', 'black');
end
for i = 1:length(out_mvmt_idx_Y)
    line((out_mvmt_idx_Y(i)*ones(1, length(axes_height))), axes_height, 'Color', 'black');
end
for i = 1:length(out_mvmt_idx_Z)
    line((out_mvmt_idx_Z(i)*ones(1, length(axes_height))), axes_height, 'Color', 'black');
end

if (get(handles.norms,'Value') == get(handles.norms,'Max'))
    for i = 1:length(out_mvmt_idx_norm)
        line((out_mvmt_idx_norm(i)*ones(1, length(axes_height))), axes_height, 'Color', 'black');
    end
end

% --- Executes on button press in mv_up.
function mv_up_Callback(hObject, eventdata, handles)
mv_Callback(hObject, eventdata, handles,1.05);
calc_all(hObject, eventdata, handles)

% --- Executes on button press in mv_down.
function mv_down_Callback(hObject, eventdata, handles)
mv_Callback(hObject, eventdata, handles,0.95);
calc_all(hObject, eventdata, handles)

% --- Executes on update of mvthresh.
function mvthresh_Callback(hObject, eventdata, handles)
mv_Callback(hObject, eventdata, handles,1.0);
calc_all(hObject, eventdata, handles)

% base callback for (re)plotting rotation graphs
function rt_Callback(hObject, eventdata, handles, incr)

%get data
rotat_thresh = str2num(get(handles.rtthresh,'String'));
mv_data = getappdata(handles.mvthresh,'mv_data');

%calc new outliers
rotat_thresh = rotat_thresh*incr;
out_rotat_idx = (find(abs(mv_data(:,4:6)) > rotat_thresh))';
out_rotat_idx_p=[];
out_rotat_idx_r=[];
out_rotat_idx_y=[];
for i =1:length(out_rotat_idx)
    if out_rotat_idx(i) <= length(mv_data)
        out_rotat_idx_p=[out_rotat_idx_p; out_rotat_idx(i)];
    elseif out_rotat_idx(i) > length(mv_data) & out_rotat_idx(i) < 2*length(mv_data)
        out_rotat_idx_r=[out_rotat_idx_r; out_rotat_idx(i)-length(mv_data)];
    else
        out_rotat_idx_y = [out_rotat_idx_y; out_rotat_idx(i)-2*length(mv_data)];
    end
end
out_rotat_idx_p = out_rotat_idx_p';
out_rotat_idx_r = out_rotat_idx_r';
out_rotat_idx_y = out_rotat_idx_y';

%find norm/sum outliers
if (get(handles.norms,'Value') == get(handles.norms,'Max'))
    normv = zeros(length(mv_data),1);
    for i=1:length(mv_data)
        normv(i) = norm(mv_data(i,4:6));
    end
    out_rt_idx_norm = find(normv>rotat_thresh);
    out_rt_idx_norm = out_rt_idx_norm';
    setappdata(handles.rtthresh,'rt_norm_outliers',out_rt_idx_norm);
end


%update text
set(handles.rtthresh,'String',num2str(rotat_thresh));
setappdata(handles.rtthresh,'rt_p_outliers',out_rotat_idx_p);
setappdata(handles.rtthresh,'rt_r_outliers',out_rotat_idx_r);
setappdata(handles.rtthresh,'rt_y_outliers',out_rotat_idx_y);

axes(handles.rotatGraph);
cla;
if (get(handles.norms,'Value') == get(handles.norms,'Max'))
    plot(normv);
else
    plot(mv_data(:,4:6));
end
l=ylabel('rotation \newline  [rad]');
set(l,'VerticalAlignment','Bottom');
xlabel('scans');
y_lim = get(gca, 'YLim');
if (get(handles.norms,'Value') == get(handles.norms,'Max'))
    l=legend('norm','Location','East');
else
    l=legend('pitch', 'roll', 'yaw', 'Location', 'East');
end
%set(l,'Position',[116.4 9.4 12 5.1]);
h = gca;
set(h,'Ygrid','on');

thresh_rt_x = 1:length(mv_data);
thresh_rt_y = rotat_thresh*ones(1,length(mv_data));
line(thresh_rt_x, thresh_rt_y, 'Color', 'black');
line(thresh_rt_x, -1*thresh_rt_y, 'Color', 'black');

y_lim = get(gca, 'YLim');
for i = 1:length(out_rotat_idx_p)
    line((out_rotat_idx_p(i)*ones(1, 2)), y_lim, 'Color', 'black');
end
for i = 1:length(out_rotat_idx_r)
    line((out_rotat_idx_r(i)*ones(1, 2)), y_lim, 'Color', 'black');
end
for i = 1:length(out_rotat_idx_y)
    line((out_rotat_idx_y(i)*ones(1, 2)), y_lim, 'Color', 'black');
end

if (get(handles.norms,'Value') == get(handles.norms,'Max'))
    for i = 1:length(out_rt_idx_norm)
        line((out_rt_idx_norm(i)*ones(1, 2)), y_lim, 'Color', 'black');
    end
end



% --- Executes on button press in rt_up.
function rt_up_Callback(hObject, eventdata, handles)
rt_Callback(hObject, eventdata, handles,1.05)
calc_all(hObject, eventdata, handles)

% --- Executes on button press in rt_down.
function rt_down_Callback(hObject, eventdata, handles)
rt_Callback(hObject, eventdata, handles,0.95)
calc_all(hObject, eventdata, handles)

% --- Executes on update of mvthresh.
function rtthresh_Callback(hObject, eventdata, handles)
rt_Callback(hObject, eventdata, handles,1.0)
calc_all(hObject, eventdata, handles)


% --- Executes when difference button uis pressed
function diffs_Callback(hObject, eventdata, handles)

%switch thresholds
tmp = get(handles.mvthresh,'String');
set(handles.mvthresh,'String',getappdata(handles.mvthresh,'altval'));
setappdata(handles.mvthresh,'altval',tmp);
tmp = get(handles.rtthresh,'String');
set(handles.rtthresh,'String',getappdata(handles.rtthresh,'altval'));
setappdata(handles.rtthresh,'altval',tmp);

%switch data used
mv_data = getappdata(handles.mvthresh,'mv_data');
tmp = mv_data(:,1:6);
mv_data(:,1:6) = mv_data(:,8:13);
mv_data(:,8:13) = tmp;

setappdata(handles.mvthresh,'mv_data',mv_data);

%plot movement
mv_Callback(hObject, eventdata, handles,1.0);

%plot rotation
rt_Callback(hObject, eventdata, handles,1.0);

%calculate all outliers and plot
calc_all(hObject, eventdata, handles);


% --- Executes when norms checkbox is changed
function norms_Callback(hObject, eventdata, handles)

%plot movement
mv_Callback(hObject, eventdata, handles,1.0);

%plot rotation
rt_Callback(hObject, eventdata, handles,1.0);

%calculate all outliers and plot
calc_all(hObject, eventdata, handles);

% --- Executes when show design checkbox is changed
function showDesign_Callback(hObject, eventdata, handles)
    %plot movement
    z_Callback(hObject, eventdata, handles,1.0);

    
    
function showCorr_Callback(hObject, eventdata, handles)
if (get(handles.showCorr,'Value') == get(handles.showCorr,'Max'))
%display correlations
    [SPM,design] = get_design(handles);
    mv_data = getappdata(handles.mvthresh,'mv_data');
    sessions = getappdata(handles.showDesign,'sessions');
    f=figure; %changed by ddw for matlab65
    %f = figure();
    setappdata(handles.showCorr,'figure',f);
    for s = sessions
        rows = SPM.Sess(s).row;
        cols = SPM.Sess(s).col;

        %create partial matrix to correlate
        %(we only want to correlate with the motion parameters within each
        %session). NOTE: This may cause weird behaviour in weird designs...
        part = [SPM.xX.X(rows,cols) mv_data(rows,1:6)];
        cm{s} = corrcoef(part);
        a = subplot(length(sessions),1,s);
        
        imagesc(cm{s}(1:end-6,end-5:end));
        colorbar;
        set(a,'XTickLabel',{'x','y','z','pitch','roll','yaw'});
        set(a,'YTick',[1:length(cols)]);
        names = {};
        for i=1:length(cols)
            names(i) = SPM.Sess(s).U(i).name;
        end
        set(a,'YTickLabel',names);
        title(sprintf('Session %d',s));
        
    end
else
    f = getappdata(handles.showCorr,'figure');
    close(f);
end



function showSpec_Callback(hObject, eventdata, handles)
if (get(handles.showSpec,'Value') == get(handles.showSpec,'Max'))
%get data and compute power spectrum
    [SPM,design] = get_design(handles);
    mv_data = getappdata(handles.mvthresh,'mv_data');
    sessions = getappdata(handles.showDesign,'sessions');
    f=figure;    %changed by ddw for matlab65
    %f = figure();
    setappdata(handles.showSpec,'figure',f);
    
    %sampling freq.
    sf = 1/SPM.xY.RT;
    
    for s = sessions
        rows = SPM.Sess(s).row;
        cols = SPM.Sess(s).col;

        %create partial design matrix which only contains relevant data
        %for curent session
        data = [SPM.xX.X(rows,cols) mv_data(rows,1:6)];
                
        cf = sf/2; %Nyquist freq.
        n = size(data,1);
        freqs = (0:cf/n:cf-cf/n)';%these are the descrete freqs used by dct

        %this is done in a loop (and not in matrix ops)
        %since dct encounters memory problems for large matrices.
        hold on
        f = zeros(size(data));
        for i=1:size(data,2)
            %calculate dct
            f(:,i) = dct(data(:,i));
            %normalize
            F(:,i) = f(:,i)/sum(abs(f(:,i)));
        end

        a = subplot(length(sessions), 1, s);
        hs = plot(freqs,f);
        names = {};
        for i=1:length(cols)
            names(i) = SPM.Sess(s).U(i).name;
        end
        names(end+1:end+6) = {'x','y','z','pitch','roll','yaw'};
        l = legend(names,'Location','EastOutside');
        pos = get(l,'Position');
        set(l,'Visible','off');
        for i = 1:length(names)
            color = get(hs(i),'Color');
            box = uicontrol('Style','checkbox','String',names(i),'ForegroundColor', color,'Callback',{@setvisibility},'Value',1,'UserData',hs(i),'Units','normalized');
            tmppos = get(box,'Position');
            tmppos(1:2) = pos(1:2);
            set(box,'Position',tmppos);
            pos(2) = pos(2)+ 0.035;
        end
            
        title(sprintf('Session %d',s));
        xlabel('Frequency [Hz]');
        ylabel('Power(normalized)');
        
        %try finding highpass freq. in SPM.xX.K(s).Hparam
        try
            cutoff = 1/SPM.xX.K(s).HParam;
        catch
            fprintf('no highpass cutoff frequency found in SPM.mat, using default (128).\n');
            cutoff = 1/128;
        end        
        %draw cutoff freq.
        ylim = get(a,'YLim');
        line(cutoff*ones(1,2),ylim,'Color','black');
    end
else
    f = getappdata(handles.showSpec,'figure');
    close(f);
end

%toggle visibility for spectrum graph
function setvisibility(handle,tmp)
    h = get(handle,'UserData');
    if (get(handle,'Value') == get(handle,'Max'))
        set(h,'Visible','on');
    else
        set(h,'Visible','off');
    end


% --- Executes on button press in savefile.
function savefile_Callback(hObject, eventdata, handles)

global num_sess
%save outliers
        out_idx = round(str2num(get(handles.all_outliers, 'String')));
        str=['outliers',num2str(num_sess),'=out_idx;'];eval(str);
        str=['save -ascii -tabs outliers_run',num2str(num_sess),'.txt outliers',num2str(num_sess),';'];       
        eval(str);
        disp(['Outliers for run ',num2str(num_sess),' saved to file outliers_run',num2str(num_sess),'.txt'])
        disp('Press any key to continue with next subject/run')
      
% --- Executes during object creation, after setting all properties.
function all_outliers_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Outputs from this function are returned to the command line.
function varargout = art_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function zthresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function mvthresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function rtthresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


