function commission_ROBOTIC_VUC

global expeTabs ...
    expePanelTabs ...
    expeNum...
    data...
    aniScreen ...
    histScreen ...
    parHistory...
    env...
    ctrl...
    kpiHistory...
    kpiHistoryData...
    K1_min K1_max K2_min K2_max K3_min K3_max K4_min K4_max bo_results ...
    
kpiHistoryData = [];

expeNum = 0;

data = struct;

data.parameters = [];
data.Xt = [];
data.Ut = [];
data.Xv = [];
data.Qv = [];
data.KPI = [];
data.noise = [];
data.success = [];

addpath('backend');

screen = get(0,'screensize');

W = screen(3);
H = screen(4);

% set main figure
main = figure(...
    'Name',         'SBO REXPEK 2022-2026 : Robotic Virtual Use Case -- Version 3. UGent-D2LAB/FM-MIRO ©.',...
    'Position',     [W/2-2*W/6 H/2-2*H/6 2*W/3 2*H/3],...
    'NumberTitle',  'off');

W = 2*W/3;
H = 2*H/3;

% set control panel
dashBoard = uipanel(...
    'Parent',       main,...
    'Position',     [0.02 0.02 .376 0.96]);

W1 = .367*W;
H1 = 0.96*H;

% set history tables
ctrlPanel = uipanel(...
    'Parent',       dashBoard,...
    'Position',     [.01 .61 .98 .39],...
    'Title',        'Settings');

hstryPanel1 = uipanel(...
    'Parent',       dashBoard,...
    'Position',     [.01 .31 .98 .3],...
    'Title',        'Parameter history');

hstryPanel2 = uipanel(...
    'Parent',       dashBoard,...
    'Position',     [.01 .01 .98 .3],...
    'Title',        'KPI history');

columnNames = cell(1,4); for i = 1:4; columnNames{i} = ['K',num2str(i)]; end
parHistory = uitable(...
    'Parent',       hstryPanel1,...
    'ColumnName',   columnNames, ...
    'Units',        'normalized', ...
    'Position',     [.01 .02 .98 .96]);

columnNames = {'Cost','Success'};
tooltip = ['Cost: A unitless value that is calculated from the trajectory points. You are optimizing this one.' newline ...
    'Success: Did we manage to push the slider into the vicinity of the red dot or not?'];
kpiHistory = uitable(...
    'Parent',       hstryPanel2,...
    'ColumnName',   columnNames, ...
    'Tooltip',      tooltip, ...
    'Units',        'normalized',...
    'Position',     [.01 .02 .98 .96]);

% set experiment panel
expePanel = uipanel(...
    'Parent',       main,...
    'Position',     [0.04+0.376 0.02 0.564 0.96],...
    'Title',        'Experiment dashboard');
expePanelTabs = uitabgroup(...
    'Parent',       expePanel);
expeTabs = struct;

W2 = 0.564*W;
H2 = 0.96*H;

K1_min = .1;
K1_max = 3;
K2_min = .1;
K2_max = 3;
K3_min = .1;
K3_max = 3;
K4_min = .1;
K4_max = 3;

% set experiment

[env,ctrl] = env1;

% set sliders and text boxes

uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'text', ...
    'Units',        'normalized',...
    'Position',     [.05 .825 .1 .1] , ...
    'String',       'K1');
slider1 = uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'slider', ...
    'Min',          K1_min, ...
    'Max',          K1_max, ...
    'Value',        K1_max, ...
    'Units',        'normalized',...
    'Position',     [.2 .825 .6 .1], ...
    'Callback',     @sliderCb1);
text1 = uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'edit', ...
    'Units',        'normalized',...
    'Position',     [.85 .825 .1 .1], ...
    'String',       num2str(get(slider1, 'Value')), ...
    'Callback',     @textCb1);

uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'text', ...
    'Units',        'normalized',...
    'Position',     [.05 .675 .1 .1] , ...
    'String',       'K2');
slider2 = uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'slider', ...
    'Min',          K2_min, ...
    'Max',          K2_max, ...
    'Value',        K2_max, ...        
    'Units',        'normalized',...
    'Position',     [.2 .675 .6 .1], ...
    'Callback',     @sliderCb2);
text2 = uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'edit', ...
    'Units',        'normalized',...
    'Position',     [.85 .675 .1 .1], ...
    'String',       num2str(get(slider2, 'Value')), ...
    'Callback',     @textCb2);

uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'text', ...
    'Units',        'normalized',...
    'Position',     [.05 .525 .1 .1] , ...
    'String',       'K3');
slider3 = uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'slider', ...
    'Min',          K3_min, ...
    'Max',          K3_max, ...
    'Value',        K3_max, ...        
    'Units',        'normalized',...
    'Position',     [.2 .525 .6 .1], ...
    'Callback',     @sliderCb3);
text3 = uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'edit', ...
    'Units',        'normalized',...
    'Position',     [.85 .525 .1 .1], ...
    'String',       num2str(get(slider3, 'Value')), ...
    'Callback',     @textCb3);

uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'text', ...
    'Units',        'normalized',...
    'Position',     [.05 .375 .1 .1] , ...
    'String',       'K4');
slider4 = uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'slider', ...
    'Min',          K4_min, ...
    'Max',          K4_max, ...
    'Value',        0.5, ...        
    'Units',        'normalized',...
    'Position',     [.2 .375 .6 .1], ...
    'Enable',       'off', ...
    'Callback',     @sliderCb4);
text4 = uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'edit', ...
    'Units',        'normalized',...
    'Position',     [.85 .375 .1 .1], ...
    'String',       num2str(get(slider4, 'Value')), ...
    'Enable',       'off', ...
    'Callback',     @textCb4);

handles.texts(1) = text1;
handles.texts(2) = text2;
handles.texts(3) = text3;
handles.texts(4) = text4;

handles.sliders(1) = slider1;
handles.sliders(2) = slider2;
handles.sliders(3) = slider3;
handles.sliders(4) = slider4;

handles.main = main;

guidata(main,handles);

% set push buttons

uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'pushbutton',...
    'Units',        'normalized',...
    'Position',     [.1 .1 .2 .2], ...
    'String',       'evaluate',...
    'Tooltip',      ['Run the experiment once with the given parameters, in order to evaluate the KPIs.' newline 'You practically press this button, then change the parameters, then press this button again, and so on.'],...
    'Callback',     @runExperiment);

if ~isempty(env.robot)
uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'pushbutton',...
    'Units',        'normalized',...
    'Position',     [.3 .1 .2 .2], ...    
    'String',       'animate',...
    'Tooltip',      'This will open up a 3D animation of the pusher-slider robot in another view. Just a fancy animation to see.',...
    'Callback',     @aniExperiment);
end

uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'pushbutton',...
    'Units',        'normalized',...
    'Position',     [.5 .1 .2 .2], ...
    'String',       'history',...
    'Tooltip',      'This will open up a plot of the history of the parameter values and KPIs in another window.', ...
    'Callback',     @visHist);

uicontrol(...
    'Parent',       ctrlPanel,...
    'Style',        'pushbutton',...
    'Units',        'normalized',...
    'Position',     [.7 .1 .2 .2], ...
    'Tooltip',      ['When you are done with the tuning, please press this.' newline 'This will save all the results to a .mat file that you can give us.'], ...
    'String',       'save&quit',...
    'Callback',     @stop);

%comment out this button here in the final version to not have BO on the UI:

uicontrol(...
    'Parent', ctrlPanel,...
    'Style', 'pushbutton',...
    'Units', 'normalized',...
    'Position', [.9 .1 0.09 .2],...
    'String', 'BO',...
    'Visible', 'off',...
    'Callback', @runBO);

end

% auxiliary function

function stop(hObject,~)
global...
    data...
    aniScreen...
    histScreen

answer = questdlg(...
        'Would you like to save the results and quit?', ...
        'SBO REXPEK 2022-2026 : Robotic Virtual Use Case -- Version 3. UGent-D2LAB/FM-MIRO ©.', ...
        'Yes.','No.','Yes.');
    switch answer
        case 'Yes.'
            name = inputdlg('Please enter your name (to be used in the filename):','SBO REXPEK 2022-2026 : Robotic Virtual Use Case -- Version 3. UGent-D2LAB/FM-MIRO ©.');
            name = name{:};
            name = regexprep(name, '[^\w\s]', '_');
            name = regexprep(name, '\s', '_');
            filename = ['Results_' name '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.mat'];
            disp(['Saving to:' filename])
            save(filename,'data');
            msgbox(['Thanks. Please send us the file: ' filename],'SBO REXPEK 2022-2026 : Robotic Virtual Use Case -- Version 3. UGent-D2LAB/FM-MIRO ©.','none');

            if exist('aniScreen') && ~isempty(aniScreen) && isvalid(aniScreen)
                close(aniScreen);
            end
            if exist('histScreen') && ~isempty(histScreen) && isvalid(histScreen)
                close(histScreen);
            end
            close(guidata(hObject).main);
        case 'No.'
    end
end

function [cost] = BO_controller_cost(params)
    global iteration ctrl env bo_iterations_data expeTabs expePanelTabs
    if isempty(iteration), iteration = 0; end
    iteration = iteration + 1;

    ctrl.par(1) = params.K1;
    ctrl.par(2) = params.K2;
    ctrl.par(3) = params.K3;
    ctrl.par(4) = 0.5;

    % Run the simulation
    [bo_iterations_data(iteration).Xt,bo_iterations_data(iteration).Ut,~,bo_iterations_data(iteration).noise,bo_iterations_data(iteration).success] = simulateEnvironment(env,ctrl);
    bo_iterations_data(iteration).KPI = postProcess(bo_iterations_data(iteration).Xt,bo_iterations_data(iteration).Ut,env,bo_iterations_data(iteration).success);
    
    % Calculate my KPIs
    N_steps = size(bo_iterations_data(iteration).Xt,2);
    differences_along_dim_2 = diff(bo_iterations_data(iteration).Xt(1:2,:),1,2);
    distances = sqrt(sum(differences_along_dim_2.^2,1));
    total_distance = sum(distances);
    final_distance = sqrt((bo_iterations_data(iteration).Xt(1,end)-env.goal.x).^2+(bo_iterations_data(iteration).Xt(2,end)-env.goal.y).^2);

    %cost1 = min(100,total_distance);
    %cost2 = min(100,(1/13.)*N_steps*env.time.dt);
    cost = bo_iterations_data(iteration).KPI(1);

    % Save results to CSV
    new_data_line = [cost; 0; ...
        0; 0; 0; 0;  ...
        total_distance; final_distance; N_steps; bo_iterations_data(iteration).success ; params.K1; params.K2; params.K3; 0.5];
    writematrix(new_data_line.', 'current_run/BO_current_run.csv', 'WriteMode', 'append');
    
    if ~isempty(fields(expeTabs))
        expeTabs.tab(end+1) = uitab(expePanelTabs,'Title',['BO' num2str(iteration)]);
        visLastExp(expeTabs.tab(end),bo_iterations_data(iteration).Xt,bo_iterations_data(iteration).Ut,env,bo_iterations_data(iteration).KPI,[params.K1,params.K2,params.K3,params.K4],bo_iterations_data(iteration).success);
        expePanelTabs.SelectedTab = expeTabs.tab(end);
    end
end


function runBO(hObject,~)
    global K1_min K1_max K2_min K2_max K3_min K3_max K4_min K4_max bo_results
    % Execute Bayesian Optimization
    diary current_run/matlab_log.txt
    % Define the parameter bounds
    params = optimizableVariable('K1',[K1_min,K1_max],'Transform','none');
    params(2) = optimizableVariable('K2',[K2_min,K2_max],'Transform','none');
    params(3) = optimizableVariable('K3',[K3_min,K3_max],'Transform','none');
    %params(4) = optimizableVariable('K4',[K4_min,K4_max],'Transform','none');
    !rm current_run/BO_current_run.csv
    !echo "cost, constraint_violation, a*KPI1, b*KPI3, KPI1, KPI3, total_distance, final_distance, N_steps, bo_iterations_data(iteration).success , params.K1, params.K2, params.K3, params.K4">current_run/BO_current_run.csv 
    bo_results = bayesopt(@(params) BO_controller_cost(params), params, 'MaxObjectiveEvaluations', 1000, 'AcquisitionFunctionName','expected-improvement-plus','IsObjectiveDeterministic',true);
    save('current_run/bo_results.mat','bo_results');
    diary off
end

% experiment & postprocessing

function runExperiment(hObject,~)

global...
    expeTabs...
    expePanelTabs...
    expeNum...
    data...
    env...
    ctrl...
    parHistory...
    kpiHistory...
    histScreen...
    kpiHistoryData

expeNum = expeNum+1;

handles = guidata(hObject);

data(expeNum).parameters = zeros(1,4);
for i = 1:4
ctrl.par(i) = get(handles.sliders(i),'value');
data(expeNum).parameters(i) = ctrl.par(i);
end

[data(expeNum).Xt,data(expeNum).Ut,~,data(expeNum).noise,data(expeNum).success] = simulateEnvironment(env,ctrl);
data(expeNum).KPI = postProcess(data(expeNum).Xt,data(expeNum).Ut,env,data(expeNum).success);

expeTabs.tab(expeNum) = uitab(expePanelTabs,'Title',num2str(expeNum));
 
parameters = zeros(expeNum,4);
colors = zeros(expeNum,3);
kpis = zeros(expeNum,4);
rowNames = cell(1,expeNum);
for n = 1:expeNum
    idx = expeNum-n+1;
    parameters(n,:) = data(idx).parameters;
    cost = data(idx).KPI(1);
    if data(idx).success
        colors(n,:) = [1 1 1];
    else
        colors(n,:) = [0.8588 0.5725 0.5490];
    end
    kpis(n,:) = [ cost, 0, 0, data(idx).success ];
    rowNames{n} = num2str(idx);
end
data(expeNum).CostTimeDistanceSuccess = kpis(1,:); %add the current experiment to it (it's the first item in the list)

% Filter for rows with success = 1
successRows = kpis(:,4) == 1;
successCosts = kpis(successRows, 1);

% Check if there are any successful rows
if any(successRows)
    % Find the index of the minimum cost among successful rows
    [~, minSuccessIndex] = min(successCosts);
    
    % Find the corresponding original index
    successIndices = find(successRows);
    minIndex = successIndices(minSuccessIndex);
    
    % Define light green color
    lightGreen = [0.7 1 0.7];
    
    % Update colors to highlight the minimum cost row with success = 1
    colors(minIndex, :) = lightGreen;
end


set(parHistory,'Data',parameters);
set(parHistory,'Rowname',rowNames);
set(parHistory,'BackgroundColor',colors);

set(kpiHistory,'Data',[kpis(:,1) kpis(:,4)]);
set(kpiHistory,'Rowname',rowNames);
set(kpiHistory,'BackgroundColor',colors);

if expeNum > 5
    kpis5 = kpis(1:5,:);
else
    kpis5 = kpis;
end

visLastExp(expeTabs.tab(expeNum),data(expeNum).Xt,data(expeNum).Ut,env,kpis5',data(expeNum).parameters,data(expeNum).success)
expePanelTabs.SelectedTab = expeTabs.tab(expeNum);

kpiHistoryData = kpis;
if exist('histScreen') && ~isempty(histScreen) && isvalid(histScreen)
    visHist
end

end

function values = postProcess(Xt,Ut,env,success)

T = size(Xt,2);
dd = zeros(2,T);
Vt = zeros(2,T-1);
for t = 1:T
    R = rot(Xt(3,t));
    dd(:,t) = R*[Xt(4,t);-env.geom.b/2-env.geom.r] + Xt(1:2,t);
    if t < T
        Vt(:,t) = R*Ut(:,t);
    end
end
distance = sqrt((env.initial.x - env.goal.x)^2+(env.initial.y - env.goal.y)^2);
energy = 1e2;

% duration
value1 = T*env.time.dt;
% slider distance
value2 = sum(sqrt(sum(diff(Xt(1:2,:),1,2).^2,1)));
% pusher distance
value3 = sum(sqrt(sum(diff(dd,1,2).^2,1)));
% estimated energy 1
value4 = sum(sum(diag([1 1 1/env.model.beta^2])*diff(Xt(1:3,:),1,2).^2))/energy;


v1 = T*env.time.dt;  
%v2 = sum(sqrt(sum(diff(Xt(1:2,:),1,2).^2,1)));  
v3 = sum(sqrt(sum(diff(dd,1,2).^2,1)));  
v4 = sum(sum(Ut(1,:).^2))*env.time.dt;  
v5 = sum(sum(Ut(2,:).^2))*env.time.dt;

obj = log(min(20,v1+2*v3+1*v5+2e1*v4+1e3*(1-success)));

values = [obj;0;0;0];

end

% visualization tools

function visLastExp(tab,Xt,Ut,env,value,parameters,success)

sgtitle(['Experiment with parameters (K1,K2,K3,K4) = (',...
    num2str(round(parameters(1),1)),',',...
    num2str(round(parameters(2),1)),',',...
    num2str(round(parameters(3),1)),',',...
    num2str(round(parameters(4),1)),')'],'parent',tab,'fontsize',12)

subplot(2,2,1,'Parent',tab);

hold on;
Xt_wrappedToPi = Xt;
Xt_wrappedToPi(3,:) = wrapToPi(Xt_wrappedToPi(3,:));
plot((0:size(Xt,2)-1)*env.time.dt,Xt_wrappedToPi');
plot([1 1]*(size(Xt,2)-1)*env.time.dt,[-1 1]*pi,'k--');
plot([0 env.time.T],[env.goal.x env.goal.x],'--','color',[0 0.4470 0.7410]);
plot([0 env.time.T],[env.goal.y env.goal.y],'--','color',[0.8500 0.3250 0.0980]);
legend('x','y','\theta','d'); xlabel('time [s]'); ylabel(['states [m],[m],[rad],[m]'])
grid on; box on; axis([0 env.time.T -pi pi]);
title('States')

subplot(2,2,2,'Parent',tab);

hold on;
plot((0:size(Ut,2)-1)*env.time.dt,Ut');
plot([1 1]*(size(Xt,2)-1)*env.time.dt,[-.5 .5],'k--');
legend('v','w'); xlabel('time [s]'); ylabel('input velocities [m/s]');
grid on; box on; axis([0 env.time.T -.5 .5]);
title('Inputs')

subplot(2,2,3,'Parent',tab);
if success == false
    set(gca, 'Color', [0.8588 0.5725 0.5490]); % RGB equivalent of #db928c
end
visTraj2D(Xt,env,[0 1 0],[-.6 .6 -.2 1]);
xlabel('x [m]'); ylabel('y [m]');
title('Trajectory')

subplot(2,2,4,'Parent',tab);
bar(value(1,:));
xticklabels({'last','-1','-2','-3','-4'});
title(['Cost' newline '(last 5)'])

end

function visTraj2D(Xt,env,color,window)
a = env.geom.a;
b = env.geom.b;
r = env.geom.r;
if nargin < 3
    color = 'g';
end
circ.x = cos(linspace(0,2*pi,20));
circ.y = sin(linspace(0,2*pi,20));
[~,T] = size(Xt);
hold on; daspect([1 1 1]);
% vis2D(Xt(1,1),Xt(2,1),Xt(3,1),a,b,color);
fill(env.goal.x+circ.x*0.05,env.goal.y+circ.y*0.05,'r','edgecolor','none');
dd = zeros(2,T);
ls = zeros(2,T);
for t = union(1:3:T,T)
    vis2D(Xt(1,t),Xt(2,t),Xt(3,t),a,b,color);
    dd(:,t) = rot(Xt(3,t))*[Xt(4,t);-b/2-r] + Xt(1:2,t);
    h = fill(circ.x*r+dd(1,t),circ.y*r+dd(2,t),'y'); set(h,'facealpha',.5);
    plot([0 ls(1,t)]+dd(1,t),[0 ls(2,t)]+dd(2,t),'k','linewidth',2);
    axis(window);
end
box on;
grid on;
end

function visHist(~,~)

global...
    histScreen...
    data...
    expeNum...
    kpiHistoryData

parameters = zeros(expeNum,4);
kpis = zeros(expeNum,4);
ticks = cell(expeNum,1);
for n = 1:expeNum
    parameters(n,:) = data(n).parameters;
    kpis(n,:) = data(n).KPI;
    idx = expeNum-n+1;
    ticks{n} = num2str(idx);
end

if isempty(histScreen)
    histScreen = figure(...
        'Name',         'SBO REXPEK 2022-2026 : Robotic Virtual Use Case -- Version 3. UGent-D2LAB/FM-MIRO ©.',...
        'NumberTitle',  'off');
elseif ~isvalid(histScreen)
    histScreen = figure(...
        'Name',         'SBO REXPEK 2022-2026 : Robotic Virtual Use Case -- Version 3. UGent-D2LAB/FM-MIRO ©.',...
        'NumberTitle',  'off');
else
    figure(histScreen);
end

clf
sgtitle('history','fontsize',12);

subplot(2,1,1,'Parent',histScreen);

plot(1:expeNum,fliplr(parameters'),'o--','linewidth',2,'Markersize',5);
legend('K1','K2','K3','K4'); xticks(1:expeNum); xticklabels(ticks);
grid on; box on;

subplot(2,1,2,'Parent',histScreen);

if ~isempty(kpiHistoryData)
    plot(1:expeNum,kpiHistoryData(:,1)','*--','linewidth',2,'Markersize',5);
    legend('cost','\Delta_s','\Delta_p','E'); xticks(1:expeNum); xticklabels(ticks);
    grid on; box on;
end
end

% animation tools

function aniExperiment(~,~)

global ...
    data... 
    expePanelTabs...
    env...
    aniScreen...
    expeNum

if expeNum > 0

Num = str2double(expePanelTabs.SelectedTab.Title);

fps = 4;

if isempty(data(Num).Xv)
    data(Num).Xv = interp1(linspace(0,1,length(data(Num).Xt)),data(Num).Xt',linspace(0,1,round(length(data(Num).Xt)*env.time.dt*fps)))';
    data(Num).Qv = taskTraj2jointTraj(data(Num).Xv,env);
end

if isempty(aniScreen)
    aniScreen = figure(...
        'Name',         'SBO REXPEK 2022-2026 : Robotic Virtual Use Case -- Version 3. UGent-D2LAB/FM-MIRO ©.',...
        'NumberTitle',  'off');
elseif ~isvalid(aniScreen)
    aniScreen = figure(...
        'Name',         'SBO REXPEK 2022-2026 : Robotic Virtual Use Case -- Version 3. UGent-D2LAB/FM-MIRO ©.',...
        'NumberTitle',  'off');
else
    figure(aniScreen);
end

aniTraj3D(data(Num).Xv,data(Num).Qv,env,[0 1 0],[-.6 .6 -.2 1 0 .7],[150 30],aniScreen); %,'bad.gif');

end

end

function aniTraj3D(Xv,Qv,env,color,window,azel,aniScreen)

r = env.geom.r;

xg = .05*cos(linspace(0,2*pi));
yg = .05*sin(linspace(0,2*pi));

[xc,yc,zc] = cylinder; xc = r*xc; yc = r*yc; zc = .1*zc;

T = size(Xv,2);
dd = zeros(2,T);
for t = 1:T
    dd(:,t) = rot(Xv(3,t))*[Xv(4,t);-env.geom.b/2-env.geom.r] + Xv(1:2,t);
end

for t = 1:size(Xv,2)
    if t == 1
        camproj('perspective');
    end
    if isvalid(aniScreen)
        figure(aniScreen);
        clf
        sgtitle('animation','fontsize',12);
        hold on; daspect([1 1 1]); box on; grid on;
        show(env.robot,Qv(:,t));
        fill(env.goal.x+xg,env.goal.y+yg,'r','edgecolor','none');
        plot3(dd(1,1:t),dd(2,1:t),zeros(1,t),'k','linewidth',1);
        plot3(Xv(1,1:t),Xv(2,1:t),zeros(1,t),'color',color,'linewidth',1);
        vis3D(Xv(1,t),Xv(2,t),Xv(3,t),env.geom.a,env.geom.b,color);
        surf(xc+dd(1,t),yc+dd(2,t),zc,'edgecolor','none');
        camproj('perspective'); view(azel);
        axis(window);
        light('style','local','position',[0 .5 1]);
        drawnow
    else
        break
    end
end
end

function Qt = taskTraj2jointTraj(Xt,env)
ik = inverseKinematics('RigidBodyTree',env.robot);
weights = [1 1 1 1 1 1];
Q0 = [pi/2 -pi/4 pi/2 -3*pi/4 -pi/2 0]';
T0 = getTransform(env.robot,Q0,'tool0');
Qt = zeros(6,size(Xt,2));
for t = 1:size(Xt,2)
    pos = [rot(Xt(3,t))*[Xt(4,t);-env.geom.b/2-env.geom.r] + Xt(1:2,t);0.1];
    T = T0; T(1:3,4) = pos;
    [Qt(:,t),~] = ik('tool0',T,weights,Q0);
    Q0 = Qt(:,t);
end
end

% slider & text callbacks

function sliderCb1(hObject,~)
handles = guidata(hObject);
newVal = get(hObject,'Value');
set(handles.texts(1),'String',num2str(newVal));
end

function textCb1(hObject,~)
global K1_min K1_max
handles = guidata(hObject);
newVal = str2double(get(hObject, 'String'));
if isnan(newVal) || newVal < K1_min || newVal > K1_max
    errordlg('Please enter a valid number within the slider range.', 'Invalid Input');
    set(hObject, 'String', num2str(get(handles.sliders(1), 'Value')));
else
    set(handles.sliders(1),'Value', newVal);
end
end

function sliderCb2(hObject,~)
    handles = guidata(hObject);
    newVal = get(hObject, 'Value');
    set(handles.texts(2),'String',num2str(newVal));
end

function textCb2(hObject,~)
global K2_min K2_max
handles = guidata(hObject);
newVal = str2double(get(hObject, 'String'));
if isnan(newVal) || newVal < K2_min || newVal > K2_max
    errordlg('Please enter a valid number within the slider range.', 'Invalid Input');
    set(hObject, 'String', num2str(get(handles.sliders(2), 'Value')));
else
    set(handles.sliders(2),'Value', newVal);
end
end

function sliderCb3(hObject,~)
handles = guidata(hObject);
newVal = get(hObject, 'Value');
set(handles.texts(3),'String',num2str(newVal));
end

function textCb3(hObject,~)
global K3_min K3_max
handles = guidata(hObject);
newVal = str2double(get(hObject, 'String'));
if isnan(newVal) || newVal < K3_min || newVal > K3_max
    errordlg('Please enter a valid number within the slider range.', 'Invalid Input');
    set(hObject, 'String', num2str(get(handles.sliders(3), 'Value')));
else
    set(handles.sliders(3),'Value', newVal);
end
end

function sliderCb4(hObject,~)
handles = guidata(hObject);
newVal = get(hObject, 'Value');
set(handles.texts(4),'String',num2str(newVal));
end

function textCb4(hObject,~)
global K4_min K4_max
handles = guidata(hObject);
newVal = str2double(get(hObject, 'String'));
if isnan(newVal) || newVal < K4_min || newVal > K4_max
    errordlg('Please enter a valid number within the slider range.', 'Invalid Input');
    set(hObject, 'String', num2str(get(handles.sliders(4), 'Value')));
else
    set(handles.sliders(4),'Value', newVal);
end
end

