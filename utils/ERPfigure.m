function hFigure = ERPfigure(varargin)

% ERPfigure set keypress functions for a figure which has or will have plots in it
% f = ERPfigure opens a new figure with several keypress functionalities and returs a handle to the figure. 
% use the new figure as you would use any figure started with the figure
% command.
% Type h when the figure is in focus to see a list of options to use with
% the figure. 
%
% ERPfigure(handle) adds those functionalities to an existing figure

% V1.0 9 Nov 2012 Leon Deouell HUJI 
% V1.1 14 Jan 2013 Leon Deouell: fixed so that legends are not affected by axes
%          manipulations
% V1.2 30 Nov 2013 Leon Deouell: copy to new window now requires a double
%          click on the axes
% v1.3 28 Mar 2016 Leon Deouell: added option to have a moving vertical
%      line with time label 
%       5 Apr 2016 Leon Deouell: only one time label in the corner of the figure, rather than 
%                  on each axes. The text can be dragged. Multiple time
%                  lines are allowed (by using the + key function again and
%                  again. To delete a line - use the righ key
% v1.4 19 Oc 2017 Edden Gerber: Added two lines of code in the
% LocalCopyToNewWindow function to also copy the original figure's
% colormap. Also added semicolons everywhere to try to suppress the output
% to the console but it didn't work. 
% v1.41 31 Oc 2017 Edden Gerber: allowed "PolarAxes" object in addition to
% "Axes". 
% 
% This function calls for FIGEXTRAS 

if nargin<1 %this is the case when we start a new figure
    hFigure = figure;
    set(hFigure,'keypressFcn','ERPfigure(get(gcf,''CurrentCharacter''));');
    set(hFigure,'windowbuttondownfcn','if strcmp(get(gcf,''selectiontype''),''open''), ERPfigure(''CopyToNewWindow''), end');
    textmessage = uicontrol('style','pushbutton', 'string', 'Press h to get options; click this text to hide',...
         'position',[2 2 300 20],...
         'tag','init',...
         'visible','off',...
         'callback','set(gco,''visible'',''off'')');
    try
        figextras;
    catch
        warning('to get the extra image functions you need the figextras function on your path');
    end
    
    return

elseif isfig(varargin{1}) 
    hFigure = varargin{1};
    set(hFigure,'keypressFcn','ERPfigure(get(gcf,''CurrentCharacter''));');
    set(hFigure,'windowbuttondownfcn','if strcmp(get(gcf,''selectiontype''),''open''), ERPfigure(''CopyToNewWindow''), end');
    try
        figextras;
    end
    return
    
elseif ischar(varargin{1})%this is the case when we pass key presses
    hFigure = gcf;
%     key = get(hFigure,'CurrentCharacter');
    key = varargin{1};
    
    if isempty(findobj(hFigure,'type','axes')) && isempty(findobj(hFigure,'type','polaraxes'))
        return
    end
    helpstring = {'Options:','',...
                 '---------',...
                 'Double-clicking an axis will open it in a new window',...
                 '',...
                  'up\down arrow   zoom in\out y scale (make waveform bigger\smaller)',...
                  'y                         set y axis limits',...
                  'x                         set x axis limits',...
                  'c                         set c (color) limits',...
                  't                         set ticks and tick labels',...
                  '| (vertical)            toggle y grid',...
                  '- (dash)                toggle x grid',...
                  'z                         toggle zero lines',...
                  'h                         start this help dialog',...
                  '+                         add time cursor (click right button on the line to delete)'};
                  

    switch key
        case 30
           uparrow(hFigure);
        case 31
           downarrow(hFigure);
        case 'y'
          setlimY(hFigure);
        case 'x'
          setlimX(hFigure);

        case 'c'
          setlimC(hFigure);
        case 't'
          setticks(hFigure);
        case '|'
           setXgrid(hFigure);
        case '-'
           setYgrid(hFigure);
        case 'z'
           setZeroLines(hFigure);
        case 'CopyToNewWindow'
            LocalCopyToNewWindow(gca);
        case 'h'
            m = msgbox(helpstring, 'ERPfigure help');
            set(findobj(m,'tag','MessageBox'),'fontsize',9, 'fontweight','bold'); % a hack to get the font larger on a messagebox
        case '+'
            setTimeline(hFigure);
    end
    
else
    help mfilename
    error('ERPfigure: unknown input argument');
end
  

function subplots = find_subplots(hFigure)
subplots = findobj(hFigure,'type','axes','-not','tag','legend');

        

function downarrow(hFigure)
subplots = find_subplots(hFigure);
sample = get(subplots(1),'ylim');
newrange = sample * 1.1;
set(subplots,'ylim', newrange);


function uparrow(hFigure)
subplots = find_subplots(hFigure);
sample = get(subplots(1),'ylim');
newrange = sample * 0.9;
set(subplots,'ylim', newrange);

function setlimY(hFigure)
subplots = find_subplots(hFigure);
sample = get(subplots(1),'ylim');
defaultanswer{2} = num2str(sample(1));
defaultanswer{1} = num2str(sample(2));
answer = inputdlg({'High Y', 'Low Y'},'Input new Y limits',1,defaultanswer);
if isempty(answer)
    return;
end
HighY = str2num(answer{1}); %#ok<*ST2NM>
LowY = str2num(answer{2});
if HighY <= LowY
    errordlg('High should be higher then Low...');
    return
end
newrange = [LowY HighY];
set(subplots,'ylim', newrange);

function setlimC(hFigure)
subplots = find_subplots(hFigure);
sample = get(subplots(1),'clim');
defaultanswer{2} = num2str(sample(1));
defaultanswer{1} = num2str(sample(2));
answer = inputdlg({'High C', 'Low C'},'Input new Y limits',1,defaultanswer);
if isempty(answer)
    return;
end
HighC = str2num(answer{1}); %#ok<*ST2NM>
LowC = str2num(answer{2});
if HighC <= LowC
    errordlg('High should be higher then Low...');
    return;
end
newrange = [LowC HighC];
set(subplots,'clim', newrange);

function setlimX(hFigure)
subplots = find_subplots(hFigure);
sample = get(subplots(1),'xlim');
defaultanswer{2} = num2str(sample(1));
defaultanswer{1} = num2str(sample(2));
answer = inputdlg({'High X', 'Low X'},'Input new X limits',1,defaultanswer);
if isempty(answer)
    return;
end
HighX = str2num(answer{1}); %#ok<*ST2NM>
LowX = str2num(answer{2});
if HighX <= LowX
    errordlg('High should be higher then Low...');
    return;
end
newrange = [LowX HighX];
set(subplots,'xlim', newrange)

function setticks(hFigure)
subplots = find_subplots(hFigure);
% sample = get(subplots(1),'xtick');
answer = inputdlg({'Xticks','Xticklabels(separated by | sign)', 'Yticks (row numbers not the actual values)', 'Yticklabels(separated by | sign)'},...
                'Enter positions of new ticks (type none to remove all)',...
                1,{'','auto', '', 'auto'});
if isempty(answer)
    return;
end

Xtick = answer{1};
Xticklabel = answer{2};

if strcmpi(Xtick,'none')
    set(subplots ,'xtick',[], 'xticklabelmode', 'auto');
else
    newticksX = str2num(Xtick);
    switch Xticklabel
        case {'none','None'}
            newXlabels = [];
            xlabelmode = 'manual';
        case {'auto', 'Auto',''}
            newXlabels = [];
            xlabelmode = 'auto';
        otherwise
            newXlabels = Xticklabel;
%             newXlabels = str2num(Xticklabel);

            xlabelmode = 'manual';
    end
    
    if ~isempty(newticksX)
        set(subplots ,'xtick',newticksX, 'xticklabel',newXlabels, 'xticklabelmode', xlabelmode);
    end
end

% Y = answer{2};

Ytick = answer{3};
Yticklabel = answer{4};

if strcmpi(Ytick,'none')
    set(subplots ,'ytick',[], 'yticklabelmode', 'auto');
else
    newticksY = str2num(Ytick);
    switch Yticklabel
        case {'none','None'}
            newYlabels = [];
            ylabelmode = 'manual';
        case {'auto', 'Auto',''}
            newYlabels = [];
            ylabelmode = 'auto';
        otherwise
            newYlabels = Yticklabel;
            ylabelmode = 'manual';
    end
    
    if ~isempty(newticksY)
        set(subplots ,'ytick',newticksY, 'yticklabel',newYlabels, 'yticklabelmode', ylabelmode);
    end
end

function setXgrid(hFigure)
subplots = find_subplots(hFigure);
sample = get(subplots(1),'xgrid');
switch sample
    case 'off'
        set(subplots, 'xgrid','on');
    case 'on'
        set(subplots, 'xgrid','off');
end

function setYgrid(hFigure)
subplots = find_subplots(hFigure);
sample = get(subplots(1),'ygrid');
switch sample
    case 'off'
        set(subplots, 'ygrid','on');
    case 'on'
        set(subplots, 'ygrid','off');
end

function setTimeline(hFigure)
subplots = find_subplots(hFigure);
currentpointer = get(gcf,'pointer');
currentwbuf = get(gcf,'WindowButtonUpFcn');
set(gcf,'pointer','cross');
set(gcf,'WindowButtonUpFcn','set(gcf,''pointer'',''arrow'')');
waitfor(gcf,'pointer');

lines = findobj(subplots,'tag','timeline');
if ~strcmp(get(gco,'type'),'line');
    return
end
ylim = get(gca,'ylim');
ymin = ylim(1);
C = get (gca, 'CurrentPoint');
x = C(1,1);

for i = 1:length(subplots)
    axes(subplots(i));
%     text(C(1,1)+.1, ymin, num2str(round(x,2)),'tag',ah=findall(gcf,'shapetype','textbox');,'backgroundcolor','y');
    hold on;
%     vlines(i) = plot([nearest_time nearest_time], get(gca,'ylim'),'tag','timeline');
    vlines(i) = plot([x x], get(gca,'ylim'),'tag','timeline');
    col = get(vlines(1),'color');
    set(vlines(i),'ButtonDownFcn',@GetTimeLine);
    set(gcf,'WindowButtonUpFcn','set(gcf,''WindowButtonMotionFcn'','''')');
end
tbox = moveableTextBox(num2str(round(x,2)),[.05 .05 .07 0.07]);
set(tbox,'tag','timeline_text','color',col);
set(vlines,'userdata',{vlines, tbox});
set(gcf,'pointer',currentpointer);
% set the context menu to delete the lines
hcmenu = uicontextmenu;
set(vlines, 'uicontextmenu', hcmenu);
uimenu(hcmenu, 'Label','Delete lines','callback','l = get(gco,''userdata'');, delete([l{1},l{2}])');




function GetTimeLine(hObject, callbackdata)
% set(gcf,'WindowButtonMotionFcn',['lines = get(gco,''userdata'');,moveline(lines)'])
set(gcf,'WindowButtonMotionFcn',{@MoveTimeLine,get(gco,'userdata')});

function MoveTimeLine(hObject, callbackdata,handles)
%lines is an array - the first element is the line handles and the second
%is the related textbox 
C = get (gca, 'CurrentPoint');
x = C(1,1);
lines = handles{1};
t = handles{2}
set(lines,'xdata', [x x]);
%t = findall(gcf,'tag','timeline_text');
% pos = get(t(1),'position');
% set(t,'string',num2str(x), 'position',[round(x,2), pos(2), pos(3)])
set(t,'string',num2str(x));


function setZeroLines(hFigure)
subplots = find_subplots(hFigure);
init = findobj(gcf,'tag','init');
if ~isempty(init) && ~isempty(get(init,'userdata')) 
    u = get(init,'userdata');
    if isfield(u,'zlines')
        if strcmp(get(u.zlines.x(1),'visible'),'on')
            set([u.zlines.x u.zlines.y],'visible','off');
        else
            set([u.zlines.x u.zlines.y],'visible','on');
        end
        return;
    end
end

x = zeros(length(subplots),1);
y = zeros(length(subplots),1);
for i = 1:length(subplots)
    axes(subplots(i));
    [x(i), y(i)] = localzerolines(gca,':');
end
u.zlines.x = x;
u.zlines.y = y;
set(init,'userdata',u);

function [xzero, yzero] = localzerolines(ax, style)
%[xzero, yzero]  = zerolines(ax) plots horizontal and vertical zerolines in the axes ax
% if no axes handle is provided, it plots in the current axes object.
%xzero is the handle to the horizontal zeroline
%yzero is the handle to the vertical zeroline
%Use these handles to delete or make the lines invisible

if nargin == 0
    ax = gca;
    style = 'k:';
elseif nargin == 1
    style = 'k:';
elseif nargin == 2
    style = ['k' style];
end

hold on;

for a = 1:length(ax)
    
    curax = ax(a);
    axes(curax);
    hold on;
    xzero = plot(get(curax,'xlim'), [ 0 0 ],style) ;%xline
    yzero = plot([ 0 0], get(curax,'ylim'),style); % yline
    
end

function f=LocalCopyToNewWindow(a)

%Copy current axes to a new figure window with all its children and prperties
%(except position)

cm = colormap(gcf); % get current colormap

if ~exist('a', 'var')
    a = gca;
end
% f = ERPfigure;
f = eval([mfilename ';']); %create a new figure using this function 
set(f,'visible','off');
new = copyobj(a, f);
set(new, 'units','normal','pos',[.1 .1 .8 .8]);
set(f, 'visible','on');

colormap(f,cm); % set same colormap as original figure

function answer=isfig(number)

%ISFIG(NUMBER) returns 1 if number is a handle to a figure, 0 otherwise

l=findobj('type','figure');
try
    if ~isempty(find(l==number)), answer = 1; else, answer = 0; end
catch
    answer = 0;
end

function tbox = moveableTextBox(string, position)

%create text box in the figure that can be moved
if ~exist('position','var')
    position = [.05 .05 .07 .07];
end
% I use annotation below and not just text because text needs an axes and
% annotation does not (actually it does, but the annotation function
% creates one and hides it in the back ground apprently so it saves some
% lines of code here)
tbox = annotation('textbox','string', string, 'LineStyle','none','tag','timeline_text');
set(tbox, 'ButtonDownFcn',@moveObjwithPointer);
set(tbox,'position',position);


function moveObjwithPointer(hObject, callbackdata)
% move the current object to where the pointer is
curWinMotionFcn = get(gcf,'WindowButtonMotionFcn'); %store for later
curWinButtonUpFcn = get(gcf,'WindowButtonUpFcn'); %store for later
curPosition = get(hObject,'position');
set(gcf,'WindowButtonMotionFcn',@setnewpos);
set(gcf,'WindowButtonUpFcn',{@restoreFcns, curWinMotionFcn, curWinButtonUpFcn});

function restoreFcns(hObject, callbackdata, motionfcn,upfcn )
set(gcf,'WindowButtonMotionFcn',motionfcn);
set(gcf,'WindowButtonUpFcn',upfcn);

function newposition = setnewpos(hObject, callbackdata)
set(gco,'units',get(gcf,'units')); %make sure we speak the same language of units
boxposition = get(gco,'position');
newposition = [get(gcf,'CurrentPoint')-boxposition(3:4)/2 boxposition(3:4)]; % keep the size as before. We only want to move
% newposition = get(gcf,'CurrentPoint') ; % keep the size as before. We only want to move
set(gco,'position',newposition);