function extras = figextras(action)
%FIGEXTRAS Add extra functionality to a  figure
%figextras with no inputs adds the following menu items:
%Plot line from image: click on an image and get a 2D plot of the row you
%                                                   clicked
%Plot column from image: click on an image and get a 2D plot of the column you
%                                                   clicked
%figextras(figure_handle) adds the same for the given figure
%Other input arguments are used by the submenus and usually are not to be
%used from the command prompt

%Leon Deouell, HUJI, 19/5/2006
% 4/6/2006 add the situation where the image doesn't have explicit x or y
%                   values 
% 16/4/2017 add the option to scroll through the rows of an image in the
%            plotrow option

if nargin == 0  | ~isstr(action)
    if nargin == 0
        fig = gcf;
    else
        fig = action;
    end
    extras = uimenu('Tag','Extras','Label','Extras');
    Plotline = uimenu('parent',extras,'Label','Plot row from image','callback','figextras(''rowplot'')');
    Plotcolumn = uimenu('parent',extras,'Label','Plot column from image','callback','figextras(''colplot'')');
    Smooth = uimenu('parent',extras,'Label','Smooth','callback','figextras(''smooth'')');
%     Makenice = uimenu('parent',extras,'Label','Make nice plots','callback','figextras(''makenice'')');
    Makenice = uimenu('parent',extras,'Label','Make nice plots');
    Makenice_normal = uimenu('parent',Makenice,'Label','normal (just smooth lines)','callback','figextras(''makenice_normal'')');
    Makenice_big = uimenu('parent',Makenice,'Label','big (smooth lines and enlarge fonts)','callback','figextras(''makenice_big'')');
    Makenice_custom = uimenu('parent',Makenice,'Label','custom','callback','figextras(''makenice_custom'')');
    LineToggle = uimenu('parent',extras,'Label','Line toggle','callback','linetoggle');
    return
end

switch action
    case 'makenice_normal'
       makenice('normal');
    case 'makenice_big'
       makenice('big');
    case 'makenice_custom'
        %get curent properties
        axhandles = findobj('parent',gcf,'type','axes');
        axesfont = num2str(get(axhandles(1),'fontsize'));
        labelsfont = num2str(get( get(axhandles(1),'xlabel'),'fontsize'));
        titlefont = num2str(get( get(axhandles(1),'title'),'fontsize'));
        lines = findobj('parent',axhandles(1),'type','line');
        linewidth = num2str(get(lines(1),'linewidth'));
        %get new properties
        answer = inputdlg({'title','labels','axes','width'},'Enter parameter of gaussian filter',1,...
                            {titlefont, labelsfont, axesfont,linewidth});
        
        if isempty (answer)
            return
        end
        titlefont = answer{1}
        labelsfont = answer{2}
        axesfont = answer{3}
        linewidth = answer{4}
            
        eval(['makenice(''title'', ' titlefont ', ',...
            '''labels'', ' labelsfont ', ',...
            '''axes'', ' axesfont '' ', ',...
            '''width'', ' linewidth ')'])
         
        
                 

    case 'smooth'
        answer = inputdlg({'Kernel size (square if only a scalar is given)','Sigma'},'Enter parameter of gaussian filter');
        
        if isempty (answer)
            return
        end
        kernel = str2num(answer{1});
        sigma = str2num(answer{2});

        if isempty(kernel) || isempty(sigma)
            return
        end
        myfilter = fspecial('gaussian',kernel,sigma);
        % layout = get(handles.Layout,'userdata');
        % h = layout.subplots;
        h = findobj(gcf,'type','image');
        for i = 1: length(h)
            parent = get(h(i),'parent');
            data = get(parent,'userdata');
            if isempty(data) %the data is untouched
                data = get(h(i),'cdata');
                set(parent,'userdata',data)
            end
            fdata = imfilter(data, myfilter);
            set(h(i),'cdata',fdata)
        end        

    case 'rowplot'        

        p = ginput(1);
        x = round(p(1))
        y = round(p(2))
        h = gco;
        if strcmp(get(h,'type'),'image') | strcmp(get(h,'type'),'surface')
            data = get(h,'cdata');
            size_data = size(data);
            xdata = get(h, 'xdata');
            ydata = get(h, 'ydata');
            store.data = data;
            
            if length(xdata) ~= size_data(2)
                xdata = 1:size_data(2); % this happens if x data was not specified when plotting the image
            end
            if length(ydata) ~= size_data(1)
                ydata = 1:size_data(1); % this happens if y data was not specified when plotting the image
            end
            store.ydata = ydata;
            %----find the closest row
            d = ydata - y;
            [m, mInd] = min(abs(d));
            y = mInd;
            %--------------------------
            name = get(h,'tag');
            ERPfigure
            plot(xdata, data(y,:),'tag','waveform','userdata',mInd)
            
            if ~isempty(name) 
                t = [' from ' name];
            else
                t = '';
            end
            hTitle = title([  num2str(ydata(mInd)), ' (Row  ', num2str(y), ')',  t ]);
        else
            error('The point clicked was not an image object')
        end
        store.title = hTitle;
        set(gca,'userdata',store)
        set(gcf,'KeyPressFcn', @stepthrough)
        helpdlg('Press right\left arrows to scroll through the rows of the image')
        
    case 'colplot'
        p = ginput(1);
        x = round(p(1));
        y = round(p(2));
        h = gco;
        if strcmp(get(h,'type'),'image') | strcmp(get(h,'type'),'surface')
            data = get(h,'cdata');
            size_data = size(data);
            xdata = get(h, 'xdata');
            ydata = get(h, 'ydata');
            if length(xdata) ~= size_data(2)
                xdata = 1:size_data(2); % this happens if x data was not specified when plotting the image
            end
            if length(ydata) ~= size_data(1)
                ydata = 1:size_data(1); % this happens if y data was not specified when plotting the image
            end
            %----find the closest column
            d = xdata - x;
            [m, mInd] = min(abs(d));
            x = mInd;
            %--------------------------
            name = get(h,'tag');
            ERPfigure
            plot(ydata, data(:,x))
            if ~isempty(name) 
                t = [' from ' name];
            else
                t = '';
            end
            title([  num2str(xdata(mInd)), ' (Column  ', num2str(x), ')',  t ])
        else
            error('The point clicked was not an image object')
        end
    case linetoggle
     

        %linetoggle   interface for toggling lines in a plot on and off.
        %linetoggle(sourcefig, lines, names)
        %   Opens a figure with checkboxes and a line legend. Currently it assumes
        %   that all the lines are initially visible. Checking off the removes the
        %   line. Checking back restores the line. The names are editable. 
        % Input:
        %  sourcefig - the figure with the lines to toggle (default, if no input is
        %                given, is current figure)
        %   lines - lines handles. If lines is empty, the lineoggle finds the lines
        %                     in the figure using findobj.
        %   names - a cell array of names for the lines. They are presented in
        %   editable fields that can be changed. If names are not given, and the
        %   line objects have string tags the names are dervied from the tags
        %
        %linetoggle with no input arguments uses the lines in the current figure
        %and leaves the names empty - the user can fill in the names. 
        %
        % Note: if the source figure has a name, the linetoggle figure shows the
        % source figure name at the top, so the user can easily link the toggle
        % figure to the source figure. 
        %


        if nargin<1
            sourcefig = gcf;
        end
        if ~exist('lines','var') || isempty(lines)
            lines = findobj(sourcefig,'type','line');
        end
        numlines = numel(lines);
        linetable = struct('handles', [],'color',[],'LineStyle',[], 'LineWidth',[]);
        f= figure('unit','normal', 'position',[.1 .001 .2 .9],'resize','off','name',get(sourcefig,'name'));
        set(f,'units','pixels');
        fsize = get(f,'position');
        fsize = fsize(3:4);
        checkheight = 20 ; %pixels
        if numlines * checkheight > fsize(2)-10
            errordlg('too many lines to work with this function')
        end
        for i = 1:numlines
            h = lines(i);
             linetable(i).handles =h;
             linetable(i).color = get(h,'color');
             linetable(i).LineStyle = get(h,'LineStyle');
             linetable(i).LineWidth = get(h,'LineWidth');
             %create the graphics
             linetable(i).checkhandle = uicontrol(f,'style','checkbox','unit','pixels',...
                 'position', [ checkheight , (fsize(2) - 20 - checkheight*(i-1)), 20 20],...
                 'backgroundcolor',get(gcf,'color'),...
                 'userdata',linetable(i).handles,...
                 'value',1,...
                 'callback',[ 'if get(gco,''value''),',...
                                                'set(get(gco,''userdata''),''visible'',''on''),',...
                                            'else,',...
                                                ' set(get(gco,''userdata''),''visible'',''off''),',...
                                             'end']);


             check_pos = get(linetable(i).checkhandle,'position');
             x(1:2) =[  checkheight *3 , checkheight * 8]/fsize(1); 
             y(1:2) = (check_pos(2)+checkheight/2)/fsize(2);

             %draw the lines with the appropriate proerties 
            linetable(i).linehandle = annotation('line', x, y,...
                 'color',linetable(i).color,...
                 'LineStyle', linetable(i).LineStyle,...
                 'LineWidth',linetable(i).LineWidth);

             %add names if provided
             checkpos = get(linetable(i).checkhandle,'position');
             textpos = checkpos + [checkheight*8, 0, 80, 0];
        %      if exist('names','var')
        %           uicontrol(f,'style','edit','unit','pixels',...
        %          'position',textpos,...
        %          'backgroundcolor',get(gcf,'color'),...
        %          'string',names{i})
        %      end

             linetable(i).name =   uicontrol(f,'style','edit','unit','pixels',...
                                                           'position',textpos,...
                                                           'backgroundcolor',get(gcf,'color'));
            if exist('names','var')
                set(linetable(i).name, 'string',names{i})
            else
                set(linetable(i).name, 'string', get(linetable(i).handles,'tag')) 
            end
        %      
        end


        guidata(f,linetable)
    otherwise
        error('Unknown option in imagelineplot')
end

function stepthrough(hObject, whatever)

line = findobj(gcf,'tag','waveform');
ax = get(line,'parent');
data = get(ax, 'userdata'); %fields: data, ydata, title
number = get(line,'userdata');

C = abs(get(gcf,'CurrentCharacter'));
number = number - round(28.5-C); %this is a shorthand for the following if statement
%%%%%%%%
% % if C == 29 % right arrow
%     number = number +1;
% elseif C == 28 % left arrow
%     number = number -1;
% end
%%%%%%%%%%

if number > 0 && number < size(data.data,1)
    set(line, 'ydata', data.data(number,:),'userdata',number)
    titlestring = get(data.title,'string');
    %parse the title
    from = strfind(titlestring,'from');
    set(data.title, 'string', [  num2str(data.ydata(number)), ' (Row  ', num2str(number), ')',  from ])
end




