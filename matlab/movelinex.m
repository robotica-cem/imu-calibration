function movelinex(varargin)
% movelinex -- move lines in a (x, y) plot along the x-axis by click-and-drag.
%  movelinex('demo') demonstrates itself, using random (x, y)
%   data plotted with EraseMode = 'xor'.
%  movelinex(h) where h is a figure handle enables moving lines in the plot
%   by clicking and dragging.
%   For a line to be possible to move, it must have non-empty, 
%   non-zero z data. 
%   This is to distinguish between lines that can be moved, 
%   and lines that must remain fixed.

% The functionality is implemented using recursive calls to this same function 
% using a second argument specifying the action to take. The different 
% second arguments are 'down', 'move', and 'up'.

% Kjartan Halvorsen
% 2001-01-22 

% Based in code from the function marksafe by Dr. Charles R. Denham, ZYDECO.

LINE_BUTTONDOWNFCN = ['movelinex(gcf,''down'')'];
W_BUTTONMOTIONFCN = ['movelinex(gcf,''move'')'];
W_BUTTONUPFCN = ['movelinex(gcf,''up'')'];
W_BUTTONDWNFCN = ['movelinex(gcf,''windowdown'')'];
LINE_ERASEMODE = 'xor';
LINE_COLOR = 'r';

% Demonstration.

if nargin > 0 & isequal(varargin{1}, 'demo')
   help(mfilename)
   x1=linspace(0, 1, 101);
   y1=sin(3*pi*x1);

   x2=[0.5 0.5];
   y2=[-1 1];
   
   hpl=plot(x1, y1, 'r', x2, y2, 'g', 'EraseMode', 'xor');
   set(hpl(2),'ZData',[1 1]); % Enables moving the vertical line.

   h=gcf;

   set(gcf, 'Name', [mfilename ' demo'])

   movelinex(h);

   return
end

%if isempty(gcbo)
if (nargin<2)
   if isempty(varargin)
      set(gcf,'WindowButtonDownFcn',W_BUTTONDWNFCN);
      lh = findobj(gcf, 'Type', 'line');
   else
      set(varargin{1},'WindowButtonDownFcn',W_BUTTONDWNFCN);
      lh = findobj(varargin{1}, 'Type', 'line');
   end

   set(lh, 'ButtonDownFcn', LINE_BUTTONDOWNFCN)
end

% Process mouse events.

if (nargin==2)
   h = varargin{1};
   action=varargin{2};

   switch lower(action)
      case 'down' 
         % Process a down click. 
         % If the line can be moved at all (controlled by having non-empty
         % and non-zero ZData.
         % If a) left-button-click (the figures SelectionType is 'normal') 
         % then the WindowButtonMotionFcn and the           
         % WindowButtonUpFcns of the current figure is set,
         % If b) middle-button-click or for windows, clicking both 
         % buttons simultanously, alternatively, holding the shift 
         % key while left-clicking 
         % (in all cases the figures SelectionType is 'extended') 
         % then the line is deleted.

         if (get(gcbo,'ZData'))
            if (strcmp(get(h,'SelectionType'),'normal'))
               set(h,'WindowButtonMotionFcn',W_BUTTONMOTIONFCN, ...
                     'WindowButtonUpFcn',W_BUTTONUPFCN)
               set(gcbo,'EraseMode','xor')
               orig.x=get(gcbo,'XData');
               orig.p=get(h,'CurrentPoint');
               set(gcbo,'UserData',orig);
               set(h,'UserData',gcbo);
            elseif (strcmp(get(h,'SelectionType'),'extend'))
               set(gcbo,'XData',[]);
               set(gcbo,'YData',[]);
               set(gcbo,'ZData',[]);
               drawnow
            end %if
         end
      case 'move'
         % Read the current point; if inside the axes then 
         % plot the line at the new position         

         cp=get(h,'CurrentPoint');

         xl=get(gca,'XLim');
         xrange = diff(xl);
   
         theOldUnits = get(gca, 'Units');
         set(gca, 'Units', 'pixels');
         pos = get(gca, 'Position');
         theWidth = pos(3);
         set(gca, 'Units', theOldUnits)

         if (cp(1)<pos(1) | cp(1)>pos(1)+pos(3) | cp(2)<pos(2) | cp(2)>pos(2)+pos(4))
            return
         end


         lh=get(h,'UserData'); 
         orig=get(lh,'UserData'); % The original 

         cp=get(h,'CurrentPoint');

         % Scale to correct units.
         dx = (cp(1)-orig.p(1))*xrange/theWidth;

         set(lh,'XData', orig.x+dx);

         drawnow

      case 'up'
         % Reset the WindowButtinMotionFcn to the empty string to stop the
         % callback.

         set(h,'WindowButtonMotionFcn','');
         set(h,'UserData',[]);


      case 'windowdown'
         % The user has clicked somewhere on the graph outside the moveable
         % lines. If the right button was clicked, then create a new
         % vertical line at the current point.

         if (strcmp(get(h,'SelectionType'),'alt'))
   
            cp=get(h,'CurrentPoint');

            xl=get(gca,'XLim');
            xrange = diff(xl);
            yl=get(gca,'YLim');
            yrange = diff(yl);
   
            theOldUnits = get(gca, 'Units');
            set(gca, 'Units', 'pixels');
            pos = get(gca, 'Position');
            theWidth = pos(3);
            theHeight = pos(4);
            set(gca, 'Units', theOldUnits);

            if (cp(1)<pos(1) | cp(1)>pos(1)+pos(3) | cp(2)<pos(2) | cp(2)>pos(2)+pos(4))
               return
            end


            % Scale to correct units.
            x = xl(1)+(cp(1)-pos(1))*xrange/theWidth*ones(2,1);
            y = yl;
            z = ones(2,1);

            set(h,'NextPlot','add'); % Same as 'hold on'   
            set(gca,'NextPlot','add'); % Same as 'hold on'   
            lh=plot(x,y,LINE_COLOR);
  
            set(lh,'ZData', z);

            set(lh, 'ButtonDownFcn', LINE_BUTTONDOWNFCN);

         else
            disp('Window down left button')
         end
            
   end
end
