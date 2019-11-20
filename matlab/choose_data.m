function  [cyclefr, figh, hsub]=choose_data(stpfr, mdf, titl, ...
						instruction, xtralines);
%
%  [cyclefr, figh]=detectstepsplot(stepfr, mdf, [titl],...
%                                  [instruction], [xtralines]);
% 
% Plots the time series mdf and adds vertical lines at positions
% given in stpfr. The user is then asked to move and/or modify the
% set of lines.
%
% Input
%    stpfr            ->   Vector with indices
%    mdf              ->   time series. If cell array, then the data
%                          (columns) in each cell is plotted in a subplot. 
%    mdf              ->   time series. If cell array, then the data
%    title            ->   String. title appearing in figure top
%    instruction      ->   String. additional instructions
%    xtralines        ->   Vector of indices giving vertical lines
%                          that are fixed.

% Kjartan Halvorsen
% 2001-02-22
%
% Revisions
% 2004-10-15  Added the possibility to have fixed vertical lines
%             Added use of multiple plot windows

if (nargin<3 | isempty(titl))
  titl = 'Adjust the vertical lines';
end
if (nargin<4)
  instruction = '';
end
if (nargin<5)
  xtralines=[];
end

fig=figure;
figh=fig;


if iscell(mdf)
  nseries = length(mdf);
  margin = 0.1;
  plotmargin = 0.03;
  plotheight = (1 -(2*margin + (nseries-1)*plotmargin))/nseries;
  hsub = zeros(nseries,1);
  for i=nseries:-1:1
    hsub(i) = subplot(nseries, 1, i, ...
		   'Parent', gcf,...
		   'Box','off',...
		   'Position', [margin, ...
		    margin+(nseries-i)*(plotheight+plotmargin) ...   
		    0.77 plotheight]);
    try
      plot(mdf{i});
    catch
      keyboard
    end
    
    xtl = get(gca,'XTickLabel');
    if i<nseries
      set(gca,'XTickLabel',[]);
    end
  end
else
  plot(mdf);
  mdf = {mdf};
end


hold on

if (~isempty(stpfr))
%  keyboard
  cyclelines=plot([stpfr';stpfr'],...
		  [min(min(mdf{1}))*ones(1,length(stpfr))
		   max(max(mdf{1}))*ones(1,length(stpfr))],...
		  'r','EraseMode','xor');
  set(cyclelines,'ZData',[1;1]);
end

if (~isempty(xtralines))
  plot([xtralines';xtralines'],...
       [min(min(mdf{1}))*ones(1,length(xtralines))
	max(max(mdf{1}))*ones(1,length(xtralines))],...
       'm');
end

title(titl, 'Interpreter','none')

axh = gca;

movelinex(fig);

khazoomsafe('on');

if (~isempty(stpfr))
  if (length(stpfr)>1)
    msg={instruction,['Drag the red lines in the plot to adjust ',...
		      'the start of the cycles.'],...
	 ['Add lines by right clicking', ...
	  ' in the figure. Remove lines by holding down ''shift''', ...
	  ' while clicking.'],...
	 ['Zoom by drawing over the region of ',...
	  ' interest. Double-click to zoom back to original.'],...
	 'Click ok when finished'};
  else
    msg={instruction,['Drag the red line in the plot to adjust ',...
		      'the start of the cycle.'],...
	 ['Zoom by drawing over the region of ',...
	  ' interest. Double-click to zoom back to original.'],...
	 'Click ok when finished'};
    
  end	  
else
  msg = {instruction,['Add lines by right clicking', ...
	  ' in the figure. Remove lines by holding down ''shift''', ...
	  ' while clicking.'],['Zoom by drawing over the region of ',...
	  ' interest. Double-click to zoom back to original.'],...
       'Click ok when finished'};
end

msg=textwrap(msg,40);
uiwait(msgbox(msg,'Instructions'));

% Check the XData values of the cyclelines
figure(fig)
lines=get(axh,'Children');

cyclefr=[];
for i=1:length(lines)
  if (get(lines(i),'ZData'))
    xd=get(lines(i),'XData');
    cyclefr=[cyclefr;xd(1)];
  end
end

cyclefr=fix(cyclefr);
cyclefr(find(cyclefr<1))=[];

cyclefr=sort(cyclefr);

% uiwait(msgbox({'Beginning of cycles found at frame #',...
%                sprintf('%d  ',cyclefr)},'Cycles'))
