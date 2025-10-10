function varargout = plotUnc(varargin)
%% Purpose:
%
%  This routine will plot the 2D uncertainty bounds about a standard plot
%  of x and y data.  The bounds are assumed to be +/- of the specified
%  uncertainty value
%
%  Example Call Sequences:
%
%  plotUnc(x,y,u,l); 
%  plotUnc(y,u,l);
%  plotUnc(h,x,y,u,l);
%
%  MATLAB Command Line Example:
%
%  plotUnc(sin(0:0.01:2*pi).*10, (0:0.01:2*pi).*0 + 2, ...
%  (0:0.01:2*pi).*0 + 2,'FaceAlpha',0.1,'EdgeColor','none');
%
%  Where 
%
%  x and y are similar to the plot.m routine,
%  u represents the upper bound (above y) and l represents the lower bound
%  (below y) for the shaded region (uncertainty).
%
%  User input properties can be specified either for the built-in plot 
%  routine or the built-in fill routine.
%
%
%% Revision History:
%  Darin C. Koblick                                         (c) 08-15-2022
%% ----------------------- Begin Code Sequence ----------------------------
if nargin == 0
     x = [1:0.01:6*pi; 1:0.01:6*pi; 1:0.01:6*pi]';
     y = [sin(x(:,1)).*3 + 20, cos(x(:,2)).*3 + 0,  cos(x(:,3)).*sin(x(:,3)).*3 + 10];
     u = ones(size(x)) + abs(normrnd(0,0.5,size(x)));
     l = u;
   figure('color',[0 0 0]);
   subplot(2,2,1);
   plotUnc(y(:,1),u(:,1),l(:,1),'r','FaceColor','none','EdgeAlpha',0.4);
   plotUnc(y(:,2),u(:,2),l(:,2),'g','FaceColor','none','EdgeAlpha',0.4);
   plotUnc(y(:,3),u(:,3),l(:,3),'b','FaceColor','none','EdgeAlpha',0.4);
   %legend('\color{white}y_1','\color{white}y_2','\color{white}y_3', ...
   %       'Location','NW','color','none','AutoUpdate','off');
   set(gca,'color','k','xColor','w','yColor','w');
   grid on;
   title('\color{white}plotUnc(y,u,l,cSpec)');
   subplot(2,2,2);
   plotUnc(y,u,l,'FaceAlpha',0.3,'EdgeAlpha',0.2);
   set(gca,'color','k','xColor','w','yColor','w');
   grid on;
   title('\color{white}plotUnc(y,u,l)');
   subplot(2,2,3);
   yyaxis left;
   plotUnc(x(:,1),y(:,1),u(:,1),l(:,1),'w','Edgecolor','w','FaceAlpha',0.4);
   ylim([-10 30]);
   yyaxis right;
   plotUnc(x(:,2:end),y(:,2:end),u(:,2:end),l(:,2:end),'b','EdgeColor','w','FaceAlpha',1);
   ylim([-10 30]);
   set(gca,'color','k','xColor','w','yColor','w');
   grid on;
   title('\color{white}plotUnc(x,y,u,l,cSpec)');
   hs = subplot(2,2,4);
   ha = plotUnc(hs,y,u,l,'y','Edgecolor','none','FaceAlpha',1);
   for tc = 1:size(ha.Children,1)
       ha.Children(tc).AmbientStrength = 0.4;
       ha.Children(tc).DiffuseStrength = 0.0;
       ha.Children(tc).SpecularStrength = 0.8;
       ha.Children(tc).SpecularExponent = 1000;
       ha.Children(tc).FaceLighting = 'gouraud';
       ha.Children(tc).BackFaceLighting  = 'unlit';
       ha.Children(tc).EdgeLighting = 'gouraud';
   end
   camlight(ha,'headlight');
   lightangle(ha,0,30);
   set(gca,'color','k','xColor','w','yColor','w');
   grid on;
   title('\color{white}plotUnc(h,y,u,l,opts\{:\})');

   return;
elseif nargin < 3
   error('Missing  y, u, and l inputs');
end
%4 Possible input sequences:
%plotUnc(h,x,y,u,l, ...)
%plotUnc(h,y,u,l, ...)
%plotUnc(x,y,u,l, ...)
%plotUnc(y,u,l, ...)
if nargin < 5
    varargin = [varargin, repmat({[]},1,5-nargin)];
end
idxRem = false(size(varargin));
%Fifth input can be either empty, a colorSpec, l, or another argument
if isempty(varargin{5})
   idxRem(5) = true;
elseif isnumeric(varargin{5}) && isnumeric(varargin{4})
           l = varargin{5};
           u = varargin{4};
 idxRem(4:5) = true;
elseif isColorSpec(varargin{5}) && isnumeric(varargin{4})
           l = varargin{4};
           u = varargin{3};
 idxRem(3:4) = true;       
end
if all(~idxRem(1:4))
    %Forth input can be either empty, a colorSpec, u, l, or another argument
    if isempty(varargin{4})
        idxRem(4) = true;
    elseif isnumeric(varargin{4})
               l = varargin{4};
               u = varargin{3};
     idxRem(3:4) = true;
    elseif isColorSpec(varargin{4}) 
               l = varargin{3};
               u = varargin{2};
     idxRem(2:3) = true;
    end
end
if all(~idxRem(1:3))
   %Third input can be either y, u, or l 
             l = varargin{3};
             u = varargin{2};
   idxRem(2:3) = true;
end
%Remove excess arguments from input:
varargin(idxRem) = [];
%Grab the property Fields of the plotting function:
            hp = plot(NaN,NaN);
    plotfNames = fieldnames(get(hp));
    delete(hp);
            hf = fill(NaN,NaN,NaN);
    fillfNames = fieldnames(get(hf));
    delete(hf);
  vararginPlot = extractVarargins(varargin(:),plotfNames);
%Plot all arguments (except uncertainty):
           h = plot(vararginPlot{:});
          ha = h.Parent;
          hold(ha,'on');
       cSpec = {h.Color}';
      nLines = size(cSpec,1);
for tl=nLines:-1:1
          x(:,tl) = permute(h(tl).XData,[2 1]);
          y(:,tl) = permute(h(tl).YData,[2 1]);
end
vararginFill = extractVarargins(varargin(:),fillfNames);
%Remove x and/or y from the variable input arguments:
count = 1;
for i=1:size(vararginFill,1)
    if ~isnumeric(vararginFill{i}) && ~ishandle(vararginFill{i}(1))
        break;
    else
       count = count+1; 
    end
end
%Add an empty array allocation for the color code, x, and y axis:
vararginFill = [{ha}; x; y; [0 0 0]; vararginFill(count:end)];
idxEdgeColor = strcmpi(vararginFill,'EdgeColor');
if any(idxEdgeColor)
   edgeColorLoc = NaN;
else
   vararginFill = [vararginFill; 'EdgeColor'; [0 0 0]];
   edgeColorLoc = size(vararginFill,1);
end
%Add the Uncertainty as a shaded region over each line:
for tl=1:nLines
       %Assume second and third arguments are for x and y vectors:
            vararginFill{2} = [x(:,tl); flipud(x(:,tl))];
            vararginFill{3} = [y(:,tl)-abs(l(:,tl)); flipud(y(:,tl)+abs(u(:,tl)))];
            vararginFill{4} = cSpec{tl};
            if ~isnan(edgeColorLoc)
                vararginFill{edgeColorLoc} = cSpec{tl};
            end
       %Add the color specification code:
       patch(vararginFill{:});
end
%Take the lines used for the color and remove them:
for tl=size(ha.Children,1):-1:1
    idx(tl,1) = strcmpi(ha.Children(tl).Type,'Line');
end
delete(ha.Children(idx));
if nargout > 0
   varargout{1} = ha; 
end
end

function outputArgs = extractVarargins(inputArgs,refArgs)
%Find those arguments that intersect:
outputArgs = inputArgs;
    idxRem = false(size(outputArgs,1),1);
for count=1:size(inputArgs,1)
    if ~isnumeric(inputArgs{count}) && ~ishandle(inputArgs{count}(1))
       break; 
    end
end

%Determine if the next argument is a colorSpec:
if isColorSpec(inputArgs{count}) && any(strcmpi(refArgs,'Color'))
   count = count + 1;
elseif isColorSpec(inputArgs{count}) && ~any(strcmpi(refArgs,'Color'))
   %Remove the colorSpec argument:
   outputArgs(count) = [];
              idxRem = false(size(outputArgs,1),1);
end

while count < size(outputArgs,1)
    
    if any(strcmpi(outputArgs{count},refArgs))
       %Keep the Argument identifer and the value from the list:
         idxRem(count) = false;
       idxRem(count+1) = false;
                 count = count+2;
    else
       %Remove from the argument identifier and value from the list:
         idxRem(count) = true;
       idxRem(count+1) = true;
                 count = count+2;
    end 
end
outputArgs(idxRem) = [];
end

function output = isColorSpec(input)
shortName = {'r','g','b','c','m','y','k','w'};
 fullName = {'red','green','blue','cyan','magenta','yellow','black','white'};
   output = false;
if iscell(input)
    if all(cell2mat(cellfun(@(x)any(strcmpi(x,fullName)),input,'UniformOutput',false))) || ...
       all(cell2mat(cellfun(@(x)any(strcmpi(x,shortName)),input,'UniformOutput',false)))
        output = true;
    end
else
    if any(strcmpi(input,fullName)) || any(strcmpi(input,shortName))
           output = true; 
    end
end
end