function boldifyPlot(h,f,lineWidth,fontWeight,titleFontWeight)
%
% use like this:
% 
%

if nargin < 5
  titleFontWeight = 14;
end

if nargin < 4
  fontWeight = 12;
end

if nargin < 3
  lineWidth = 2;
end

for i = 1:numel(f)
  if ~isnan(f(i))
    set(f(i),'LineWidth',lineWidth);
  end
end
set(get(h,'CurrentAxes'),'LineWidth',lineWidth);
set(get(h,'CurrentAxes'),'FontSize',fontWeight);
set(get(h,'CurrentAxes'),'FontWeight','bold');
set(get(get(h,'CurrentAxes'),'YLabel'),'FontSize',fontWeight);
set(get(get(h,'CurrentAxes'),'YLabel'),'FontWeight','bold');
set(get(get(h,'CurrentAxes'),'XLabel'),'FontSize',fontWeight);
set(get(get(h,'CurrentAxes'),'XLabel'),'FontWeight','bold');
set(get(get(h,'CurrentAxes'),'Title'),'FontSize',titleFontWeight);
set(get(get(h,'CurrentAxes'),'Title'),'FontWeight','bold');
