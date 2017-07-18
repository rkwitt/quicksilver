function boldifyBoxPlot(h,f)
lineWidth = 2;
fontWeight = 16;
titleFontWeight = 18;

for i = 1:numel(h)
  if ~isnan(h(i))
    set(h(i),'LineWidth',lineWidth);
  end
end
set(get(f,'CurrentAxes'),'LineWidth',lineWidth);
set(get(f,'CurrentAxes'),'FontSize',fontWeight);
set(get(f,'CurrentAxes'),'FontWeight','bold');
set(get(get(f,'CurrentAxes'),'YLabel'),'FontSize',fontWeight);
set(get(get(f,'CurrentAxes'),'YLabel'),'FontWeight','bold');
set(get(get(f,'CurrentAxes'),'XLabel'),'FontSize',fontWeight);
set(get(get(f,'CurrentAxes'),'XLabel'),'FontWeight','bold');
set(get(get(f,'CurrentAxes'),'Title'),'FontSize',titleFontWeight);
set(get(get(f,'CurrentAxes'),'Title'),'FontWeight','bold');
