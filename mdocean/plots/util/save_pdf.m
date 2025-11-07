function save_pdf(fig,filename)
%SAVE_PDF Save figure as pdf to specified filename

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,filename,'-dpdf','-r0')

end
