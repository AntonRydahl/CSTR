function fig_2_pdf()
folder = pwd;
filePattern = fullfile(folder, '*.fig');
figFiles = dir(filePattern);
for k = 1:length(figFiles)
  figName = figFiles(k).name;
  filePath = fullfile(folder, figName);
  fig = openfig(filePath,'invisible');
  pdfName = replace(figName,'.fig','');
  set(gcf,'Units','inches');
  screenposition = get(gcf,'Position');
  set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
  print(pdfName,'-dpdf')
end