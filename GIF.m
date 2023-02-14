clear
clc

hold off, offset = 0;
subplot('position',[.10 .25 .70 .70]);
LW = 'Linewidth'; lw = 0.7;
MS = 'MarkerSize'; ms = 12;
FS = 'FontSize'; fs = 10;
IN = 'interpret'; LX = 'latex';
L = 12; m = 6; m2 = m/2;
close all
set(gcf , 'Color' , 'w')

% LeftALmostBandMatrix
lbc = fill(offset+[0 2 6 6 0],[6 6 2 0 6]-m,.69*[1 1 1]); set(lbc, 'EdgeColor', 'None'), hold on
plot(offset+[0 0 6 6 0],[0 6 6 0 0]-m,'k',LW,lw)

axis([0 17 -17 0]),axis square, axis off
 
% LeftVector
offset = offset + 6 + 0.5;
plot(offset+[0 0 1 1 0],[0 6 6 0 0]-m,'k',LW,lw)
text(offset+.3,-m2,'$\bf y_{j}$',FS,fs,IN,LX)

% Equals sign
offset = offset + 1 + 0.5;
plot(offset+[.2 .8],(-m2+.1)*[1 1],'k',LW,lw)
plot(offset+[.2 .8],(-m2-.1)*[1 1],'k',LW,lw)

% RightALmostBandMatrix
offset = offset + 1 + 0.5;
rbc = fill(offset+[0 0 2 6 6 5 0],[5 6 6 2 0 0 5]-m,.69*[1 1 1]); set(rbc, 'EdgeColor', 'None')
plot(offset+[0 0 6 6 0],[0 6 6 0 0]-m,'k',LW,lw)

% RightVector
offset = offset + 6 + 0.5;
plot(offset+[0 0 1 1 0],[0 6 6 0 0]-m,'k',LW,lw)
% text(offset+.3,-m2,'$\bf u^{k}$',FS,fs,IN,LX)