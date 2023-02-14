clear
clc
%%
% A1显示方法以及RK

hold off, offset = 0;
subplot('position',[.10 .25 .70 .70]);
LW = 'Linewidth'; lw = 0.7;
MS = 'MarkerSize'; ms = 12;
FS = 'FontSize'; fs = 5;
IN = 'interpret'; LX = 'latex';
L = 12; m = 6; m2 = m/2;
close all
set(gcf , 'Color' , 'w')

% LeftALmostBandMatrix
lbc = fill(offset+[0 0 6 6 0],[5 6 6 5 5]-m,.69*[1 1 1]); set(lbc, 'EdgeColor', 'None'), hold on
rbc = fill(offset+[0 0 2 6 6 5 0],[5 5 5 1 0 0 5]-m,.69*[1 1 1]); set(rbc, 'EdgeColor', 'None')
plot(offset+[0 0 6 6 0],[0 6 6 0 0]-m,'k',LW,lw)

axis([0 30 -30 0]),axis square, axis off

% LeftVector
offset = offset + 6 + 0.5;
plot(offset+[0 0 1 1 0],[0 6 6 0 0]-m,'k',LW,lw)
text(offset+.2,-m2,'$\bf u^{k+r}$',FS,fs,IN,LX)

% Equals sign
offset = offset + 1 + 0.5;
plot(offset+[.2 .8],(-m2+.1)*[1 1],'k',LW,lw)
plot(offset+[.2 .8],(-m2-.1)*[1 1],'k',LW,lw)

% RightALmostBandMatrix
offset = offset + 1 + 0.5;
lbc = fill(offset+[0 0 6 6 0],[5 6 6 5 5]-m,.69*[1 1 1]); set(lbc, 'EdgeColor', 'None'), hold on
rbc = fill(offset+[0 0 2 6 6 4 0],[4 5 5 1 0 0 4]-m,.69*[1 1 1]); set(rbc, 'EdgeColor', 'None')
plot(offset+[0 0 6 6 0],[0 6 6 0 0]-m,'k',LW,lw)

% RightVector
offset = offset + 6 + 0.5;
plot(offset+[0 0 1 1 0],[0 6 6 0 0]-m,'k',LW,lw)
text(offset+.1,-m2,'$\bf u^{k+r-1}$',FS,fs,IN,LX)

% RightPlus
offset = offset + 1 + 0.5;
plot(offset+[.2 .8],(-m2)*[1 1],'k',LW,lw)
plot(offset+[.5 .5],[-m2+0.3 -m2-0.3],'k',LW,lw)

% ...
offset = offset + 1 + 0.5;
plot(offset+[.2 0.5 .8],(-m2)*[1 1 1],'.k',LW,lw)

% RightPlus
offset = offset + 1 + 0.5;
plot(offset+[.2 .8],(-m2)*[1 1],'k',LW,lw)
plot(offset+[.5 .5],[-m2+0.3 -m2-0.3],'k',LW,lw)

% RightALmostBandMatrix
offset = offset + 1 + 0.5;
lbc = fill(offset+[0 0 6 6 0],[5 6 6 5 5]-m,.69*[1 1 1]); set(lbc, 'EdgeColor', 'None'), hold on
rbc = fill(offset+[0 0 2 6 6 4 0],[4 5 5 1 0 0 4]-m,.69*[1 1 1]); set(rbc, 'EdgeColor', 'None')
plot(offset+[0 0 6 6 0],[0 6 6 0 0]-m,'k',LW,lw)

% RightVector
offset = offset + 6 + 0.5;
plot(offset+[0 0 1 1 0],[0 6 6 0 0]-m,'k',LW,lw)
text(offset+.25,-m2,'$\bf u^{k}$',FS,fs,IN,LX)

