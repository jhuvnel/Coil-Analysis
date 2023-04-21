%%
RL = load('R:\Morris, Brian\Monkey Data\Nancy\20191118\Cycles\Nancy-NA-20191118-ElectricalStim-PulseTrains-LARP-stim6ref5-200pps-phase1Dur200-phase2Dur200-IPG0-phase1Amp250-phase2Amp250_CycleAvg.mat');
RZ = load('R:\Morris, Brian\Monkey Data\Nancy\20191118\Cycles\Nancy-NA-20191118-ElectricalStim-PulseTrains-LHRH-stim7ref15-200pps-phase1Dur200-phase2Dur200-IPG0-phase1Amp250-phase2Amp250_CycleAvg.mat');
RR = load('R:\Morris, Brian\Monkey Data\Nancy\20191118\Cycles\Nancy-NA-20191118-ElectricalStim-PulseTrains-RALP-stim3ref2-200pps-phase1Dur200-phase2Dur200-IPG0-phase1Amp250-phase2Amp250_CycleAvg.mat');
%%
RL = load('R:\Morris, Brian\Monkey Data\Opal\20220617\Cycles\Opal-NA-20220617-ElectricalStim-PulseTrains-LARP-stim15ref11-200pps-phase1Dur200-phase2Dur200-IPG0-phase1Amp250-phase2Amp250_CycleAvg.mat');
RZ = load('R:\Morris, Brian\Monkey Data\Opal\20220617\Cycles\Opal-NA-20220617-ElectricalStim-PulseTrains-LHRH-stim6ref11-200pps-phase1Dur200-phase2Dur200-IPG0-phase1Amp250-phase2Amp250_CycleAvg.mat');
RR = load('R:\Morris, Brian\Monkey Data\Opal\20220617\Cycles\Opal-NA-20220617-ElectricalStim-PulseTrains-RALP-stim3ref11-200pps-phase1Dur200-phase2Dur200-IPG0-phase1Amp250-phase2Amp250_CycleAvg.mat');
%%
LL = load('R:\Morris, Brian\Monkey Data\GiGi\20191206\cycles\GiGi-NA-20191206-ElectricalStim-PulseTrains-LARP-stim5ref11-200pps-phase1Dur200-phase2Dur200-IPG0-phase1Amp150-phase2Amp150_CycleAvg.mat');
LR = load('R:\Morris, Brian\Monkey Data\GiGi\20191206\cycles\GiGi-NA-20191206-ElectricalStim-PulseTrains-RALP-stim3ref11-200pps-phase1Dur200-phase2Dur200-IPG0-phase1Amp250-phase2Amp250_CycleAvg.mat');
LZ = load('R:\Morris, Brian\Monkey Data\GiGi\20191206\cycles\GiGi-NA-20191206-ElectricalStim-PulseTrains-LHRH-stim8ref11-200pps-phase1Dur200-phase2Dur200-IPG0-phase1Amp250-phase2Amp250_CycleAvg.mat');
%%
H_Gigi_ringXYZ_re_G = ... % note the apostrophe to transpose rows to columns; re_G means re image system coords
 [-0.207131	 0.978284	-0.00756488; % +X yellow Gigi head ring X (unlike Nancy CT)
   -0.967023	 -0.205906	 -0.149897;   %+Y (neated from A green in ARS coords) for Gigi head ring Y
  -0.148199	-0.0237328	 0.988673]'; % red for Gigi head ring Z
%H = H_Gigi_ringXYZ_re_G
C_Gigi_leftcanalLRZ_re_G = ... % note the apostrophe to transpose rows to columns; re_G mean re image system coords
 [ 0.676684	0.734721	 0.0477903; % green +LA
  -0.726268	0.676745	-0.120628;  % yellow +LP
  -0.12097	0.0469186	 0.991547]';% red +Z
%C=C_Gigi_leftcanalLRZ_re_G

H_Nancy_ringXYZ_re_G = ... % note the apostrophe to transpose rows to columns; re_G means re image system gantry coords
  [-0.0196192	0.999762	-0.00953597;
   -0.999273	 -0.0192961	 0.03287;
    0.0326782	0.0101739	 0.999414]'; % transposed for consistency with C:\Users\cdell\OneDrive - Johns Hopkins\VNEL1DRV\Manuscripts-VNEL1DRV\Morris1-StimRefPosn\  3DSlicerAnalaysisForeeVORCanalVsFieldCoilAlignment.xlsx

C_Nancy_leftcanalLRZ_re_G = ... % note the apostrophe to transpose rows to columns; re_G means re image system gantry coords
  [ 0.73676	  0.608407	-0.295002;
   -0.664495	0.732182	-0.149518;
    0.125028	0.306186	 0.943726]';

H_Opal_ringXYZ_re_G = ... % note the apostrophe to transpose rows to columns; re_G means re image system gantry coords
  [-0.950912	 0.294825   0.0940482;
   -0.2935    -0.955551   0.027941;
    0.0981056	-0.0010337	0.995175]';
  
C_Opal_leftcanalLRZ_re_G = ... % note the apostrophe to transpose rows to columns; re_G means re image system gantry coords  
  [ -0.573072	 0.818776	-0.0345426;
    -0.809681	-0.559193	 0.178101;
     0.126509	 0.130033	 0.983406]';
 
%P_CLRZfromFXYZ = inv(C)*H
%P_CLRZfromFXYZ = C\H % x = A\b is computed differently than x = inv(A)*b and is recommended for solving systems of linear equations
%P_CLRZfromFLRZ = P_CLRZfromFXYZ*rotZ3deg(-45)

P_CLRZfromFXYZ_Gigi = C_Gigi_leftcanalLRZ_re_G\H_Gigi_ringXYZ_re_G; 
P_CLRZfromFLRZ_Gigi = P_CLRZfromFXYZ_Gigi*rotZ3deg(-45);

P_CLRZfromFXYZ_Nancy = C_Nancy_leftcanalLRZ_re_G\H_Nancy_ringXYZ_re_G;
P_CLRZfromFLRZ_Nancy = P_CLRZfromFXYZ_Nancy*rotZ3deg(-45);

P_CLRZfromFXYZ_Opal = C_Opal_leftcanalLRZ_re_G\H_Opal_ringXYZ_re_G; 
P_CLRZfromFLRZ_Opal = P_CLRZfromFXYZ_Opal*rotZ3deg(-45);
%%
x =filtfilt(ones(1,7)/7,1,RR.Results.segmentData.RE_Velocity_X);
y =filtfilt(ones(1,7)/7,1,RR.Results.segmentData.RE_Velocity_Y);
z =filtfilt(ones(1,7)/7,1,RR.Results.segmentData.RE_Velocity_Z);
l =filtfilt(ones(1,7)/7,1,RR.Results.segmentData.RE_Velocity_LARP);
r =filtfilt(ones(1,7)/7,1,RR.Results.segmentData.RE_Velocity_RALP);
lrz = [l r z];
to = -45;
rot = [cosd(to) -sind(to) 0; sind(to) cosd(to) 0; 0 0 1]*lrz';
figure
plot(x,'k','LineWidth',4)
hold on
plot(y,'color',[0.4 0.4 0.4'],'LineWidth',4)
plot(z,'color',[0.7 0.7 0.7'],'LineWidth',4)
plot(l,'color',[0 0 1],'LineWidth',2)
plot(r,'color',[0 0 .5],'LineWidth',2)
plot(rot(1,:),'color',[1 0 0],'LineWidth',.5)
plot(rot(2,:),'color','g','LineWidth',.5)
plot(rot(3,:),'color',[.4 0 0],'LineWidth',.5)
hold off
%%
% 1. Pull Raw LARP, RALP, Z Velocity data and filter using original
% parameters
% 2. Pull values at positions used previously
% 3. Calculate misalignment values for used cycles
% 4. Apply 45 degree Z rotation to cycle velocity values LRZ->XYZ
% 5. Apply desired rotation about Y axis
% 6. Apply -45 degree rotation about Z XYZ-> LRZ
% 7. Recalculate misalignment values

toRot = -16;

% 1.
% Larp = filtfilt(ones(1,3)/3,1,RL.Results.segmentData.RE_Velocity_LARP);
% Ralp = filtfilt(ones(1,3)/3,1,RL.Results.segmentData.RE_Velocity_RALP);
% Z = filtfilt(ones(1,3)/3,1,RL.Results.segmentData.RE_Velocity_Z);

% 2.
rpr = repmat([1 0 0],length(RL.Results.usedpullIndsR_NystagCorr),1);
% LRZ = [Larp(RL.Results.usedpullIndsR) Ralp(RL.Results.usedpullIndsR) Z(RL.Results.usedpullIndsR)];
LRZL = RL.Results.usedMisalign3DR_NystagCorr;
% 3.
misVals = acosd(dot(rpr,LRZL,2)./(vecnorm(rpr,2,2).*vecnorm(LRZL,2,2)));
[misVals RL.Results.usedMisalignR_NystagCorr'];

% 4.
XYZVL = [cosd(-45) -sind(-45) 0; sind(-45) cosd(-45) 0; 0 0 1]*LRZL';

% 5.

XYZRVL = [cosd(toRot) 0 sind(toRot); 0 1 0; -sind(toRot) 0 cosd(toRot)]*XYZVL;

% 6.
RLLRZV = [cosd(45) -sind(45) 0; sind(45) cosd(45) 0; 0 0 1]*XYZRVL;

% 7.
misVals2 = acosd(dot(rpr,RLLRZV',2)./(vecnorm(rpr,2,2).*vecnorm(RLLRZV',2,2)));
% misVals and misVals2 are equal if toRot = 0
[misVals2 misVals RL.Results.usedMisalignR_NystagCorr'];
Lold = mean(misVals);
Lnew = mean(misVals2);


% 1. Pull Raw LARP, RALP, Z Velocity data and filter using original
% parameters
% 2. Pull values at positions used previously
% 3. Calculate misalignment values for used cycles
% 4. Apply 45 degree Z rotation to cycle velocity values LRZ->XYZ
% 5. Apply desired rotation about Y axis
% 6. Apply -45 degree rotation about Z XYZ-> LRZ
% 7. Recalculate misalignment values

% 1.
% Larp = filtfilt(ones(1,3)/3,1,RR.Results.segmentData.RE_Velocity_LARP);
% Ralp = filtfilt(ones(1,3)/3,1,RR.Results.segmentData.RE_Velocity_RALP);
% Z = filtfilt(ones(1,3)/3,1,RR.Results.segmentData.RE_Velocity_Z);

% 2.
rpr = repmat([0 1 0],length(RR.Results.usedpullIndsR_NystagCorr),1);
% LRZ = [Larp(RR.Results.usedpullIndsR) Ralp(RR.Results.usedpullIndsR) Z(RR.Results.usedpullIndsR)];
LRZR = RR.Results.usedMisalign3DR_NystagCorr;
% 3.
misVals = acosd(dot(rpr,LRZR,2)./(vecnorm(rpr,2,2).*vecnorm(LRZR,2,2)));
[misVals RR.Results.usedMisalignR_NystagCorr'];

toRot = -16;
r1a = [1 0 0; 0 cosd(toRot) -sind(toRot); 0 sind(toRot) cosd(toRot)]*LRZR';
r2a = [cosd(-toRot) 0 sind(-toRot); 0 1 0; -sind(-toRot) 0 cosd(-toRot)]*r1a;

% 4.
XYZVR = [cosd(-45) -sind(-45) 0; sind(-45) cosd(-45) 0; 0 0 1]*LRZR';
% 5.
XYZRVR = [cosd(toRot) 0 sind(toRot); 0 1 0; -sind(toRot) 0 cosd(toRot)]*XYZVR;

% 6.
RRLRZV = [cosd(45) -sind(45) 0; sind(45) cosd(45) 0; 0 0 1]*XYZRVR;

% 7.
misVals2 = acosd(dot(rpr,RRLRZV',2)./(vecnorm(rpr,2,2).*vecnorm(RRLRZV',2,2)));
[misVals2 RR.Results.usedMisalignR_NystagCorr'];
Rold = mean(misVals);
Rnew = mean(misVals2);

% misVals and misVals2 are equal if toRot = 0

% 1. Pull Raw LARP, RALP, Z Velocity data and filter using original
% parameters
% 2. Pull values at positions used previously
% 3. Calculate misalignment values for used cycles
% 4. Apply 45 degree Z rotation to cycle velocity values LRZ->XYZ
% 5. Apply desired rotation about Y axis
% 6. Apply -45 degree rotation about Z XYZ-> LRZ
% 7. Recalculate misalignment values

% 1.
% Larp = filtfilt(ones(1,3)/3,1,RR.Results.segmentData.RE_Velocity_LARP);
% Ralp = filtfilt(ones(1,3)/3,1,RR.Results.segmentData.RE_Velocity_RALP);
% Z = filtfilt(ones(1,3)/3,1,RR.Results.segmentData.RE_Velocity_Z);

% 2.
rpr = repmat([0 0 -1],length(RZ.Results.usedpullIndsR_NystagCorr),1);
% LRZ = [Larp(RR.Results.usedpullIndsR) Ralp(RR.Results.usedpullIndsR) Z(RR.Results.usedpullIndsR)];
LRZZ = RZ.Results.usedMisalign3DR_NystagCorr;
% 3.
misVals = acosd(dot(rpr,LRZZ,2)./(vecnorm(rpr,2,2).*vecnorm(LRZZ,2,2)));
[misVals RZ.Results.usedMisalignR_NystagCorr'];

% 4.
XYZVZ = [cosd(-45) -sind(-45) 0; sind(-45) cosd(-45) 0; 0 0 1]*LRZZ';
% 5.
XYZRVZ = [cosd(toRot) 0 sind(toRot); 0 1 0; -sind(toRot) 0 cosd(toRot)]*XYZVZ;

% 6.
RZLRZV = [cosd(45) -sind(45) 0; sind(45) cosd(45) 0; 0 0 1]*XYZRVZ;

% 7.
misVals2 = acosd(dot(rpr,RZLRZV',2)./(vecnorm(rpr,2,2).*vecnorm(RZLRZV',2,2)));
[misVals2 RZ.Results.usedMisalignR_NystagCorr'];
Zold = mean(misVals);
Znew = mean(misVals2);

% misVals and misVals2 are equal if toRot = 0
[Lold Lnew]
[Rold Rnew]
[Zold Znew]
%%
avgMisalignPlot3D = figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle(avgMisalignPlot3D,{'3D Angle of Misalignment'},'FontSize', 22, 'FontWeight', 'Bold');

avgMisalign3D = axes('Parent', avgMisalignPlot3D);
avgMisalign3D.XGrid = 'on';
avgMisalign3D.YGrid = 'on';
avgMisalign3D.FontSize = 13.5;

hold(avgMisalign3D,'on');
h=plot3vect([1;0;0],'LARP Axis',[0 1 0],2);
set(h,'LineStyle','--','Marker','o');
h=plot3vect([0;1;0],'RALP Axis',[0 0 1],2);
set(h,'LineStyle','--','Marker','o');
h=plot3vect([0;0;-1],'Yaw Axis',[1 0 0],2);
set(h,'LineStyle','--','Marker','o');


[x,y,z]=sphere();
h=surf(0.5*x,0.5*y,0.5*z);
set(h,'FaceColor','white')
avgMisalign3D.View = [135 -15];
axis vis3d
axis equal
box on;
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
%%
VN = vecnorm(LRZR,2,2);
for i = 1:length(RR.Results.usedpullIndsR_NystagCorr)
    lPlotLD = plot3([0 LRZR(i,1)/VN(i)]',[0 LRZR(i,2)/VN(i)]',[0 LRZR(i,3)/VN(i)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','k')
end


% VN = vecnorm(LRZR,2,2);
% for i = 1:length(RR.Results.usedpullIndsR_NystagCorr)
%     qq=rotate3DpointAround3DAxisByThetaDEG(LRZR(i,:),[-1 1 0],-16);
%     lPlotLD = plot3([0 qq(1)/VN(i)]',[0 qq(2)/VN(i)]',[0 qq(3)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','g')
% end

% 
% VN = vecnorm(XYZVR',2,2);
% for i = 1:length(RR.Results.usedpullIndsR_NystagCorr)
%     lPlotLD = plot3([0 XYZVR(1,i)/VN(i)]',[0 XYZVR(2,i)/VN(i)]',[0 XYZVR(3,i)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','y')
% end
% 
% VN = vecnorm(XYZRVR',2,2);
% for i = 1:length(RR.Results.usedpullIndsR_NystagCorr)
%     lPlotLD = plot3([0 XYZRVR(1,i)/VN(i)]',[0 XYZRVR(2,i)/VN(i)]',[0 XYZRVR(3,i)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','b')
% end
% 
VN = vecnorm(RRLRZV',2,2);
for i = 1:length(RR.Results.usedpullIndsR_NystagCorr)
    lPlotLD = plot3([0 RRLRZV(1,i)/VN(i)]',[0 RRLRZV(2,i)/VN(i)]',[0 RRLRZV(3,i)/VN(i)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','g')
end

VN = vecnorm(LRZL,2,2);
for i = 1:length(RL.Results.usedpullIndsR_NystagCorr)
    lPlotLD = plot3([0 LRZL(i,1)/VN(i)]',[0 LRZL(i,2)/VN(i)]',[0 LRZL(i,3)/VN(i)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','k')
end

% VN = vecnorm(LRZL,2,2);
% for i = 1:length(RL.Results.usedpullIndsR_NystagCorr)
%     qq=rotate3DpointAround3DAxisByThetaDEG(LRZL(i,:),[-1 1 0],-16);
%     lPlotLD = plot3([0 qq(1)/VN(i)]',[0 qq(2)/VN(i)]',[0 qq(3)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','g')
% end

% VN = vecnorm(r2',2,2);
% for i = 1:length(RL.Results.usedpullIndsR_NystagCorr)
%     lPlotLD = plot3([0 r2(1,i)/VN(i)]',[0 r2(2,i)/VN(i)]',[0 r2(3,i)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','y')
% end
% 
% VN = vecnorm(XYZRVL',2,2);
% for i = 1:length(RL.Results.usedpullIndsR_NystagCorr)
%     lPlotLD = plot3([0 XYZRVL(1,i)/VN(i)]',[0 XYZRVL(2,i)/VN(i)]',[0 XYZRVL(3,i)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','b')
% end
% 
VN = vecnorm(RLLRZV',2,2);
for i = 1:length(RL.Results.usedpullIndsR_NystagCorr)
    lPlotLD = plot3([0 RLLRZV(1,i)/VN(i)]',[0 RLLRZV(2,i)/VN(i)]',[0 RLLRZV(3,i)/VN(i)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','g')
end
% 
% 
VN = vecnorm(LRZZ,2,2);
for i = 1:length(RZ.Results.usedpullIndsR_NystagCorr)
    lPlotLD = plot3([0 LRZZ(i,1)/VN(i)]',[0 LRZZ(i,2)/VN(i)]',[0 LRZZ(i,3)/VN(i)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','k')
end
% 
% VN = vecnorm(XYZVZ',2,2);
% for i = 1:length(RZ.Results.usedpullIndsR_NystagCorr)
%     lPlotLD = plot3([0 XYZVZ(1,i)/VN(i)]',[0 XYZVZ(2,i)/VN(i)]',[0 XYZVZ(3,i)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','y')
% end
% 
% VN = vecnorm(XYZRVZ',2,2);
% for i = 1:length(RZ.Results.usedpullIndsR_NystagCorr)
%     lPlotLD = plot3([0 XYZRVZ(1,i)/VN(i)]',[0 XYZRVZ(2,i)/VN(i)]',[0 XYZRVZ(3,i)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','b')
% end
% 
VN = vecnorm(RZLRZV',2,2);
for i = 1:length(RZ.Results.usedpullIndsR_NystagCorr)
    lPlotLD = plot3([0 RZLRZV(1,i)/VN(i)]',[0 RZLRZV(2,i)/VN(i)]',[0 RZLRZV(3,i)/VN(i)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','g')
end
hold(avgMisalign3D,'off');




%%
% 1. Pull Raw LARP, RALP, Z Velocity data and filter using original
% parameters
% 2. Pull values at positions used previously
% 3. Calculate misalignment values for used cycles
% 4. Apply 45 degree Z rotation to cycle velocity values LRZ->XYZ
% 5. Apply desired rotation about Y axis
% 6. Apply -45 degree rotation about Z XYZ-> LRZ
% 7. Recalculate misalignment values

toRot = 8.8;

% 1.
% Larp = filtfilt(ones(1,3)/3,1,RL.Results.segmentData.RE_Velocity_LARP);
% Ralp = filtfilt(ones(1,3)/3,1,RL.Results.segmentData.RE_Velocity_RALP);
% Z = filtfilt(ones(1,3)/3,1,RL.Results.segmentData.RE_Velocity_Z);

% 2.
rpr = repmat([1 0 0],length(LL.Results.usedpullIndsL_NystagCorr),1);
% LRZ = [Larp(LL.Results.usedpullIndsR) Ralp(LL.Results.usedpullIndsR) Z(LL.Results.usedpullIndsR)];
LRZL = LL.Results.usedMisalign3DL_NystagCorr;
% 3.
misVals = acosd(dot(rpr,LRZL,2)./(vecnorm(rpr,2,2).*vecnorm(LRZL,2,2)));
[misVals LL.Results.usedMisalignL_NystagCorr'];

% 4.
XYZVL = [cosd(-45) -sind(-45) 0; sind(-45) cosd(-45) 0; 0 0 1]*LRZL';

% 5.

XYZRVL = [cosd(toRot) 0 sind(toRot); 0 1 0; -sind(toRot) 0 cosd(toRot)]*XYZVL;

% 6.
LLLRZV = [cosd(45) -sind(45) 0; sind(45) cosd(45) 0; 0 0 1]*XYZRVL;

% 7.
misVals2 = acosd(dot(rpr,LLLRZV',2)./(vecnorm(rpr,2,2).*vecnorm(LLLRZV',2,2)));
% misVals and misVals2 are equal if toRot = 0
[misVals2 misVals LL.Results.usedMisalignL_NystagCorr'];
Lold = mean(misVals)
Lnew = mean(misVals2)

% 1. Pull Raw LARP, RALP, Z Velocity data and filter using original
% parameters
% 2. Pull values at positions used previously
% 3. Calculate misalignment values for used cycles
% 4. Apply 45 degree Z rotation to cycle velocity values LRZ->XYZ
% 5. Apply desired rotation about Y axis
% 6. Apply -45 degree rotation about Z XYZ-> LRZ
% 7. Recalculate misalignment values

% 1.
% Larp = filtfilt(ones(1,3)/3,1,RR.Results.segmentData.RE_Velocity_LARP);
% Ralp = filtfilt(ones(1,3)/3,1,RR.Results.segmentData.RE_Velocity_RALP);
% Z = filtfilt(ones(1,3)/3,1,RR.Results.segmentData.RE_Velocity_Z);

% 2.
rpr = repmat([0 1 0],length(LR.Results.usedpullIndsL_NystagCorr),1);
% LRZ = [Larp(RR.Results.usedpullIndsR) Ralp(RR.Results.usedpullIndsR) Z(RR.Results.usedpullIndsR)];
LRZR = LR.Results.usedMisalign3DL_NystagCorr;
% 3.
misVals = acosd(dot(rpr,LRZR,2)./(vecnorm(rpr,2,2).*vecnorm(LRZR,2,2)));
[misVals LR.Results.usedMisalignL_NystagCorr'];

% 4.
XYZVR = [cosd(-45) -sind(-45) 0; sind(-45) cosd(-45) 0; 0 0 1]*LRZR';
% 5.
XYZRVR = [cosd(toRot) 0 sind(toRot); 0 1 0; -sind(toRot) 0 cosd(toRot)]*XYZVR;
% XYZRVR(1,:) = XYZRVR(1,:)/.7;
% 6.
LRLRZV = [cosd(45) -sind(45) 0; sind(45) cosd(45) 0; 0 0 1]*XYZRVR;

% 7.
misVals2 = acosd(dot(rpr,LRLRZV',2)./(vecnorm(rpr,2,2).*vecnorm(LRLRZV',2,2)));
[misVals2 LR.Results.usedMisalignL_NystagCorr'];
Rold = mean(misVals)
Rnew = mean(misVals2)

% misVals and misVals2 are equal if toRot = 0

% 1. Pull Raw LARP, RALP, Z Velocity data and filter using original
% parameters
% 2. Pull values at positions used previously
% 3. Calculate misalignment values for used cycles
% 4. Apply 45 degree Z rotation to cycle velocity values LRZ->XYZ
% 5. Apply desired rotation about Y axis
% 6. Apply -45 degree rotation about Z XYZ-> LRZ
% 7. Recalculate misalignment values

% 1.
% Larp = filtfilt(ones(1,3)/3,1,RR.Results.segmentData.RE_Velocity_LARP);
% Ralp = filtfilt(ones(1,3)/3,1,RR.Results.segmentData.RE_Velocity_RALP);
% Z = filtfilt(ones(1,3)/3,1,RR.Results.segmentData.RE_Velocity_Z);

% 2.
rpr = repmat([0 0 -1],length(LZ.Results.usedpullIndsL_NystagCorr),1);
% LRZ = [Larp(RR.Results.usedpullIndsR) Ralp(RR.Results.usedpullIndsR) Z(RR.Results.usedpullIndsR)];
LRZZ = LZ.Results.usedMisalign3DL_NystagCorr;
% 3.
misVals = acosd(dot(rpr,LRZZ,2)./(vecnorm(rpr,2,2).*vecnorm(LRZZ,2,2)));
[misVals LZ.Results.usedMisalignL_NystagCorr'];

% 4.
XYZVZ = [cosd(-45) -sind(-45) 0; sind(-45) cosd(-45) 0; 0 0 1]*LRZZ';
% 5.
XYZRVZ = [cosd(toRot) 0 sind(toRot); 0 1 0; -sind(toRot) 0 cosd(toRot)]*XYZVZ;

% 6.
LZLRZV = [cosd(45) -sind(45) 0; sind(45) cosd(45) 0; 0 0 1]*XYZRVZ;

% 7.
misVals2 = acosd(dot(rpr,LZLRZV',2)./(vecnorm(rpr,2,2).*vecnorm(LZLRZV',2,2)));
[misVals2 LZ.Results.usedMisalignL_NystagCorr'];
Zold = mean(misVals)
Znew = mean(misVals2)

% misVals and misVals2 are equal if toRot = 0
%%
avgMisalignPlot3D = figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle(avgMisalignPlot3D,{'3D Angle of Misalignment'},'FontSize', 22, 'FontWeight', 'Bold');

avgMisalign3D = axes('Parent', avgMisalignPlot3D);
avgMisalign3D.XGrid = 'on';
avgMisalign3D.YGrid = 'on';
avgMisalign3D.FontSize = 13.5;

hold(avgMisalign3D,'on');
h=plot3vect([1;0;0],'LARP Axis',[0 1 0],2);
set(h,'LineStyle','--','Marker','o');
h=plot3vect([0;1;0],'RALP Axis',[0 0 1],2);
set(h,'LineStyle','--','Marker','o');
h=plot3vect([0;0;-1],'Yaw Axis',[1 0 0],2);
set(h,'LineStyle','--','Marker','o');



[x,y,z]=sphere();
h=surf(0.5*x,0.5*y,0.5*z);
set(h,'FaceColor','white')
avgMisalign3D.View = [135 -15];
axis vis3d
axis equal
box on;
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])

VN = vecnorm(LRZR,2,2);
for i = 1:length(LR.Results.usedpullIndsL_NystagCorr)
    lPlotLD = plot3([0 LRZR(i,1)/VN(i)]',[0 LRZR(i,2)/VN(i)]',[0 LRZR(i,3)/VN(i)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','k')
end

% VN = vecnorm(XYZVR',2,2);
% for i = 1:length(LR.Results.usedpullIndsL_NystagCorr)
%     lPlotLD = plot3([0 XYZVR(1,i)/VN(i)]',[0 XYZVR(2,i)/VN(i)]',[0 XYZVR(3,i)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','y')
% end
% 
% VN = vecnorm(XYZRVR',2,2);
% for i = 1:length(LR.Results.usedpullIndsL_NystagCorr)
%     lPlotLD = plot3([0 XYZRVR(1,i)/VN(i)]',[0 XYZRVR(2,i)/VN(i)]',[0 XYZRVR(3,i)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','b')
% end

VN = vecnorm(LRLRZV',2,2);
for i = 1:length(LR.Results.usedpullIndsL_NystagCorr)
    lPlotLD = plot3([0 LRLRZV(1,i)/VN(i)]',[0 LRLRZV(2,i)/VN(i)]',[0 LRLRZV(3,i)/VN(i)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','g')
end

VN = vecnorm(LRZL,2,2);
for i = 1:length(LL.Results.usedpullIndsL_NystagCorr)
    lPlotLD = plot3([0 LRZL(i,1)/VN(i)]',[0 LRZL(i,2)/VN(i)]',[0 LRZL(i,3)/VN(i)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','k')
end

% VN = vecnorm(XYZVL',2,2);
% for i = 1:length(LL.Results.usedpullIndsL_NystagCorr)
%     lPlotLD = plot3([0 XYZVL(1,i)/VN(i)]',[0 XYZVL(2,i)/VN(i)]',[0 XYZVL(3,i)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','y')
% end
% 
% VN = vecnorm(XYZRVL',2,2);
% for i = 1:length(LL.Results.usedpullIndsL_NystagCorr)
%     lPlotLD = plot3([0 XYZRVL(1,i)/VN(i)]',[0 XYZRVL(2,i)/VN(i)]',[0 XYZRVL(3,i)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','b')
% end

VN = vecnorm(LLLRZV',2,2);
for i = 1:length(LL.Results.usedpullIndsL_NystagCorr)
    lPlotLD = plot3([0 LLLRZV(1,i)/VN(i)]',[0 LLLRZV(2,i)/VN(i)]',[0 LLLRZV(3,i)/VN(i)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','g')
end


VN = vecnorm(LRZZ,2,2);
for i = 1:length(LZ.Results.usedpullIndsL_NystagCorr)
    lPlotLD = plot3([0 LRZZ(i,1)/VN(i)]',[0 LRZZ(i,2)/VN(i)]',[0 LRZZ(i,3)/VN(i)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','k')
end

% VN = vecnorm(XYZVZ',2,2);
% for i = 1:length(LZ.Results.usedpullIndsL_NystagCorr)
%     lPlotLD = plot3([0 XYZVZ(1,i)/VN(i)]',[0 XYZVZ(2,i)/VN(i)]',[0 XYZVZ(3,i)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','y')
% end
% 
% VN = vecnorm(XYZRVZ',2,2);
% for i = 1:length(LZ.Results.usedpullIndsL_NystagCorr)
%     lPlotLD = plot3([0 XYZRVZ(1,i)/VN(i)]',[0 XYZRVZ(2,i)/VN(i)]',[0 XYZRVZ(3,i)/VN(i)]');
%     set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','b')
% end

VN = vecnorm(LZLRZV',2,2);
for i = 1:length(LZ.Results.usedpullIndsL_NystagCorr)
    lPlotLD = plot3([0 LZLRZV(1,i)/VN(i)]',[0 LZLRZV(2,i)/VN(i)]',[0 LZLRZV(3,i)/VN(i)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','g')
end
hold(avgMisalign3D,'off');
%%
avgMisalignPlot3D = figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle(avgMisalignPlot3D,{'3D Angle of Misalignment'},'FontSize', 22, 'FontWeight', 'Bold');

avgMisalign3D = axes('Parent', avgMisalignPlot3D);
avgMisalign3D.XGrid = 'on';
avgMisalign3D.YGrid = 'on';
avgMisalign3D.FontSize = 13.5;

hold(avgMisalign3D,'on');
h=plot3vect([1;0;0],'LARP Axis',[0 1 0],2);
set(h,'LineStyle','--','Marker','o');
h=plot3vect([0;1;0],'RALP Axis',[0 0 1],2);
set(h,'LineStyle','--','Marker','o');
h=plot3vect([0;0;-1],'Yaw Axis',[1 0 0],2);
set(h,'LineStyle','--','Marker','o');



[x,y,z]=sphere();
h=surf(0.5*x,0.5*y,0.5*z);
set(h,'FaceColor','white')
avgMisalign3D.View = [135 -15];
axis vis3d
axis equal
box on;
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
%%
    lPlotLD = plot3([0 tempS(1).plotM3DL_NystagCorr(1)]',[0 tempS(1).plotM3DL_NystagCorr(2)]',[0 tempS(1).plotM3DL_NystagCorr(3)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','k')


VN = vecnorm(tempS(1).M3DL_NystagCorr,2,2);
for i = 1:length(VN)
    lPlotLD = plot3([0 tempS(1).M3DL_NystagCorr(i,1)/VN(i)]',[0 tempS(1).M3DL_NystagCorr(i,2)/VN(i)]',[0 tempS(1).M3DL_NystagCorr(i,3)/VN(i)]');
    set(lPlotLD,'LineWidth',3.5,'DisplayName','','Color','k')
end
%%
x = Results.segmentData.RE_Position_X;%(Results.usedpullIndsR);
y = Results.segmentData.RE_Position_Y;%(Results.usedpullIndsR);
z = Results.segmentData.RE_Position_Z;%(Results.usedpullIndsR);
XYZ = [x';y';z'];
toRot = 16;
XYZRot = [cosd(toRot) 0 sind(toRot); 0 1 0; -sind(toRot) 0 cosd(toRot)]*XYZ;

Data_In = struct();
Data_In.Data_LE_Pos_X = XYZRot(1,:)';
Data_In.Data_LE_Pos_Y = XYZRot(1,:)';
Data_In.Data_LE_Pos_Z = XYZRot(1,:)';

Data_In.Data_RE_Pos_X = XYZRot(1,:)';
Data_In.Data_RE_Pos_Y = XYZRot(2,:)';
Data_In.Data_RE_Pos_Z = XYZRot(3,:)';

Data_In.Fs = 1000;
data_rot = 1;
OutputFormat = [];
DAQ_code = 7;
[New_Ang_Vel] = voma__processeyemovements([],[],[],[],0,data_rot,DAQ_code,OutputFormat,Data_In);
L = New_Ang_Vel.RE_Vel_LARP(Results.usedpullIndsR);
R = New_Ang_Vel.RE_Vel_RALP(Results.usedpullIndsR);
Z2 = New_Ang_Vel.RE_Vel_Z(Results.usedpullIndsR);

rpr = repmat([1 0 0],length(Z2),1);
LRZTransformed2 = [L R Z2];
misVals2 = acosd(dot(rpr,LRZTransformed2,2)./(vecnorm(rpr,2,2).*vecnorm(LRZTransformed2,2,2)));
%%
Xv = Results.segmentData.RE_Velocity_X(Results.usedpullIndsR);
Yv = Results.segmentData.RE_Velocity_Y(Results.usedpullIndsR);
Zv = Results.segmentData.RE_Velocity_Z(Results.usedpullIndsR);

XYZV = [Xv';Yv';Zv'];
toRot = 16;
XYZRV = [cosd(toRot) 0 sind(toRot); 0 1 0; -sind(toRot) 0 cosd(toRot)]*XYZV;
d = [cosd(-45) -sind(-45) 0; sind(-45) cosd(-45) 0; 0 0 1]*XYZRV
%%
rpr = repmat([1 0 0],length(Z),1);
LRZTransformed = [Larp Ralp Z];
misVals = acosd(dot(rpr,LRZTransformed,2)./(vecnorm(rpr,2,2).*vecnorm(LRZTransformed,2,2)));

%%
z =filtfilt(ones(1,3)/3,1,LR.Results.segmentData.LE_Velocity_Z);
l =filtfilt(ones(1,3)/3,1,LR.Results.segmentData.LE_Velocity_LARP);
r =filtfilt(ones(1,3)/3,1,LR.Results.segmentData.LE_Velocity_RALP);
figure
plot(l,'color','g')
hold on
plot(LR.Results.usedpullIndsL,l(LR.Results.usedpullIndsL),'k*')

plot(r,'color','b')
plot(LR.Results.usedpullIndsL,r(LR.Results.usedpullIndsL),'k*')
plot(z,'color','r')
plot(LR.Results.usedpullIndsL,z(LR.Results.usedpullIndsL),'k*')