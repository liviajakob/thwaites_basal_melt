%% Read data

moor = csvread('./e_Moorings_0_8C_depth.csv',1,0) ;
moor(1108,:) = [] ; moor(1930:1945,:) = [] ; % Remove duplicate and no values
sglV = csvread('./2d_Subglacial_lakes_volume_change.csv',1,0) ;
melt = csvread('./e_Thwaites_basal_melt_anomalies.csv',1,0) ;
meltPIG = csvread('./e_Pine_Island_basal_melt_anomalies.csv',1,0) ;

%% Resample to melt time sampling

moorR = interp1(moor(:,1),moor(:,2),melt(:,1)) ;
sglVR = interp1(sglV(:,1),sglV(:,2),melt(:,1)) ;

%% Compute flux

sglFluxR = -gradient(sglVR,melt(:,1))*1e9/365.25/24/3600 ;

%% Build data matrix
 
Data=[moorR sglFluxR melt(:,2) meltPIG(:,2)];

%% Perform correlation
figure;[R,Pv,H]=corrplot(Data) ; %,'TestR','on') ;

H(2,1).Parent.XLim=H(1,1).Parent.XLim; H(3,1).Parent.XLim=H(1,1).Parent.XLim; H(4,1).Parent.XLim=H(1,1).Parent.XLim;
H(3,2).Parent.XLim=H(2,2).Parent.XLim; H(4,2).Parent.XLim=H(2,2).Parent.XLim;
H(4,3).Parent.XLim=H(3,3).Parent.XLim;

subplot(4,4,1) ;  delete(gca) ;  subplot(4,4,2) ;  delete(gca) ;  subplot(4,4,3);delete(gca) ;  subplot(4,4,4);delete(gca) ; 
subplot(4,4,6) ;  delete(gca) ; subplot(4,4,7) ;  delete(gca) ;  subplot(4,4,8);delete(gca) ;  
subplot(4,4,11) ;  delete(gca) ;  subplot(4,4,12) ;  delete(gca) ;
subplot(4,4,16) ;  delete(gca) ;

H(2,1).Parent.YLabel.String = 'SGL flux (m^3 s^-^1)' ;  H(3,1).Parent.YLabel.String = 'T. m.a. (m yr^-^1)' ;  H(4,1).Parent.YLabel.String = 'PIG m.a. (m yr^-^1)' ;
H(4,1).Parent.XLabel.String = 'Therm. d. (m)' ;  H(4,2).Parent.XLabel.String = 'SGL flux (m^3 s^-^1)' ;  H(4,3).Parent.XLabel.String = 'T. m.a. (m yr^-^1)' ;

%% Lag + correlation 
[xcf,lags] = crosscorr(Data(:,2),Data(:,3)) ;
sglFluxRLag = circshift(sglFluxR,lags(find(xcf==max(xcf)))) ;
DataLag = Data ; DataLag(:,2) = sglFluxRLag ;
fhl = figure;[RLag,PvLag,HLag]=corrplot([sglFluxRLag melt(:,2)]); 

subplot(2,2,1);delete(gca) ; subplot(2,2,2);delete(gca) ; subplot(2,2,4) ;delete(gca) ;
HLag(2,1).Parent.XLabel.String = 'SGL flux (m^3 s^-^1)' ;
HLag(2,1).Parent.YLabel.String = 'T. m.a. (m yr^1)' ;

fluxvect = 0:40:590 ;  % Resample to get equal distribution along the flux vector and avoid overfitting by data cluster around low subglacial flux

M = interp1(sglFluxRLag,melt(:,2),fluxvect) ;
f=fit(fluxvect',M','b*x^m');

figure(fhl);hold on;plot(f,'k--');legend('off')
HLag(2,1).Parent.XLabel.String = 'SGL flux (m^3 s^-^1)' ;
HLag(2,1).Parent.YLabel.String = 'T. m.a. (m yr^1)' ;

%% Correlation between PIG and Thwaites melt after the transient melt increase i.e. post year 2015

fhl = figure;[Rm,Pm,Hm]=corrplot([Data(35:end,3) Data(35:end,4)]) ;
subplot(2,2,1);delete(gca) ; subplot(2,2,2);delete(gca) ; subplot(2,2,4) ;delete(gca) ;
Hm(2,1).Parent.XLabel.String = 'T. m.a. (m yr^1)' ;
Hm(2,1).Parent.YLabel.String = 'PIG m.a. (m yr^1)' ;

