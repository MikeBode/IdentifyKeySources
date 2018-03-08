clear all

PT = 0; % 0 = no plotting; 1 = plotting.

% =====================================================================================
% ======== Start by matching the reef identities in the two biophysical models ========
% =====================================================================================

load HockData             % Load the spatial data from the supplemental material in Hock et al. 2017 (H)
% This MAT file contains (1) the coordinates of each reef in the CONNIE model,
% and a list of the key source reefs identified in that paper.

load original_centroids   % Load the equivalent spatial data from Bode et al. 2012 (A)
centroid = [lg lt];

% Go through the reefs in H one-by-one. A reef Hi is matched with reef Aj if
% (i)  d(Hi,Aj) <= d(Hi,Ak) for all k
% (ii) d(Hi,Aj) <  TD, where TD is a threshold distance
TD = 4/105; % Convert 4 km to degrees at latitude 19

A_Match = nan.*ones(length(centroid),1);
H_Match = nan.*ones(length(ReefLat),1);
for i = 1:length(centroid); A_Match_cell{i,1} = []; end

for i = 1:length(ReefLon)
   if mod(i,500) == 1; disp(['Completed: ' num2str(i)]); end
   
   % Compute the Euclidean distance between all reefs in A and reef Hi
   D = pdist2([ReefLon(i) ReefLat(i)],centroid);
   
   % Assign Hi to Aj if they're closest, and close enough, and record the match
   if min(D) < TD
      F = find(D == min(D));
      H_Match(i) = F(1);
      A_Match(F(1)) = i;
      AMC = [A_Match_cell{F(1)} i]; 
      A_Match_cell{F(1)} = AMC;
   end
end

% Generate summary statistics about the number of matches between the reefs in the two models
Missing_H_Reefs = sum(isnan(H_Match))
Missing_H_Reefs_prop = sum(isnan(H_Match))./length(H_Match)
Missing_A_Reefs = sum(isnan(A_Match))
Missing_A_Reefs_prop = sum(isnan(A_Match))./length(A_Match)

% =====================================================================================
% ========= Identify key sources in model A using the methods in Hock et al.  =========
% =====================================================================================

% Load connectivity matrices for model A
for y = 1996:2002; 
   eval(['load ConnectivityMatrices_Model_A_' num2str(y) '_P5;']); 
   eval(['load ConnectivityMatrices_Model_A_' num2str(y) '_P7;']); 
end
W = whos('*psurv*');
NumReefs_A = size(eval(W(1).name),1); % Number of reefs in model A
Edges = 1:NumReefs_A;

% Discard any reefs in model A that are outside the domain of the GBR Marine Park (to match model H)
DiscardReefs = find(centroid(:,2) > -10.3355);

% Go through the matrices one-by-one
for w = 1:length(W)
   
   % Choose a matrix (a "scenario") from the set
   W(w).name
   eval(['ThisMat = ' W(w).name ';']);
   
   % CRITERION 1: Calculate the best out-degree reefs
   OutDegree = sum(ThisMat > 0, 2);
   OutDegree(DiscardReefs) = -inf;
   V = OutDegree; V(DiscardReefs) = [];
   Best_criterion1 = find(OutDegree > median(V));
   Worst_criterion1 = find(OutDegree <= median(V));
   
   % CRITERION 2: Calculate the best total output reefs
   TotalOutput = sum(ThisMat,2);
   TotalOutput(DiscardReefs) = -inf;
   V = TotalOutput; V(DiscardReefs) = [];
   Best_criterion2 = find(TotalOutput > median(V));
   Worst_criterion2 = find(TotalOutput <= median(V));
   
   % CRITERION 3: Calculate the reefs that provide the strongest (>10%) links to other reefs
   RelativeSupply = ThisMat./repmat(sum(ThisMat),NumReefs_A,1); % How much does each reef contribute to the other reefs?
   RelativeSupply = RelativeSupply > 0.1; % Are any of those connections larger than the stated threshold (10%)?
   StrongLinks = sum(RelativeSupply,2); % Add up all the strong out connections from each reef
   StrongLinks(DiscardReefs) = -inf;
   V = StrongLinks; V(DiscardReefs) = [];
   Best_criterion3 = find(StrongLinks > median(V));
   Worst_criterion3 = find(StrongLinks <= median(V));
   
   % CRITERION 4: Calculate source reefs that are linked to other source reefs (identified by criteria 1,2,3)
   % NB: We don't need to discard reefs outside the GBRMP boundary from this criterion, because they wouldn't be intermediate sources
   SourcesSoFar = intersect(Best_criterion3,...
      intersect(Best_criterion1,Best_criterion2));
   NotSourcesSoFar = setdiff(1:NumReefs_A,SourcesSoFar); % These are not sources according to criteria 1,2,3 (call them intermediate sources)
   TM = ThisMat; % Duplicate the connectivity matrix
   TM(NotSourcesSoFar,:) = 0; TM(:,NotSourcesSoFar) = 0; % Remove all connections that don't come from, or go to, an intermediate source
   SourceSource_connections = sum(TM > 0,2); % Sum how many connections out of each intermediate source go to another intermediate source
   Best_criterion4 = find(SourceSource_connections > median(SourceSource_connections));
   Worst_criterion4 = find(SourceSource_connections <= median(SourceSource_connections));
   
   % CRITERION 5: Calculate the number of reefs in the network that can be reached via a directed path from each reef
   % NB: We have not calculated this criteria, because all connectivity matrices are effectively strongly connected. That is,
   %     it is possible to reach any reef in the system by a directed route from any other reef in the system. This means that
   %     all the reefs would be in the top (or bottom) 50% percentile.
   
   % Amalgamate the criteria for this connectivity matrix. Search for reefs that are in the top 50% of every criteria
   % Also search for reefs that are the worst sources (in the bottom 50% of every criteria).
   AA = intersect(Best_criterion4,...
      intersect(Best_criterion3,...
      intersect(Best_criterion1, Best_criterion2)));
   KeySources_EachMatrix(:,w) = histc(AA,Edges); % Store the best sources for this connectivity matrix
   
   BB = intersect(Worst_criterion4,...
      intersect(Worst_criterion3,...
      intersect(Worst_criterion1, Worst_criterion2)));
   BadSources_EachMatrix(:,w) = histc(BB,Edges); % Store the worst sources for this connectivity matrix
end

% Identify key sources in model A by looking across the different scenarios (matrices).
% Key sources are reefs that fit all 5 criteria in an above median number of scenarios.
% We identify the "worst sources" using the inverse definition.
KS_model_A = sum(KeySources_EachMatrix,2);
A_Sources = find(KS_model_A > median(KS_model_A));

BS_model_A = sum(BadSources_EachMatrix,2);
A_Worst = find(BS_model_A > median(BS_model_A));

% ======================================================================================
% ========= Check if Key Sources in model A are also Key Sources in model H  =========
% ======================================================================================

% Key sources in A that are also in H
KeySourcesA_DefH = [];
for a = 1:length(A_Sources)  % Note that a single key source in model A may match multiple reefs in model H
   KeySourcesA_DefH = [KeySourcesA_DefH; A_Match_cell{A_Sources(a)}'];
end
SourcesInBoth = unique(intersect(KeySourcesA_DefH,KeySources));
Num_SourcesInBoth = length(SourcesInBoth);
Prop_SourcesInBoth_AinH = Num_SourcesInBoth./length(A_Sources)

% Which of the Key Sources in H are the worst sources in A?
KeySourcesH_DefA = H_Match(KeySources);
WorstA_KeyH = unique(intersect(A_Worst,KeySourcesH_DefA));
Num_WorstA_KeyH = length(WorstA_KeyH);
Prop_WorstA_KeyH = Num_WorstA_KeyH./length(KeySources)

% ======================================================
% ================== Plot the results ==================
% ======================================================
load original_centroids; centroid = [lg lt];
load HockData Reef*

Th = 0.66; Rot = [cos(Th) -sin(Th); sin(Th) cos(Th)]; % Create a rotation matrix
figure(1), clf, subplot('position',[0.06 0.06 0.89 0.89]); hold on, box on; MS = 10; FS = 11;
CL = get(gca,'colororder'); CL1 = CL(1,:); CL2 = CL(3,:); CL3 = CL(7,:);
A_base = [145.5 -12.5];

LandCol = [70 140 67]./243;
PlotAustralianOutline
centroid = centroid(:,1:2)*Rot;
A = [ReefLon,ReefLat]*Rot; ReefLon = A(:,1); ReefLat = A(:,2);
plot(ReefLon,ReefLat,'.','markersize',MS-5,'color',0.5.*ones(1,3))
plot(ReefLon(KeySources),ReefLat(KeySources),'.','markersize',MS,'color',CL1)
plot(centroid(A_Sources,1),centroid(A_Sources,2),'.','markersize',MS,'color',CL2)
plot(ReefLon(SourcesInBoth),ReefLat(SourcesInBoth),'.','markersize',MS,'color',CL3)
ylim([min(ReefLat) max(ReefLat)])
ylim([-113 -95])
xlim([103.5 108.5])

X = 106;
Y = -102.5; plot(X,Y,'.','markersize',MS*2,'color',CL1); text(X+0.5,Y,'Model H','fontsize',FS)
Y = -103.5; plot(X,Y,'.','markersize',MS*2,'color',CL2); text(X+0.5,Y,'Model A','fontsize',FS)
Y = -104.5; plot(X,Y,'.','markersize',MS*2,'color',CL3) ; text(X+0.5,Y,'Both','fontsize',FS)
set(gca,'xtick',[],'ytick',[])

set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 11 35]*0.6)
set(gcf, 'renderer', 'painters')
print('-dtiff','-r400','KeySourcePerformance.tiff')
















