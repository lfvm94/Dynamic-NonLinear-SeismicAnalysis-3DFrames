function [Coordxyz,PointsElm]=ExtractTopologySAP2000(ProgramPath, ...
                              APIDLLPath, ModelPath)
% [Coordxyz,PointsElm]=ExtractTopologySAP2000(ProgramPath,APIDLLPath,...
%                                             ModelPath)
%---------------------------------------------------------------------
%    PURPOSE
%     To extract the node connectivity and node coordinates from a 
%     SAP2000 model with the aid of the SM toolbox.
% 
%    INPUT:  ProgramPath:       location of the SAP2000 software in disc
%
%            APIDLLPath:        location of the SAP DLL file in disc
%
%            ModelPath:         path of the existing SAP2000 file from
%                               which to extract the information
%
%    OUTPUT: Coordxyz :         node coordinates [x,y,z]
%
%            PointsElm :        node connectivity of each element
%
%--------------------------------------------------------------------
%    Notes: 
%           To download the SM Toolbox visit: 
%           https://github.com/RJ-Soft/SM-Toolbox/releases/tag/7.0.2
%--------------------------------------------------------------------

% LAST MODIFIED: L.F.Verduzco    2023-05-01
% Copyright (c)  Faculty of Engineering.
%                Autonomous University of Queretaro
%--------------------------------------------------------------------

%% Determine the type of application and its version

SM.App( 'sap' );

SM.Ver( '22' );

%% Create Sap2000 object

[ Sobj ] = SM.Helper.CreateObject( ProgramPath,APIDLLPath );

%% Create SapModel object

[ Smdl ]=SM.SapModel();

%% Start Sap2000 application

SM.ApplicationStart;

%% Initialize model

ret1 = SM.InitializeNewModel;

% Option 2: open an existing file
ret = SM.File.OpenFile(ModelPath);

%% Create the analysis model
ret4 = SM.Analyze.CreateAnalysisModel;

%% Get point element names

nnodes = SM.PointElm.Count;
% Get cartesian point element coordinates
for i=1:nnodes
    node=num2str(i);
    [ret,x,y,z]= SM.PointElm.GetCoordCartesian(node);
    Coordxyz(i,:)=[x,y,z];
end
%% Get frame object names

[ret8,NumberElm,MyNameElm]=SM.FrameObj.GetNameList();

for i=1:NumberElm
    elm=num2str(i);
    
    % get names of points of each element
    [ret3,p1,p2]=SM.FrameObj.GetPoints(elm);
    PointsElm(i,1)=str2num(p1);
    PointsElm(i,2)=str2num(p2);
    
end

%% Close Sap2000
SM.ApplicationExit(false);

% ------------------------------ End ----------------------------------