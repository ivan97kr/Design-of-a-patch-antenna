%%Create pcbStack object
pcbobj = pcbStack;

%%Create board shape
    %Creating BoardShape metal layer.
        %Creating Rectangle1 shape.
        Rectangle1 = antenna.Rectangle;
        Rectangle1.Name = "Rectangle1";
        Rectangle1.Center = [0 0];
        Rectangle1.Length = 0.0218;
        Rectangle1.Width = 0.0467;
        Rectangle1 = rotate(Rectangle1,0,[0 0 -1],[0 0 1]);
    BoardShape = Rectangle1;
pcbobj.BoardShape = BoardShape;

%%Create Stackup
    %Creating MetalLayer2 metal layer.
        %Creating Rectangle3 shape.
        Rectangle3 = antenna.Rectangle;
        Rectangle3.Name = "Rectangle3";
        Rectangle3.Center = [0 0];
        Rectangle3.Length = 0.0121;
        Rectangle3.Width = 0.0371;
        Rectangle3 = rotate(Rectangle3,0,[0 0 -1],[0 0 1]);
    MetalLayer2 = Rectangle3;
    %Creating DielectricLayer1 dielectric layer.
    DielectricLayer1 = dielectric("Name",'Custom',"EpsilonR",4.4,"LossTangent",0.02,"Thickness",0.0016);
    %Creating MetalLayer1 metal layer.
        %Creating Rectangle4 shape.
        Rectangle4 = antenna.Rectangle;
        Rectangle4.Name = "Rectangle4";
        Rectangle4.Center = [0 0];
        Rectangle4.Length = 0.0218;
        Rectangle4.Width = 0.0467;
        Rectangle4 = rotate(Rectangle4,0,[0 0 -1],[0 0 1]);
    MetalLayer1 = Rectangle4;

%%Create Feed
feedloc = [[-0.0025 0 3 1];...
    ];

%%Create Via
vialoc = [[-0.005 0.0129 3 1];...
[-0.005 -0.0129 3 1];...
[-0.005 0.0157 3 1];...
[-0.005 -0.0165 3 1];...
[-0.005 0 3 1];...
[-0.005 0.0057 3 1];...
[-0.005 -0.0057 3 1];...
[-0.005 0.0093 3 1];...
[-0.005 -0.0093 3 1]
    ];

%%Create Metal
metalobj = metal;
metalobj.Name = 'PEC';
metalobj.Conductivity = Inf;
metalobj.Thickness = 0; % 0 mils

pcbobj.Conductor = metalobj;

%%Assign properties
pcbobj.BoardThickness = 0.0016;
pcbobj.Layers = {MetalLayer2,DielectricLayer1,MetalLayer1,};
pcbobj.FeedLocations = feedloc;
pcbobj.FeedDiameter = 0.0005;
pcbobj.ViaLocations = vialoc;
pcbobj.ViaDiameter = 0.0006;
pcbobj.FeedViaModel = 'strip';
pcbobj.FeedVoltage = 1;
pcbobj.FeedPhase = 0;
