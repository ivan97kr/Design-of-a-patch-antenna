%Starting point, Balanis formulations adaptated to folded patch
f = 2.7e9;  
c = physconst('Lightspeed');
lambda_0 = c/f;
k_0 = (2*pi)/lambda_0;
e_r = 4.4;
W = (lambda_0/2) * sqrt( (2 / (e_r+1) ) );
h_max = (0.3 * lambda_0) / ( 2*pi*sqrt(e_r) );
h = 0.0016;
e_e = ( (e_r+1) / 2 ) + ( ( (e_r - 1) / 2 ) * ( 1 + (12*(h/W)))^-0.5 );
Delta_L = 0.412 * h *  ( (( (e_e+0.3) * ( (W/h) + 0.264) )) / ( (e_e - 0.258) * ( (W/h) + 0.8)) );
R_r = (120*lambda_0)/(W*(1-((1/24)*(k_0*h)^2)));
R_in = 50;
L_eff = (lambda_0/(4*sqrt(e_e)));
L = L_eff - (Delta_L);
beta = pi/(2*L);
L_1 = (1/beta)*acos(sqrt(R_in/R_r));

%Lunghezza del ground del pifa
Lg = L + (6*h);
%Larghezza del gorund del pifa
Wg = W + (6*h);
%punto di alimentazion nel sistema di riferimento di matlab
FeedOffSet = L/2 - L_1;
freq = linspace(2e9,3.3e9,14);
% !!! pifa gia caricato manualemnte !!! tramite antennaDesigner e le dim.
% ottenute dalle formule


% s1 = sparameters(pf,freq);
% rfplot(s1)



D_0 = (1/30)*((W/lambda_0)^2)*R_r;
% BWe_rad = 2 * acos(sqrt((7.03*lambda_0^2)/(4*((3*(L_eff^2))+ h^2)*(pi^2))))
% BWh_rad = 2 * acos(sqrt (  1  / (2*(k_0*W) ) ))
% 
% BWe = (BWe_rad *180)/pi
% BWh = (BWh_rad *180)/pi



val=zeros(1,14);
frequencies= zeros(1,14);

%scelgo la mesh
for i = 1:90
meshData = mesh (pifa_to_mesh,'MaxEdgeLength',lambda_0/(i));
%returnLoss(antennaObject,[2430e6:27e6:2970e6])
s=sparameters(pifa_to_mesh,freq);

[val(i),indice]=min(20*log10(abs(s.Parameters)));
frequencies(i) = freq(indice);
%rfplot(s);
eta_R(i)=20*log10(abs(s.Parameters(8)));

count = i
end
plot(eta_R)
plot(frequencies)
grid on;

%scelgo mesh 30 samples per lambda
% meshData = mesh (pf,'MaxEdgeLength',lambda_0/30);
% s1 = sparameters(pf,freq);


 

% Optimization
% 
% %feedoffset di -0.0017 è il risultato migliore, provo adesso a variare L e
% %W
% 
% 
% L tra 0.0100 e 0.0160 40 valori mentre W tra 0.0300 e 0.0376 sempre 40
% valori
%
%
Loffset = linspace(0.0100,0.0160,40);
Woffset = linspace(0.0300,0.0376,40);
S11_coarse = zeros(numel(Loffset),numel(Woffset));
for m = 1:length(Loffset)
    m
    pf.Length = Loffset(m);
    pf.GroundPlaneLength = Loffset(m)+(6*h)+0.0001;
    pf.FeedOffset = [((-Loffset(m)/2)+0.0036) 0];

    for n = 1:length(Woffset)    
        n

        pf.GroundPlaneWidth = Woffset(n)+(6*h);
        pf.Width =  Woffset(n);
        pf.ShortPinWidth = Woffset(n);
        S = sparameters(pf,2.7e9);
        S11_fine(m,n) = 20*log10(abs(S.Parameters));
    end
end

[X,Y] = meshgrid(Loffset,Woffset);
S11_coarse_transpose = transpose(S11_fine)
contour(X,Y,S11_coarse_transpose)




%refine a grana sottile
Loffset = linspace(0.0118,0.0125,40);
Woffset = linspace(0.0350,0.0376,40);
S11_fine = zeros(numel(Loffset),numel(Woffset));
for m = 1:length(Loffset)
    m
    pf.Length = Loffset(m);
    pf.GroundPlaneLength = Loffset(m)+(6*h)+0.0001;
    pf.FeedOffset = [((-Loffset(m)/2)+0.0036) 0];

    for n = 1:length(Woffset)    
        n

        pf.GroundPlaneWidth = Woffset(n)+(6*h);
        pf.Width =  Woffset(n);
        pf.ShortPinWidth = Woffset(n);
        S = sparameters(pf,2.7e9);
        S11_fine(m,n) = 20*log10(abs(S.Parameters));
    end
end

[X,Y] = meshgrid(Loffset,Woffset);
S11_fine_transpose = transpose(S11_fine)
contour(X,Y,S11_fine_transpose)




%1D refining 1 dimensione alla volta... meno efficente (prima su L poi su
%W)
% coarse refining L
Loffset = linspace((-L/2+L),(L+L/2),40);
S11_L = zeros(1,numel(Loffset));
for m =1:length(Loffset)
    pf.Length = Loffset(m);
    pf.GroundPlaneLength = Loffset(m)+(6*h);
    pf.FeedOffset = [((-Loffset(m)/2)+0.0038) 0];
    S = sparameters(pf,freq);
    S11_L(m) = 20*log10(abs(S.Parameters(8)))
end
plot(Loffset,S11_L)
[S11_min_L,index_s11_L] = min(S11_L);
L_opt = Loffset(index_s11_L)
pf.Length = L_opt;
pf.FeedOffset = [(-L_opt/2)+0.0038 0]
pf.GroundPlaneLength = L_opt+(6*h);

% thin refining 
Loffset = linspace(0.0100,0.0130,40);
S11_L = zeros(1,numel(Loffset));
for m =1:length(Loffset)
    pf.Length = Loffset(m);
    pf.GroundPlaneLength = Loffset(m)+(6*h);
    pf.FeedOffset = [((-Loffset(m)/2)+0.0038) 0];
    S = sparameters(pf,freq);
    S11_L(m) = 20*log10(abs(S.Parameters(8)))
end
plot(Loffset,S11_L)
[S11_min_L,index_s11_L] = min(S11_L);
L_opt = Loffset(index_s11_L)
pf.Length = L_opt;
pf.FeedOffset = [(-L_opt/2)+0.0038 0]
pf.GroundPlaneLength= L_opt+(6*h);


% coarse refining W
Woffset = linspace(((-W/2)+W),((W/2)+W),40);
S11_W = zeros(1,numel(Loffset));
for m =1:length(Loffset)
    pf.GroundPlaneWidth = Woffset(m)+(6*h);
    pf.Width =  Woffset(m);
    pf.ShortPinWidth = Woffset(m);
    S = sparameters(pf,freq);
    S11_W(m) = 20*log10(abs(S.Parameters(8)))
end
plot(Woffset,S11_W)
[S11_min_W,index_s11_W] = min(S11_W);
W_opt = Woffset(index_s11_W)
pf.Width = W_opt;
pf.GroundPlaneWidth = W_opt + (6*h);
pf.ShortPinWidth = W_opt;

%  fine refining
Woffset = linspace(0.0310,0.038,40);
S11_W = zeros(1,numel(Loffset));
for m =1:length(Loffset)
    pf.GroundPlaneWidth = Woffset(m)+(6*h);
    pf.Width =  Woffset(m);
    pf.ShortPinWidth = Woffset(m);
    S = sparameters(pf,freq);
    S11_W(m) = 20*log10(abs(S.Parameters(8)))
end
plot(Woffset,S11_W)
[S11_min_W,index_s11_W] = min(S11_W);
W_opt = Woffset(index_s11_W)
pf.Width = W_opt;
pf.GroundPlaneWidth = W_opt + (6*h);
pf.ShortPinWidth = W_opt;



