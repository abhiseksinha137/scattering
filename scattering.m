Reset;
L=100e-9;
x=[-L:0.1:L]*1e-9;
y=x;

Rf=L;                   % Detector  circle at radius Rf
xfList=linspace(-L,L,100);  %  x points of the detector
yfList=Rf^2-xfList.^2;      %  corresponding y points of the detector
% xf=0;
% yf=L;

[X,Y]=meshgrid(x,y);
R=sqrt(X.^2+y.^2);

m=9.1093837e-31;
hbar=1.05457182e-34;
eV=1.6e-19;


E=linspace(0,20,200)*eV;   
E01=3*eV;    % Electron peak 1
E02=9*eV;    % Electron peak 2
E03=15*eV;   % Electron peak 3
sigma=3*eV;  % Width of each electron peak; 
P0=sqrt(3*exp(-(E-E01).^2/sigma^2)+  2*exp(-(E-E02).^2/sigma^2) + 1*exp(-(E-E03).^2/sigma^2)); % Initial energy spectrum
SC=sqrt(1./E); % Scattering Cross sectionis assumed to be inversely proportional to the energy of the electron;

P=zeros(size(E));  % Initilize the probabilities
for j=1:length(xfList)
    xf=xfList(j); % A point on the detector;
    yf=yfList(j);
    for i=1:length(E)
        k=sqrt(2*m*E(i)/hbar); % Wave vector of the electron of the corresponding energy
        % The probabiity Amplitude is the sum of the scattered wave + Initial wave 
        A=SC(i)*P0(i)*(1./R.*exp(1i*k*R+ 1i*k*sqrt((xf-X).^2+(yf-Y).^2))) + P0(i)*1/L*exp(1i*k*Rf);
    
        W(i)= sum(sum(A));
    end
    P=P+abs(W).^2;
    P(isinf(P))=0;
end
plot(E/eV, P0/max(P0), 'r', 'Linewidth', 2); hold on;
plot(E/eV,P/max(P), 'k','Linewidth', 2)
legend('Spectrum without Scattering', 'Spectrum with Scattering')
xlabel('Energy eV')
exportgraphics(gcf, 'Spectrum.png', Resolution=400);


