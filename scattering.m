Reset;
L=100e-9;
x=[-L:0.1:L]*1e-9;
y=x;
xf=0;
yf=L;

[X,Y]=meshgrid(x,y);
R=sqrt(X.^2+y.^2);

m=9.1093837e-31;
hbar=1.05457182e-34;
eV=1.6e-19;


E=linspace(0,10,200)*eV;
E0=3*eV; sigma=5*eV;
P0=exp(-(E-E0).^2/sigma^2);
SC=1./E;
for i=1:length(E)
    k=sqrt(2*m*E(i)/hbar);
    A=SC(i)*P0(i)*(1./R.*exp(1i*k*R+ 1i*k*sqrt((xf-X).^2+(yf-Y).^2))) + 1/L*exp(1i*k*L);
    W(i)=sum(sum(A));
    P=abs(W);
    P(isinf(P))=0;
end
plot(E, P0/max(P0), 'r', 'Linewidth', 2); hold on;
plot(E,P/max(P), 'k','Linewidth', 2)
legend('Spectrum without Scattering', 'Spectrum with Scattering')
exportgraphics(gcf, 'Spectrum.png', Resolution=400);
