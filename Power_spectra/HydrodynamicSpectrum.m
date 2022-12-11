function spec = HydrodynamicSpectrum(f, para);

fc = para(1);
D = para(2);
fv = para(3);
fm = para(4);
f = abs(f);  %Kludge!
spec = D/pi^2*(1+sqrt(f/fv))./((fc - f.*sqrt(f./fv) - (f.^2)/fm).^2 + (f + f.*sqrt(f./fv)).^2);