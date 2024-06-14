function E=fita(initialparameter,freq,xm,ym)
% Used in: fit_sine_wave.m
gs=initialparameter(1)*sin(2*pi*freq*xm+initialparameter(2))+initialparameter(3);
E=sum(abs(gs-ym).^2);
ya=initialparameter(1)*sin(2*pi*freq*xm+initialparameter(2))+initialparameter(3);
hold off
plot(xm,ym,'.-','linewidth',1);
hold on;
plot(xm,ya,'r');
axis([min(xm) max(xm) min(ym) max(ym)]);
drawnow;
end